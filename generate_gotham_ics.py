"""Generate a GOTHAM run directory from its athinput; all output tables use code units."""

from argparse import ArgumentParser
from configparser import ConfigParser
from pathlib import Path
import json
import re
import sys

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicHermiteSpline

sys.path.insert(0, str(Path(__file__).resolve().parent / "pysrc"))
import HaloPotential_new as Halo
import WiersmaCooling as Cool
import solve_ode as CF
import unit_system as us

# Match src/units/units.hpp, including its historical constants.
G_CGS = 6.67408e-8
KB_CGS = 1.3806488e-16


def read_input(path):
    text = Path(path).read_text().split("<par_end>", 1)[0]
    text = re.sub(r"^\s*<([^>]+)>\s*(?:#.*)?$", r"[\1]", text, flags=re.MULTILINE)
    config = ConfigParser(interpolation=None, inline_comment_prefixes=("#",), strict=False)
    config.optionxform = str
    config.read_string(text)
    return config


class CodeUnits:
    def __init__(self, config):
        self.length = config.getfloat("units", "length_cgs")
        self.mass = config.getfloat("units", "mass_cgs")
        self.time = config.getfloat("units", "time_cgs")
        self.mu = config.getfloat("units", "mu")
        if not np.all(np.isfinite([self.length, self.mass, self.time, self.mu])) or min(
                self.length, self.mass, self.time, self.mu) <= 0:
            raise ValueError("code units must be finite and positive")
        if not np.isclose(self.mu, CF.GC.MU, rtol=0, atol=1e-12):
            raise ValueError("the cooling-flow solver assumes mu=0.62")
        self.kpc = self.length / us.KPC_TO_CM
        self.kms = self.length / self.time / us.KM_TO_CM
        self.density = self.mass / self.length**3
        self.temperature = (self.length / self.time)**2 * self.mu * us.M_P_G / KB_CGS


class GalaxyPotential:
    """GOTHAM's full cylindrical potential; the 1D flow uses its angular mean."""

    def __init__(self, config, units):
        p = {name: config.getfloat("potential", name) for name in
             ("r_scale", "rho_scale", "mass_gal", "scale_gal", "z_gal", "r_200", "rho_mean")}
        if not np.all(np.isfinite(list(p.values()))) or min(p.values()) <= 0:
            raise ValueError("potential parameters must be finite and positive")
        self.rs = p["r_scale"] * units.kpc
        self.r_floor = 1e-10 * units.kpc
        rho_s = p["rho_scale"] * units.density * us.KPC_TO_CM**3 / us.MSUN_TO_G
        self.nfw_mass = 4 * np.pi * rho_s * self.rs**3
        self.disk = Halo.MiyamotoNagaiPotential(p["mass_gal"] * units.mass / us.MSUN_TO_G,
                                              p["scale_gal"] * units.kpc,
                                              p["z_gal"] * units.kpc)
        self.outer = Halo.OuterHaloPotential(
            p["rho_mean"] * units.density * us.KPC_TO_CM**3 / us.MSUN_TO_G,
            p["r_200"] * units.kpc)
        self.G = G_CGS * us.MSUN_TO_G / us.KPC_TO_CM / us.KM_TO_CM**2
        self.g_ratio = self.G / us.G_KPC_KM2_S2_PER_MSUN
        nodes, weights = np.polynomial.legendre.leggauss(64)
        self.cos_theta, self.weights = (nodes + 1) / 2, weights / 2

    def phi_cylindrical_kms2(self, R_kpc, z_kpc):
        radius = np.maximum(np.hypot(R_kpc, z_kpc), self.r_floor)
        x = radius / self.rs
        return (-self.G * self.nfw_mass / self.rs * np.log1p(x) / x
                + self.g_ratio * (self.disk.phi_kms2(R_kpc, z_kpc) + self.outer.phi_kms2(radius)))

    def _spherical_vc2(self, radius):
        x = np.asarray(radius) / self.rs
        return (self.G * self.nfw_mass / radius * (np.log1p(x) - x / (1 + x))
                + self.g_ratio * self.outer.vc_kms(radius)**2)

    def equatorial_vc_kms(self, radius):
        return np.sqrt(self._spherical_vc2(radius) + self.g_ratio * self.disk.vc_kms(radius)**2)

    def phi_kms2(self, radius_kpc):
        radius = np.asarray(radius_kpc)[..., None]
        return np.sum(self.weights * self.phi_cylindrical_kms2(
            radius * np.sqrt(1 - self.cos_theta**2), radius * self.cos_theta), axis=-1)

    def vc_kms(self, radius_kpc):
        radius = np.asarray(radius_kpc)[..., None]
        a, b = self.disk.a_kpc, self.disk.b_kpc
        height = np.sqrt((radius * self.cos_theta)**2 + b**2)
        distance2 = radius**2 + a*a + b*b + 2*a*height
        disk_vc2 = self.G * self.disk.M_Msun * radius**2 * (
            1 + a * self.cos_theta**2 / height) / distance2**1.5
        return np.sqrt(self._spherical_vc2(np.asarray(radius_kpc))
                       + np.sum(self.weights * disk_vc2, axis=-1))

    def dln_vc_dln_r(self, radius_kpc):
        step = 1e-4
        return (np.log(self.vc_kms(np.asarray(radius_kpc) * np.exp(step)))
                - np.log(self.vc_kms(np.asarray(radius_kpc) * np.exp(-step)))) / (2 * step)


class TabulatedPotential:
    """One Hermite interpolant for the spherical potential and its derivatives."""

    def __init__(self, potential, r_min_kpc, r_max_kpc, points=2049):
        if not 0 < r_min_kpc < r_max_kpc or points < 4:
            raise ValueError("require 0 < r_min < r_max and at least four table points")
        # Padding permits endpoint roundoff without extrapolating the physical domain.
        logr = np.linspace(np.log(r_min_kpc) - 1e-6, np.log(r_max_kpc) + 1e-6, points)
        radius = np.exp(logr)
        self._phi = CubicHermiteSpline(logr, potential.phi_kms2(radius),
                                      potential.vc_kms(radius)**2, extrapolate=False)

    def phi_kms2(self, radius_kpc):
        return self._phi(np.log(radius_kpc))

    def vc_kms(self, radius_kpc):
        return np.sqrt(self._phi(np.log(radius_kpc), 1))

    def dln_vc_dln_r(self, radius_kpc):
        logr = np.log(radius_kpc)
        return 0.5 * self._phi(logr, 2) / self._phi(logr, 1)


def outer_mesh_radius(config, units):
    ghost = config.getint("mesh", "nghost")
    outer_coordinates = []
    for axis in (1, 2, 3):
        lo, hi = (config.getfloat("mesh", f"x{axis}{side}") for side in ("min", "max"))
        if not lo < 0 < hi:
            raise ValueError("the GOTHAM IC mesh must straddle the galaxy center")
        nx = config.getint("mesh", f"nx{axis}")
        if nx <= 0 or ghost < 0:
            raise ValueError("mesh cell counts must be positive and nghost nonnegative")
        outer_coordinates.append(max(abs(lo), abs(hi)) + (ghost + 0.5) * (hi - lo) / nx)
    return np.linalg.norm(outer_coordinates) * units.kpc


def cooling_from_input(config):
    fluid = "mhd" if config.has_section("mhd") else "hydro"
    if not np.isclose(config.getfloat(fluid, "gamma"), CF.GC.GAMMA, rtol=1e-8, atol=0):
        raise ValueError("the cooling-flow solver assumes gamma=5/3")
    metallicity = config.getfloat("problem", "metallicity", fallback=1/3)
    # GOTHAM gives the CGM D_floor; d_z_init is only applied to disk material.
    dust_fraction = config.getfloat("dust", "D_floor", fallback=1e-20)
    if not 0 <= dust_fraction <= 1:
        raise ValueError("D_floor must be between zero and one")
    return Cool.Wiersma_Cooling(metallicity * (1 - dust_fraction), 0.0)


def disk_profile(potential, units, zmin_kpc, zmax_kpc, mass_Msun=2.5e9, temperature_K=1e4):
    radius = np.geomspace(1e-3, 20.0, 512)
    scale = potential.disk.a_kpc
    # Normalize the truncated exponential disk to the requested mass in the mesh's z range.
    fraction = 1 - (1 + radius[-1] / scale) * np.exp(-radius[-1] / scale)
    sigma0 = mass_Msun / (2 * np.pi * scale**2 * fraction)
    cs2 = KB_CGS * temperature_K / (units.mu * us.M_P_G) / us.KM_TO_CM**2
    density = []
    for r in radius:
        phi0 = potential.phi_cylindrical_kms2(r, 0.0)
        def vertical(z):
            return np.exp(-(potential.phi_cylindrical_kms2(r, z) - phi0) / cs2)
        points = [z for z in (-10 * potential.disk.b_kpc, -potential.disk.b_kpc, 0,
                              potential.disk.b_kpc, 10 * potential.disk.b_kpc)
                  if zmin_kpc < z < zmax_kpc]
        thickness, _ = quad(vertical, zmin_kpc, zmax_kpc, points=points, epsabs=1e-10)
        density.append(sigma0 * np.exp(-r / scale) / thickness)
    density = np.asarray(density) * us.MSUN_TO_G / us.KPC_TO_CM**3
    return np.column_stack((radius / units.kpc, np.full_like(radius, temperature_K / units.temperature),
                            density / units.density, np.zeros_like(radius)))


def generate(athinput, output_dir, sonic_radius_kpc=0.7, seed=42,
             f_T=1.0, boundary_radius_kpc=None, cooling_to_flow_ratio=None, no_disk=False):
    from ClusterGenerator import ClusterParameters, StarClusterGenerator

    output_dir = Path(output_dir)
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    config = read_input(athinput)
    initialize_disk = config.getboolean("problem", "initialize_disk", fallback=True) and not no_disk
    units = CodeUnits(config)
    potential = GalaxyPotential(config, units)
    if boundary_radius_kpc is None:
        boundary_radius_kpc = 1.5 * config.getfloat("potential", "r_200") * units.kpc
    if not np.isfinite(f_T) or f_T <= 0:
        raise ValueError("f_T must be finite and positive")
    if not np.isfinite(boundary_radius_kpc) or not 0.1 < sonic_radius_kpc < boundary_radius_kpc:
        raise ValueError("require 0.1 kpc < sonic radius < boundary radius")
    if cooling_to_flow_ratio is not None:
        if not np.isfinite(cooling_to_flow_ratio) or cooling_to_flow_ratio <= 0:
            raise ValueError("cooling_to_flow_ratio must be finite and positive")
        if f_T != 1.0:
            raise ValueError("choose either f_T or cooling_to_flow_ratio")
    boundary_temperature_K = float(f_T * units.mu * us.M_P_G
                                   * (potential.vc_kms(boundary_radius_kpc) * us.KM_TO_CM)**2
                                   / (CF.GC.GAMMA * KB_CGS))
    cooling = cooling_from_input(config)
    metallicity = config.getfloat("problem", "metallicity", fallback=1/3)
    dust_fraction = config.getfloat("dust", "D_floor", fallback=1e-20)
    mesh_radius = outer_mesh_radius(config, units)
    if initialize_disk:
        for axis in (1, 2):
            lo, hi = (config.getfloat("mesh", f"x{axis}{side}") for side in ("min", "max"))
            if min(-lo, hi) * units.kpc < 20:
                raise ValueError("the mesh must contain the 20 kpc disk in both horizontal directions")
    outer_radius = max(mesh_radius, boundary_radius_kpc)
    flow_potential = TabulatedPotential(potential, 0.1, outer_radius)
    boundary = "temperature" if cooling_to_flow_ratio is None else "timescale"
    boundary_value = boundary_temperature_K if boundary == "temperature" else cooling_to_flow_ratio
    # Boundary shooting needs tighter tolerances than legacy bound/unbound selection.
    solution = CF.shoot_from_R_sonic_outer_boundary(
        flow_potential, cooling, sonic_radius_kpc, boundary_radius_kpc, boundary_value,
        outer_radius, 0.1, max_step=np.inf, epsilon=1e-4, atol=1e-10, rtol=1e-10,
        boundary=boundary, verbose=True)
    if solution is None or not np.isclose(solution.radius_kpc[-1], outer_radius, rtol=1e-8):
        raise RuntimeError("no cooling flow reaches the outer mesh, including ghost cells")
    if not np.isclose(solution.radius_kpc[0], 0.1, rtol=1e-8):
        raise RuntimeError("the inward solution did not reach 0.1 kpc")
    solution.sample_profile()
    if (np.any(solution.nH_cm3 < 1e-8) or np.any(solution.temperature_K < 100)
            or np.any(solution.temperature_K > 9.1612998e8)):
        raise ValueError("the flow extends beyond GOTHAM's cooling-table domain")
    cgm = np.column_stack((solution.radius_kpc / units.kpc, solution.temperature_K / units.temperature,
                           solution.density_cgs / units.density, solution.velocity_kms / units.kms))
    disk = None
    if initialize_disk:
        disk = disk_profile(potential, units, config.getfloat("mesh", "x3min") * units.kpc,
                            config.getfloat("mesh", "x3max") * units.kpc)
    generator = StarClusterGenerator(ClusterParameters(n_clusters=1500, radius_kpc=0.5))
    clusters = generator.generate_clusters(seed=seed)
    generator.set_velocity_profile(lambda r_pc: potential.equatorial_vc_kms(r_pc / 1000))
    # Preserve the historical 10 Myr formation interval and approximately 10 km/s vertical dispersion.
    stars = np.column_stack((clusters[["x_pc", "y_pc", "z_pc"]].to_numpy() / (1000 * units.kpc),
                             clusters[["vx", "vy"]].to_numpy() / units.kms,
                             clusters.vz.to_numpy() * 10 / units.kms,
                             clusters.t_min.to_numpy() * 1e8 * us.YR_TO_S / units.time,
                             clusters.mass_msun.to_numpy()))
    for data in ((cgm, disk, stars) if initialize_disk else (cgm, stars)):
        if not np.all(np.isfinite(data)):
            raise ValueError("nonfinite values in generated ICs")
    if np.any(cgm[:, 1:3] <= 0) or (initialize_disk and np.any(disk[:, 1:3] <= 0)):
        raise ValueError("profile temperature and density must be positive")

    config.set("problem", "profile_file", "tables/cgm_profiles.txt")
    config.set("problem", "initialize_disk", str(initialize_disk).lower())
    if initialize_disk:
        config.set("problem", "disk_profile_file", "tables/disk_profiles.txt")
    else:
        config.remove_option("problem", "disk_profile_file")
    rcirc = config.getfloat("problem", "r_circ") * units.kpc
    if not np.isfinite(rcirc) or rcirc < 0:
        raise ValueError("r_circ must be finite and nonnegative")
    v_circ = float(potential.equatorial_vc_kms(rcirc) / units.kms) if rcirc else 0.0
    config.set("problem", "v_circ", str(v_circ))
    config.set("particles", "star_particle_file", "tables/clusters.txt")
    # exclusive creation prevents overwriting a run while numerical work is in progress
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "tables").mkdir()
    np.savetxt(output_dir / "tables/cgm_profiles.txt", cgm, header="r T rho v_inflow [code units]")
    if initialize_disk:
        np.savetxt(output_dir / "tables/disk_profiles.txt", disk, header="R T rho_midplane v_unused [code units]")
    np.savetxt(output_dir / "tables/clusters.txt", stars,
               header="x y z vx vy vz t_creation [code units]; mass [Msun]")
    with (output_dir / "athinput.gotham").open("w") as handle:
        for section in config.sections():
            handle.write(f"<{section}>\n")
            for name, value in config.items(section):
                handle.write(f"{name} = {value}\n")
            handle.write("\n")
    metadata = dict(source_input=str(Path(athinput).resolve()), sonic_radius_kpc=sonic_radius_kpc,
                    mass_flow_Msun_per_yr=solution.mass_flow_rate_Msun_per_yr,
                    outer_radius_kpc=outer_radius, gas_metallicity_solar=metallicity * (1 - dust_fraction),
                    mesh_radius_kpc=mesh_radius, boundary_radius_kpc=boundary_radius_kpc,
                    outer_boundary=boundary,
                    seed=seed, cgm_gravity="angular mean of the GOTHAM potential",
                    cooling="optically thin z=0 Wiersma; no local heating or shielding",
                    disk_mass_Msun=2.5e9 if initialize_disk else 0, disk_temperature_K=1e4,
                    disk_vertical_interval_kpc=[config.getfloat("mesh", "x3min") * units.kpc,
                                                config.getfloat("mesh", "x3max") * units.kpc])
    if boundary == "temperature":
        metadata.update(boundary_temperature_K=boundary_temperature_K, f_T=f_T)
    else:
        metadata["cooling_to_flow_ratio"] = cooling_to_flow_ratio
    (output_dir / "ic_parameters.json").write_text(json.dumps(metadata, indent=2) + "\n")
    return solution


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("athinput", type=Path)
    parser.add_argument("output_dir", type=Path, help="new run directory; existing directories are never overwritten")
    parser.add_argument("--r-sonic", type=float, default=0.7, help="sonic radius in kpc (default: 0.7)")
    parser.add_argument("--r-boundary", type=float,
                        help="matching radius in kpc (default: 1.5 times input r_200)")
    boundary_options = parser.add_mutually_exclusive_group()
    boundary_options.add_argument("--f-t", type=float, default=1.0,
                        help="outer c_s^2/v_c^2, with adiabatic sound speed (default: 1)")
    boundary_options.add_argument("--cooling-to-flow", type=float, nargs="?", const=1.0,
                                  help="outer net t_cool/t_flow (default when selected: 1)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--no-disk", action="store_true", help="omit the initial gas disk; retain stellar gravity")
    args = parser.parse_args()
    solution = generate(args.athinput, args.output_dir, args.r_sonic, args.seed,
                        f_T=args.f_t, boundary_radius_kpc=args.r_boundary,
                        cooling_to_flow_ratio=args.cooling_to_flow, no_disk=args.no_disk)
    print(f"Wrote {args.output_dir}: Mdot={solution.mass_flow_rate_Msun_per_yr:.8g} Msun/yr")


if __name__ == "__main__":
    main()
