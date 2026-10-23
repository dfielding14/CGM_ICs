import json

import numpy as np
import pytest
from scipy.integrate import quad, trapezoid

import generate_gotham_ics as IC
import unit_system as us


@pytest.fixture(scope="module")
def generated(tmp_path_factory):
    root = tmp_path_factory.mktemp("gotham")
    source = root / "athinput"
    source.write_text("""<units>
length_cgs = 3.0856775809623245e21
mass_cgs = 3.036951775493658e39
time_cgs = 3.15576e15
mu = 0.62
<potential>
r_scale = 10
rho_scale = 2
mass_gal = 6500
scale_gal = 1.5
z_gal = 0.15
r_200 = 150
rho_mean = 2.5e-5
<mesh>
nx1 = 64
nx2 = 64
nx3 = 64
x1min = -50
x1max = 50
x2min = -50
x2max = 50
x3min = -50
x3max = 50
nghost = 4
<hydro>
gamma = 1.66666667
T_max = 5e8 # retain case when writing Athena parameters
<problem>
r_circ = 10
metallicity = 0.3333333333333333
<dust>
d_z_init = 0.2
D_floor = 1e-7
<particles>
particle_type = star
""")
    config = IC.read_input(source)
    units = IC.CodeUnits(config)
    potential = IC.GalaxyPotential(config, units)
    destination = root / "run"
    solution = IC.generate(source, destination)
    return source, destination, config, units, potential, solution


def test_export_roundtrip_and_no_overwrite(generated):
    source, destination, config, units, potential, solution = generated
    output = IC.read_input(destination / "athinput.gotham")
    metadata = json.loads((destination / "ic_parameters.json").read_text())
    boundary = 1.5 * config.getfloat("potential", "r_200") * units.kpc
    target = units.mu * us.M_P_G * (potential.vc_kms(boundary) * us.KM_TO_CM)**2 / (IC.CF.GC.GAMMA * IC.KB_CGS)
    assert metadata["boundary_radius_kpc"] == pytest.approx(boundary)
    assert metadata["boundary_temperature_K"] == pytest.approx(target)
    assert metadata["f_T"] == 1.
    assert solution.radius_kpc[-1] >= boundary * (1 - 1e-12)
    temperature = np.exp(np.interp(np.log(boundary), np.log(solution.radius_kpc),
                                   np.log(solution.temperature_K)))
    assert temperature == pytest.approx(target, rel=2e-5)
    assert output.getfloat("hydro", "T_max") == 5e8
    expected_cooling = IC.Cool.Wiersma_Cooling((1 - 1e-7) / 3, 0)
    assert solution.cooling.lambda_cgs(1e6, 1e-4) == pytest.approx(
        expected_cooling.lambda_cgs(1e6, 1e-4), rel=1e-12, abs=0)
    assert output.getfloat("problem", "v_circ") == pytest.approx(potential.equatorial_vc_kms(10 * units.kpc) / units.kms)
    for name in ("cgm_profiles", "disk_profiles"):
        table = np.loadtxt(destination / f"tables/{name}.txt")
        assert table.shape[1] == 4
        assert np.all(np.isfinite(table))
        assert np.all(np.diff(table[:, 0]) > 0)
        assert np.all(table[:, 1:3] > 0)
    table = np.loadtxt(destination / "tables/cgm_profiles.txt")
    mdot = us.mass_flow_rate_from_density_velocity_radius(
        table[:, 2] * units.density, table[:, 3] * units.kms, table[:, 0] * units.kpc)
    np.testing.assert_allclose(mdot, solution.mass_flow_rate_Msun_per_yr, rtol=1e-12)
    np.testing.assert_allclose(table[:, 1] * units.temperature, solution.temperature_K, rtol=1e-14)
    with pytest.raises(FileExistsError):
        IC.generate(source, destination)


def test_monopole_force_and_potential_agree(generated):
    potential = generated[4]
    radius = np.geomspace(0.001, 400, 80)
    step = 1e-4
    vc2 = (potential.phi_kms2(radius * np.exp(step))
           - potential.phi_kms2(radius * np.exp(-step))) / (2 * step)
    np.testing.assert_allclose(potential.vc_kms(radius)**2, vc2, rtol=2e-7)
    for r in (0.01, 0.15, 0.7, 10, 100, 400):
        expected, _ = quad(lambda mu: potential.phi_cylindrical_kms2(r * np.sqrt(1 - mu*mu), r * mu),
                           0, 1, epsabs=1e-8, epsrel=1e-10)
        assert potential.phi_kms2(r) == pytest.approx(expected, rel=1e-7)


def test_exporter_passes_configurable_outer_temperature(generated, tmp_path, monkeypatch):
    source, _, _, units, potential, _ = generated
    boundary, f_T = 100., .5
    target = f_T * units.mu * us.M_P_G * (potential.vc_kms(boundary) * us.KM_TO_CM)**2 / (IC.CF.GC.GAMMA * IC.KB_CGS)

    def inspect_boundary(potential_arg, cooling, rs, rb, tb, rmax, rmin, **kwargs):
        assert rs == 1.
        assert rb == boundary
        assert tb == pytest.approx(target)
        assert rmax >= boundary
        raise RuntimeError("checked boundary parameters")

    monkeypatch.setattr(IC.CF, "shoot_from_R_sonic_outer_temperature", inspect_boundary)
    with pytest.raises(RuntimeError, match="checked boundary parameters"):
        IC.generate(source, tmp_path / "run", sonic_radius_kpc=1.,
                    f_T=f_T, boundary_radius_kpc=boundary)
    assert not (tmp_path / "run").exists()


def test_disk_has_requested_mass_in_full_potential(generated):
    _, destination, config, units, potential, _ = generated
    table = np.loadtxt(destination / "tables/disk_profiles.txt")
    radius = table[:, 0] * units.kpc
    sigma = []
    for r, temp, rho in zip(radius, table[:, 1], table[:, 2]):
        cs2 = temp * units.kms**2
        phi0 = potential.phi_cylindrical_kms2(r, 0)
        column, _ = quad(lambda z: np.exp(-(potential.phi_cylindrical_kms2(r, z) - phi0) / cs2),
                         0, 50 * units.kpc, points=[0.15, 1.5], epsabs=1e-10)
        sigma.append(2 * column * rho * units.density * us.KPC_TO_CM**3 / us.MSUN_TO_G)
    mass = trapezoid(2 * np.pi * radius * sigma, radius)
    assert mass == pytest.approx(2.5e9, rel=1e-4)


def test_clusters_have_correct_units_and_circular_speeds(generated):
    _, destination, _, units, potential, _ = generated
    stars = np.loadtxt(destination / "tables/clusters.txt")
    assert stars.shape == (1500, 8)
    radius = np.hypot(stars[:, 0], stars[:, 1]) * units.kpc
    assert np.max(radius) <= 0.5
    assert np.max(abs(stars[:, 2]) * units.kpc) <= 0.01
    assert np.max(stars[:, 6] * units.time / us.YR_TO_S) <= 1e7
    assert np.all((stars[:, 7] >= 1e4) & (stars[:, 7] <= 1e6))
    speed2 = (stars[:, 3]**2 + stars[:, 4]**2) * units.kms**2
    step = 1e-4
    force = (potential.phi_cylindrical_kms2(radius * np.exp(step), 0)
             - potential.phi_cylindrical_kms2(radius * np.exp(-step), 0)) / (2 * step)
    np.testing.assert_allclose(speed2, force, rtol=1e-7)
    np.testing.assert_allclose(stars[:, 0] * stars[:, 3] + stars[:, 1] * stars[:, 4], 0, atol=1e-15)


def test_flow_satisfies_momentum_and_energy_equations(generated):
    solution = generated[5]
    logr = np.log(solution.radius_kpc)
    v = solution.velocity_kms
    cs2 = solution.sound_speed_squared_kms2
    dlnv = np.gradient(np.log(v), logr, edge_order=2)
    dlnP = np.gradient(np.log(solution.pressure_dyne_cm2), logr, edge_order=2)
    vc2 = solution.potential.vc_kms(solution.radius_kpc)**2
    momentum_residual = v*v*dlnv + cs2 / IC.CF.GC.GAMMA * dlnP + vc2
    energy_derivative = np.gradient(solution.bernoulli_kms2, logr, edge_order=2)
    energy_rhs = cs2 / (IC.CF.GC.GAMMA * IC.CF.GC.GAMMA_M1) * solution.flow_time_Gyr / solution.cooling_time_Gyr
    keep = (abs(solution.radius_kpc / 0.7 - 1) > 0.02)
    keep[:2] = keep[-2:] = False
    assert np.max(abs(momentum_residual[keep]) / vc2[keep]) < 2e-3
    # Differentiating the exported samples loses accuracy at cooling-table slope changes.
    assert np.max(abs(energy_derivative[keep] - energy_rhs[keep]) / vc2[keep]) < 1e-2


def test_radial_resolution_convergence(generated):
    coarse = generated[5]
    metadata = json.loads((generated[1] / "ic_parameters.json").read_text())
    fine = IC.CF.shoot_from_R_sonic_outer_temperature(
        coarse.potential, coarse.cooling, 0.7, metadata["boundary_radius_kpc"],
        metadata["boundary_temperature_K"], metadata["outer_radius_kpc"], 0.1,
        max_step=0.001, epsilon=1e-4)
    assert fine.mass_flow_rate_Msun_per_yr == pytest.approx(coarse.mass_flow_rate_Msun_per_yr, rel=1e-4)
    for name in ("temperature_K", "density_cgs", "velocity_kms"):
        reference = np.interp(np.log(coarse.radius_kpc), np.log(fine.radius_kpc), np.log(getattr(fine, name)))
        assert np.max(abs(np.expm1(reference - np.log(getattr(coarse, name))))) < 2e-4
