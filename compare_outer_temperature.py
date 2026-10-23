"""Compare sonic radii and outer temperatures in GOTHAM's spherical mean gravity.

Example, from the IC repository:
    python compare_outer_temperature.py --output-dir ../../build/outer-temperature
"""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import generate_gotham_ics as gotham


def compare(athinput, output_dir, boundary_radius_kpc=None):
    config = gotham.read_input(athinput)
    units = gotham.CodeUnits(config)
    potential = gotham.GalaxyPotential(config, units)
    cooling = gotham.Cool.Wiersma_Cooling(
        config.getfloat("problem", "metallicity", fallback=1/3)
        * (1 - config.getfloat("dust", "D_floor", fallback=1e-20)), 0.0)
    r200 = config.getfloat("potential", "r_200") * units.kpc
    rb = 1.5 * r200 if boundary_radius_kpc is None else boundary_radius_kpc
    if not np.isfinite(rb) or rb < r200 or rb <= 20:
        raise ValueError("outer boundary must be finite, at or beyond R200, and beyond 20 kpc")
    ghost = config.getint("mesh", "nghost")
    coordinates = []
    for axis in (1, 2, 3):
        lo, hi = (config.getfloat("mesh", f"x{axis}{side}") for side in ("min", "max"))
        coordinates.append(max(abs(lo), abs(hi)) + (ghost + 0.5)
                           * (hi - lo) / config.getint("mesh", f"nx{axis}"))
    mesh_radius = float(np.linalg.norm(coordinates) * units.kpc)
    rmax = max(mesh_radius, rb)
    # Use the same temperature normalization as the production exporter.
    reference_temperature = float(units.mu * gotham.us.M_P_G
        * (potential.vc_kms(rb) * gotham.us.KM_TO_CM)**2 / (gotham.CF.GC.GAMMA * gotham.KB_CGS))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "source_athinput.gotham").write_text(Path(athinput).read_text())
    metadata = dict(source_input=str(Path(athinput).resolve()), r200_kpc=r200,
                    boundary_radius_kpc=rb, mesh_radius_kpc=mesh_radius,
                    reference_temperature_K=reference_temperature,
                    temperature_definition="cs^2(Rb) = f_T * vc^2(Rb); adiabatic cs",
                    gas_mass_definition="integral from 0.1 kpc to R200; no initial gas disk")
    (output_dir / "parameters.json").write_text(json.dumps(metadata, indent=2) + "\n")
    figure, axes = plt.subplots(4, 2, figsize=(10, 11), sharex=True, sharey="row")
    colors = {0.5: "#0072B2", 1.0: "#009E73", 2.0: "#D55E00"}
    rows = []
    for column, rs in enumerate((1.0, 20.0)):
        for factor in (0.5, 1.0, 2.0):
            target = factor * reference_temperature
            row = dict(sonic_radius_kpc=rs, f_T=factor, target_temperature_K=target,
                       boundary_radius_kpc=rb, requested_outer_radius_kpc=rmax,
                       status="", extension_error="")
            kwargs = dict(potential=potential, cooling=cooling, R_sonic_kpc=rs,
                          R_boundary_kpc=rb, T_boundary_K=target, R_max_kpc=rmax,
                          R_min_kpc=0.1, max_step=0.01, epsilon=1e-4)
            try:
                solution = gotham.CF.shoot_from_R_sonic_outer_temperature(**kwargs)
                row["status"] = "full coverage"
            except (ValueError, RuntimeError) as exc:
                row["extension_error"] = str(exc)
                kwargs["R_max_kpc"] = rb
                try:
                    solution = gotham.CF.shoot_from_R_sonic_outer_temperature(**kwargs)
                    row["status"] = "boundary matched; mesh extension failed"
                except (ValueError, RuntimeError) as match_exc:
                    row["status"] = "boundary matching failed"
                    row["matching_error"] = str(match_exc)
                    rows.append(row)
                    print(json.dumps(row), flush=True)
                    continue
            radius = solution.radius_kpc
            assert np.all(np.isfinite(radius)) and np.all(np.diff(radius) > 0)
            assert radius[0] <= 0.100001 and radius[-1] >= rb * (1 - 1e-8)
            log_radius = np.log(radius)
            def at_radius(values, r):
                return float(np.interp(np.log(r), log_radius, values))
            def positive_at_radius(values, r):
                return float(np.exp(at_radius(np.log(values), r)))
            temperature = solution.temperature_K
            tcool = solution.cooling_time_Gyr
            tff = solution.free_fall_time_Gyr
            ratio = tcool / tff
            actual = positive_at_radius(temperature, rb)
            row.update(matched_temperature_K=actual, fractional_temperature_residual=actual/target-1,
                       sonic_temperature_K=positive_at_radius(temperature, rs),
                       mass_flow_rate_Msun_per_yr=solution.mass_flow_rate_Msun_per_yr,
                       gas_mass_within_R200_Msun=at_radius(solution.gas_mass_Msun, r200),
                       nH_at_R200_cm3=positive_at_radius(solution.nH_cm3, r200),
                       temperature_at_R200_K=positive_at_radius(temperature, r200),
                       tcool_to_tff_at_R200=at_radius(ratio, r200),
                       minimum_radius_kpc=float(radius[0]), maximum_radius_kpc=float(radius[-1]))
            np.savez_compressed(output_dir / f"profile_Rs{rs:g}_fT{factor:g}.npz",
                                radius_kpc=radius, temperature_K=temperature,
                                nH_cm3=solution.nH_cm3, density_cgs=solution.density_cgs,
                                velocity_kms=solution.velocity_kms, mach_number=solution.mach_number,
                                cooling_time_Gyr=tcool, free_fall_time_Gyr=tff,
                                cooling_to_free_fall_time=ratio, gas_mass_Msun=solution.gas_mass_Msun)
            for index, values in enumerate((temperature, solution.nH_cm3,
                                             solution.mach_number, np.where(ratio > 0, ratio, np.nan))):
                axes[index, column].loglog(radius, values, color=colors[factor],
                                          label=rf"$f_T={factor:g}$")
            rows.append(row)
            print(json.dumps(row), flush=True)
        axes[0, column].set_title(rf"$R_{{\rm sonic}}={rs:g}$ kpc")
        for ax in axes[:, column]:
            ax.axvline(r200, color="0.5", ls=":", lw=1)
            ax.axvline(rb, color="0.5", ls="--", lw=1)
            ax.grid(alpha=0.18)
        axes[2, column].axhline(1, color="0.5", lw=0.8)
        axes[-1, column].set_xlabel("Radius [kpc]")
    for ax, label in zip(axes[:, 0], ("Temperature [K]", r"$n_H$ [cm$^{-3}$]",
                                    "Mach number", r"$t_{\rm cool}/t_{\rm ff}$")):
        ax.set_ylabel(label)
    axes[0, 0].legend(loc="best", frameon=False)
    figure.suptitle(r"Spherical CGM flow: $c_s^2(R_b)=f_T v_c^2(R_b)$"
                   + f"\n$R_b={rb:.2f}$ kpc; dotted: $R_{{200}}$; dashed: $R_b$", fontsize=12)
    if any(row["status"] == "boundary matched; mesh extension failed" for row in rows):
        figure.text(0.5, 0.005, "Profiles that fail mesh extension are shown only to the matching radius.",
                    ha="center", fontsize=9)
    figure.tight_layout(rect=(0, 0.025, 1, 0.95))
    figure.savefig(output_dir / "outer_temperature_comparison.png", dpi=180)
    figure.savefig(output_dir / "outer_temperature_comparison.pdf")
    plt.close(figure)
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    if any(row["status"] == "boundary matching failed" for row in rows):
        raise RuntimeError("some boundary targets failed; see summary.csv")
    return rows


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--athinput", type=Path,
                        default=Path(__file__).resolve().parents[2] / "inputs/gotham/athinput.gotham")
    parser.add_argument("--output-dir", type=Path, required=True, help="new output directory")
    parser.add_argument("--r-boundary-kpc", type=float, help="default: 1.5 R200")
    args = parser.parse_args()
    compare(args.athinput, args.output_dir, args.r_boundary_kpc)
