"""Compare outer temperature and cooling/flow-time conditions in GOTHAM gravity.

From the IC repository:
    python compare_outer_boundaries.py --output-dir ../../build/outer-boundaries
"""

from argparse import ArgumentParser
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.integrate import OdeSolution
from scipy.optimize import brentq

import generate_gotham_ics as gotham


def ambiguous_timescale_flows(potential, cooling, rs, rb, rmax, target):
    """Expose every detected bracket for comparison; never choose a production branch."""
    def flow(x, inward=False):
        conditions = gotham.CF._calculate_sonic_point_conditions(x, rs, potential, cooling, None)
        if conditions["dlnM_dlnR"] >= 0:
            raise ValueError("wrong sonic branch")
        result = gotham.CF._integrate_from_sonic_point(
            conditions, potential, cooling, rs, rb, 0.1, 1e-4, False, inward,
            2e4, np.inf, False, "RK45", atol=1e-10, rtol=1e-10, check_unbound=False)
        if result.stop_reason != gotham.CF.StopReason.MAX_RADIUS:
            raise RuntimeError("trial did not reach matching radius")
        return result
    def residual(x):
        trial = flow(x)
        q = trial.cooling_time_Gyr[-1] / trial.flow_time_Gyr[-1]
        if not np.isfinite(q) or q <= 0:
            raise ValueError("boundary requires positive finite net cooling time")
        return float(np.log(q / target))
    grid = np.unique(np.concatenate((np.linspace(1e-5, 1 - 1e-5, 33),
                                    np.geomspace(1e-5, 1 - 1e-5, 33))))
    values = []
    for x in grid:
        try:
            with np.errstate(over="raise", divide="raise", invalid="raise"):
                values.append(residual(x))
        except (gotham.CF.NoValidDensityError, gotham.CF.NoTranssonicSolutionError,
                ValueError, RuntimeError, FloatingPointError):
            values.append(np.nan)
    roots = [brentq(residual, a, b, xtol=1e-12, rtol=1e-12)
             for a, b, fa, fb in zip(grid[:-1], grid[1:], values[:-1], values[1:])
             if np.isfinite(fa) and np.isfinite(fb) and fa * fb < 0]
    assert len(roots) > 1, "ambiguous production target must yield multiple diagnostic roots"
    results = []
    for index, x in enumerate(roots):
        assert abs(residual(x)) < 1e-5
        solution = flow(x, True)
        row = dict(branch=f"diagnostic {index + 1}", sonic_x=x,
                   status="diagnostic ambiguous branch; full coverage", extension_error="")
        if rmax > rb:
            try:
                extension = gotham.CF.IntegrateFlowEquations(
                    solution.mass_flow_rate_Msun_per_yr, solution.temperature_K[-1],
                    solution.density_cgs[-1], potential, cooling, 1, rb, rmax,
                    max_step=np.inf, check_unbound=False, atol=1e-10, rtol=1e-10)
                if extension.stop_reason != gotham.CF.StopReason.MAX_RADIUS:
                    raise RuntimeError(f"extension stops at {extension.radius_kpc[-1]:g} kpc: "
                                       f"{extension.stop_reason.value}")
                solution.result.t = np.concatenate((solution.result.t, extension.result.t[1:]))
                solution.result.y = np.concatenate((solution.result.y, extension.result.y[:, 1:]), axis=1)
                solution.result.sol = OdeSolution(
                    np.concatenate((solution.result.sol.ts, extension.result.sol.ts[1:])),
                    solution.result.sol.interpolants + extension.result.sol.interpolants)
            except (ValueError, RuntimeError) as exc:
                row.update(status="diagnostic ambiguous branch; boundary only", extension_error=str(exc))
        results.append((solution, row))
    return results


def compare(athinput, output_dir, boundary_radius_kpc=None):
    config = gotham.read_input(athinput)
    units = gotham.CodeUnits(config)
    exact_potential = gotham.GalaxyPotential(config, units)
    cooling = gotham.cooling_from_input(config)
    r200 = config.getfloat("potential", "r_200") * units.kpc
    rb = 1.5 * r200 if boundary_radius_kpc is None else boundary_radius_kpc
    concentration = r200 / exact_potential.rs
    halo_mass = exact_potential.nfw_mass * (np.log1p(concentration) - concentration / (1 + concentration))
    if not np.isfinite(rb) or rb < r200 or rb <= 20:
        raise ValueError("outer boundary must be finite, at or beyond R200, and beyond 20 kpc")
    mesh_radius = gotham.outer_mesh_radius(config, units)
    rmax = max(mesh_radius, rb)
    potential = gotham.TabulatedPotential(exact_potential, 0.1, rmax)
    reference_temperature = float(units.mu * gotham.us.M_P_G
        * (exact_potential.vc_kms(rb) * gotham.us.KM_TO_CM)**2
        / (gotham.CF.GC.GAMMA * gotham.KB_CGS))
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "source_athinput.gotham").write_text(Path(athinput).read_text())
    metadata = dict(source_input=str(Path(athinput).resolve()), r200_kpc=r200,
                    boundary_radius_kpc=rb, mesh_radius_kpc=mesh_radius,
                    reference_temperature_K=reference_temperature,
                    NFW_mass_within_R200_Msun=halo_mass,
                    temperature_definition="cs^2(Rb) = f_T * vc^2(Rb); adiabatic cs",
                    timescale_definition="q = net cooling time / (r / abs(v)) at Rb",
                    gas_mass_definition="integral from 0.1 kpc to R200; no initial gas disk",
                    atol=1e-10, rtol=1e-10, max_step="adaptive", epsilon=1e-4,
                    cooling_temperature_range_K=np.exp(cooling._logT[[0, -1]]).tolist(),
                    cooling_nH_range_cm3=np.exp(cooling._logn[[0, -1]]).tolist())
    (output_dir / "parameters.json").write_text(json.dumps(metadata, indent=2) + "\n")
    figure, axes = plt.subplots(4, 2, figsize=(11, 12), sharex=True, sharey="row")
    colors = {0.5: "#0072B2", 1.0: "#009E73", 2.0: "#D55E00"}
    rows = []
    for column, rs in enumerate((1.0, 20.0)):
        for mode in ("temperature", "timescale"):
            for factor in (0.5, 1.0, 2.0):
                target = factor * reference_temperature if mode == "temperature" else factor
                row = dict(sonic_radius_kpc=rs, boundary=mode, boundary_factor=factor,
                           boundary_target=target, boundary_radius_kpc=rb,
                           requested_outer_radius_kpc=rmax, status="", extension_error="")
                kwargs = dict(potential=potential, cooling=cooling, R_sonic_kpc=rs,
                              R_boundary_kpc=rb, boundary_value=target, boundary=mode,
                              R_max_kpc=rmax, R_min_kpc=0.1, max_step=np.inf,
                              epsilon=1e-4, atol=1e-10, rtol=1e-10)
                solutions = []
                try:
                    solution = gotham.CF.shoot_from_R_sonic_outer_boundary(**kwargs)
                    solutions = [(solution, dict(status="full coverage"))]
                except (ValueError, RuntimeError) as exc:
                    row["extension_error"] = str(exc)
                    if mode == "timescale" and factor == 0.5 and "Multiple" in str(exc):
                        solutions = ambiguous_timescale_flows(potential, cooling, rs, rb, rmax, factor)
                    else:
                        kwargs["R_max_kpc"] = rb
                        try:
                            solution = gotham.CF.shoot_from_R_sonic_outer_boundary(**kwargs)
                            solutions = [(solution, dict(status="boundary matched; mesh extension failed"))]
                        except (ValueError, RuntimeError) as match_exc:
                            row["status"] = "boundary matching failed"
                            row["matching_error"] = str(match_exc)
                            rows.append(row)
                            print(json.dumps(row), flush=True)
                for solution, details in solutions:
                    row = dict(row, **details)
                    # Evaluate the actual boundary before joining/resampling can move its knot.
                    log_rb = np.log(rb)
                    boundary_state = solution.result.sol(log_rb)
                    boundary_temperature = float(np.exp(boundary_state[0]))
                    boundary_density = float(np.exp(boundary_state[1]))
                    boundary_nH = boundary_density * gotham.CF.GC.X / gotham.us.M_P_G
                    boundary_velocity = (solution.mass_flow_rate_Msun_per_yr
                        * gotham.us.MSUN_PER_YR_TO_G_PER_S
                        / (4 * np.pi * (rb * gotham.us.KPC_TO_CM)**2 * boundary_density)
                        / gotham.us.KM_TO_CM)
                    boundary_tflow = float(gotham.us.flow_time_Gyr(rb, boundary_velocity))
                    boundary_tcool = float(gotham.us.cooling_time_Gyr(
                        boundary_density,
                        gotham.us.sound_speed_squared_kms2(boundary_temperature, gotham.CF.GC.GAMMA,
                                                           gotham.CF.GC.MU),
                        boundary_nH, cooling.lambda_cgs(boundary_temperature, boundary_nH),
                        gotham.CF.GC.GAMMA, gotham.CF.GC.GAMMA_M1))
                    actual_q = boundary_tcool / boundary_tflow
                    actual = boundary_temperature if mode == "temperature" else actual_q
                    residual = actual / target - 1
                    assert np.isfinite(residual) and abs(residual) < 1e-5, (mode, rs, factor, residual)
                    solution.sample_profile()
                    radius = solution.radius_kpc
                    assert np.all(np.isfinite(radius)) and np.all(np.diff(radius) > 0)
                    assert radius[0] <= 0.100001 and radius[-1] >= rb * (1 - 1e-8)
                    if row["status"].endswith("full coverage"):
                        assert radius[-1] >= rmax * (1 - 1e-8)
                    temperature, nH = solution.temperature_K, solution.nH_cm3
                    tcool, tflow = solution.cooling_time_Gyr, solution.flow_time_Gyr
                    ratio = tcool / tflow
                    for values in (temperature, nH, solution.density_cgs, solution.mach_number, tflow):
                        assert np.all(np.isfinite(values)) and np.all(values > 0)
                    def at_radius(values, r):
                        assert radius[0] <= r <= radius[-1]
                        return float(np.interp(np.log(r), np.log(radius), values))
                    outside = ((temperature < np.exp(cooling._logT[0]))
                               | (temperature > np.exp(cooling._logT[-1]))
                               | (nH < np.exp(cooling._logn[0]))
                               | (nH > np.exp(cooling._logn[-1])))
                    row.update(matched_temperature_K=boundary_temperature, matched_timescale_ratio=actual_q,
                               fractional_boundary_residual=residual,
                               equivalent_f_T=boundary_temperature / reference_temperature,
                               cooling_time_at_boundary_Gyr=boundary_tcool,
                               flow_time_at_boundary_Gyr=boundary_tflow,
                               sonic_temperature_K=float(np.exp(at_radius(np.log(temperature), rs))),
                               mass_flow_rate_Msun_per_yr=solution.mass_flow_rate_Msun_per_yr,
                               gas_mass_within_R200_Msun=at_radius(solution.gas_mass_Msun, r200),
                               gas_to_NFW_mass_within_R200=at_radius(solution.gas_mass_Msun, r200) / halo_mass,
                               nH_at_R200_cm3=float(np.exp(at_radius(np.log(nH), r200))),
                               minimum_radius_kpc=float(radius[0]), maximum_radius_kpc=float(radius[-1]),
                               minimum_temperature_K=float(np.min(temperature)),
                               maximum_temperature_K=float(np.max(temperature)),
                               minimum_nH_cm3=float(np.min(nH)), maximum_nH_cm3=float(np.max(nH)),
                               outside_cooling_table=bool(np.any(outside)),
                               fraction_samples_outside_cooling_table=float(np.mean(outside)))
                    np.savez_compressed(output_dir / f"profile_Rs{rs:g}_{mode}{factor:g}{row.get('branch', '').replace(' ', '_')}.npz",
                                        radius_kpc=radius, temperature_K=temperature, nH_cm3=nH,
                                        density_cgs=solution.density_cgs, velocity_kms=solution.velocity_kms,
                                        mach_number=solution.mach_number, cooling_time_Gyr=tcool,
                                        flow_time_Gyr=tflow, cooling_to_flow_time=ratio,
                                        free_fall_time_Gyr=solution.free_fall_time_Gyr,
                                        gas_mass_Msun=solution.gas_mass_Msun)
                    for index, values in enumerate((temperature, nH, solution.mach_number,
                                                    np.where(ratio > 0, ratio, np.nan))):
                        if "branch" not in row:
                            axes[index, column].loglog(radius, values, color=colors[factor],
                                                      ls="--" if mode == "temperature" else "-", lw=1.5)
                    rows.append(row)
                    print(json.dumps(row), flush=True)
        axes[0, column].set_title(rf"$R_{{\rm sonic}}={rs:g}$ kpc")
        for ax in axes[:, column]:
            ax.axvline(r200, color="0.5", ls=":", lw=1)
            ax.axvline(rb, color="0.5", ls="-.", lw=1)
            ax.grid(alpha=0.18)
        for index in (2, 3):
            axes[index, column].axhline(1, color="0.5", lw=0.8)
        axes[-1, column].set_xlabel("Radius [kpc]")
    for ax, label in zip(axes[:, 0], ("Temperature [K]", r"$n_H$ [cm$^{-3}$]",
                                    "Mach number", r"$t_{\rm cool}/t_{\rm flow}$")):
        ax.set_ylabel(label)
    legend = [Line2D([], [], color=colors[f], label=rf"$f_T$ or $q={f:g}$") for f in colors]
    legend += [Line2D([], [], color="0.2", ls="--", label="Temperature BC"),
               Line2D([], [], color="0.2", ls="-", label="Timescale BC")]
    figure.legend(handles=legend, loc="upper center", bbox_to_anchor=(0.5, 0.955),
                  ncol=3, frameon=False)
    figure.suptitle(f"Outer boundary comparison: $R_b={rb / r200:g}R_{{200}}={rb:.2f}$ kpc\n"
                   r"Dotted: $R_{200}$; dash-dot: $R_b$", fontsize=12)
    figure.text(0.5, 0.005,
                "Ambiguous q=0.5 branches appear in the separate diagnostic figure. Incomplete profiles end at the matching radius.",
                ha="center", fontsize=9)
    figure.subplots_adjust(left=0.095, right=0.985, bottom=0.065, top=0.865, hspace=0.08, wspace=0.045)
    figure.savefig(output_dir / "outer_boundary_comparison.png", dpi=180)
    figure.savefig(output_dir / "outer_boundary_comparison.pdf")
    plt.close(figure)
    diagnostic, diagnostic_axes = plt.subplots(4, 2, figsize=(11, 12), sharex=True, sharey="row")
    for row in rows:
        if row["status"] == "boundary matching failed" or (row["boundary"] != "temperature" and "branch" not in row):
            continue
        rs, factor, mode = row["sonic_radius_kpc"], row["boundary_factor"], row["boundary"]
        branch = row.get("branch", "")
        with np.load(output_dir / f"profile_Rs{rs:g}_{mode}{factor:g}{branch.replace(' ', '_')}.npz") as data:
            column = 0 if rs == 1 else 1
            for index, name in enumerate(("temperature_K", "nH_cm3", "mach_number", "cooling_to_flow_time")):
                diagnostic_axes[index, column].loglog(data["radius_kpc"],
                    np.where(data[name] > 0, data[name], np.nan),
                    color=colors[factor] if mode == "temperature" else ("#9467bd" if branch.endswith("1") else "#e377c2"),
                    ls="--" if mode == "temperature" else "-", lw=1.5)
    for column, rs in enumerate((1, 20)):
        diagnostic_axes[0, column].set_title(rf"$R_{{\rm sonic}}={rs}$ kpc")
        for ax in diagnostic_axes[:, column]:
            ax.axvline(r200, color="0.5", ls=":", lw=1)
            ax.axvline(rb, color="0.5", ls="-.", lw=1)
            ax.grid(alpha=0.18)
        for index in (2, 3):
            diagnostic_axes[index, column].axhline(1, color="0.5", lw=0.8)
        diagnostic_axes[-1, column].set_xlabel("Radius [kpc]")
    for ax, label in zip(diagnostic_axes[:, 0], ("Temperature [K]", r"$n_H$ [cm$^{-3}$]",
                                              "Mach number", r"$t_{\rm cool}/t_{\rm flow}$")):
        ax.set_ylabel(label)
    diagnostic.legend(handles=[Line2D([], [], color=colors[f], ls="--", label=rf"Temperature $f_T={f:g}$")
                               for f in colors] +
                      [Line2D([], [], color="#9467bd", label=r"$q=0.5$ root 1 (hotter)"),
                       Line2D([], [], color="#e377c2", label=r"$q=0.5$ root 2 (cooler)")],
                      loc="upper center", bbox_to_anchor=(0.5, 0.955), ncol=3, frameon=False)
    diagnostic.suptitle(r"Diagnostic ambiguous $q=0.5$ branches alongside temperature BCs"
                        + f"\n$R_b={rb:.2f}$ kpc; these branches are not selected by the exporter", fontsize=12)
    diagnostic.text(0.5, 0.005, "Cooling rates clamp outside the source table: consult the domain flags in the CSV.",
                    ha="center", fontsize=9)
    diagnostic.subplots_adjust(left=0.095, right=0.985, bottom=0.065, top=0.865, hspace=0.08, wspace=0.045)
    diagnostic.savefig(output_dir / "outer_boundary_ambiguous.png", dpi=180)
    diagnostic.savefig(output_dir / "outer_boundary_ambiguous.pdf")
    plt.close(diagnostic)
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with (output_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--athinput", type=Path,
                        default=Path(__file__).resolve().parents[2] / "inputs/gotham/athinput.gotham")
    parser.add_argument("--output-dir", type=Path, required=True, help="new output directory")
    parser.add_argument("--r-boundary-kpc", type=float, help="default: 1.5 R200")
    args = parser.parse_args()
    compare(args.athinput, args.output_dir, args.r_boundary_kpc)
