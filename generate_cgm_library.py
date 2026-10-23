"""Generate cached, spherical CGM profiles; no gas disk or star clusters."""

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from pathlib import Path
import json
import multiprocessing
import os
import time

import numpy as np

import generate_gotham_ics as IC


def _initialize_worker(potential, cooling, settings):
    global _potential, _cooling, _settings
    _potential, _cooling, _settings = potential, cooling, settings


def _solve_point(task):
    radius, x_guess = task
    start = time.perf_counter()
    entry = dict(sonic_radius_kpc=radius, x_guess=x_guess, cache_hit=False)
    try:
        solution = IC.CF.shoot_from_R_sonic(
            _potential, _cooling, radius, _settings["r_max_kpc"], _settings["r_min_kpc"],
            max_step=np.inf, epsilon=1e-4, atol=1e-9, rtol=1e-9, x_guess=x_guess)
    except (RuntimeError, ValueError) as error:
        entry.update(status="numerical_failure", detail=str(error))
        return entry, None
    finally:
        entry["solve_seconds"] = time.perf_counter() - start
    if solution is None:
        entry.update(status="no_accepted_shoot")
        return entry, None
    if not (np.isclose(solution.radius_kpc[0], _settings["r_min_kpc"], rtol=1e-8)
            and np.isclose(solution.radius_kpc[-1], _settings["r_max_kpc"], rtol=1e-8)):
        entry.update(status="incomplete_radial_coverage")
        return entry, None
    solution.sample_profile()
    arrays = {name: getattr(solution, name) for name in
              ("radius_kpc", "temperature_K", "density_cgs", "velocity_kms",
               "mach_number", "cooling_time_Gyr", "flow_time_Gyr")}
    if any(not np.all(np.isfinite(arrays[name])) or np.any(arrays[name] <= 0)
           for name in ("radius_kpc", "temperature_K", "density_cgs", "velocity_kms")):
        entry.update(status="invalid_profile")
        return entry, None
    entry.update(status="accepted", shooting_x=solution.shooting_x,
                 shooting_trials=solution.shooting_trials,
                 mass_flow_Msun_per_yr=solution.mass_flow_rate_Msun_per_yr,
                 covers_mesh=bool(_settings["r_max_kpc"] >= _settings["mesh_outer_radius_kpc"]),
                 within_cooling_table=bool(np.all(solution.nH_cm3 >= 1e-8)
                     and np.all(solution.temperature_K >= 100)
                     and np.all(solution.temperature_K <= 9.1612998e8)))
    return entry, arrays


def run_survey(athinput, output_dir, sonic_radii=(0.7, 35.0), r_min_kpc=0.1,
               r_max_kpc=None, workers=1, warm_start=False):
    """Reuse exact profiles by default; nearby shooting guesses are explicitly opt-in.

    A warm guess can select a different accepted flow under the current outer
    condition. Cache identities include the model, numerical settings and source.
    """
    started = time.perf_counter()
    config = IC.read_input(athinput)
    units = IC.CodeUnits(config)
    mesh_outer = IC.outer_mesh_radius(config, units)
    r_max_kpc = mesh_outer if r_max_kpc is None else float(r_max_kpc)
    radii = sorted(set(float(r) for r in sonic_radii))
    if (not radii or not np.all(np.isfinite([r_min_kpc, r_max_kpc, *radii]))
            or not 0 < r_min_kpc < min(radii) * (1 - 1e-4)
            or not max(radii) * (1 + 1e-4) < r_max_kpc or workers < 1):
        raise ValueError("require positive radii, r_min < every sonic point < r_max, and workers >= 1")
    potential = IC.TabulatedPotential(IC.GalaxyPotential(config, units), r_min_kpc, r_max_kpc)
    cooling = IC.cooling_from_input(config)
    settings = dict(r_min_kpc=r_min_kpc, r_max_kpc=r_max_kpc,
                    mesh_outer_radius_kpc=float(mesh_outer), warm_start=warm_start,
                    epsilon=1e-4, max_step="adaptive", output_max_log_step=0.002,
                    atol=1e-9, rtol=1e-9, potential_points=2049, terminal_unbound=True,
                    numpy_version=np.__version__, scipy_version=IC.CF.scipy.__version__)
    digest = sha256(Path(athinput).read_bytes() + json.dumps(settings, sort_keys=True).encode())
    root = Path(__file__).resolve().parent
    for name in ("generate_cgm_library.py", "generate_gotham_ics.py", "pysrc/solve_ode.py",
                 "pysrc/WiersmaCooling.py", "pysrc/unit_system.py", "pysrc/HaloPotential_new.py"):
        digest.update((root / name).read_bytes())
    for data in (cooling._rates, cooling._logT, cooling._logn):
        digest.update(data.tobytes())
    fingerprint = digest.hexdigest()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "parameters.json"
    if manifest_path.exists():
        if json.loads(manifest_path.read_text())["fingerprint"] != fingerprint:
            raise ValueError("cached model, solver or settings changed; use a new output directory")
    else:
        manifest_path.write_text(json.dumps(dict(fingerprint=fingerprint, settings=settings,
            source_input=str(Path(athinput).resolve()), input_text=Path(athinput).read_text(),
            gas_disk=False, gravity="spherical average of GOTHAM", cooling="z=0 signed Wiersma"),
            indent=2) + "\n")

    cached = {}
    for path in sorted(output_dir.glob("rsonic_*.npz")):
        with np.load(path, allow_pickle=False) as data:
            entry = json.loads(str(data["metadata"]))
            if entry["fingerprint"] != fingerprint:
                raise ValueError(f"profile does not match this model: {path}")
            entry["profile"] = str(path.resolve())
            cached[entry["sonic_radius_kpc"]] = entry
    entries = {r: dict(cached[r], cache_hit=True) for r in radii if r in cached}
    pending = [r for r in radii if r not in entries]
    _initialize_worker(potential, cooling, settings)
    pool = (ProcessPoolExecutor(max_workers=min(workers, len(pending)),
            mp_context=multiprocessing.get_context("spawn"), initializer=_initialize_worker,
            initargs=(potential, cooling, settings)) if workers > 1 and pending else None)
    try:
        # Each round can reuse completed shoots; points within a round are independent.
        for start in range(0, len(pending), workers):
            batch = pending[start:start + workers]
            hints = [cached[min(cached, key=lambda old: abs(np.log(old / r)))]["shooting_x"]
                     if warm_start and cached else None for r in batch]
            tasks = list(zip(batch, hints))
            results = pool.map(_solve_point, tasks) if pool else map(_solve_point, tasks)
            for entry, arrays in results:
                radius = entry["sonic_radius_kpc"]
                entry["fingerprint"] = fingerprint
                if arrays is not None:
                    stem = f"rsonic_{radius!r}"
                    entry["profile"] = str((output_dir / (stem + ".npz")).resolve())
                    temporary = output_dir / f"{stem}.{os.getpid()}.tmp.npz"
                    np.savez_compressed(temporary, metadata=json.dumps(entry), **arrays)
                    temporary.replace(entry["profile"])
                    if entry["within_cooling_table"]:
                        table = np.column_stack((arrays["radius_kpc"] / units.kpc,
                            arrays["temperature_K"] / units.temperature,
                            arrays["density_cgs"] / units.density, arrays["velocity_kms"] / units.kms))
                        np.savetxt(output_dir / (stem + ".txt"), table,
                                   header="r T rho v_inflow [GOTHAM code units]; CGM only")
                    cached[radius] = entry
                entries[radius] = entry
    finally:
        if pool:
            pool.shutdown()
    summary = dict(elapsed_seconds=time.perf_counter() - started, workers=workers,
                   results=[entries[r] for r in radii])
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("athinput", type=Path)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--r-sonic", nargs="+", type=float, default=[0.7, 35.0])
    parser.add_argument("--r-min", type=float, default=0.1)
    parser.add_argument("--r-max", type=float, help="default: outer mesh corner including ghost cells")
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--warm-start", action="store_true",
                        help="reuse nearby shooting guesses; can change which accepted flow is selected")
    args = parser.parse_args()
    result = run_survey(args.athinput, args.output_dir, args.r_sonic, args.r_min,
                        args.r_max, args.workers, args.warm_start)
    for entry in result["results"]:
        print(f"r_sonic={entry['sonic_radius_kpc']:g} kpc: {entry['status']}, "
              f"Mdot={entry.get('mass_flow_Msun_per_yr', float('nan')):.8g}, "
              f"cached={entry['cache_hit']}")
    print(f"Elapsed: {result['elapsed_seconds']:.3f} s")


if __name__ == "__main__":
    main()
