# CGM library and fiducial performance checks

2026-10-23. These checks use the current GOTHAM potential's spherical average,
the configured gas-phase metallicity, and signed redshift-zero Wiersma cooling.
The fixed stellar potential remains present; the library creates no initial gas
disk or star clusters. The two paper fiducials are $r_{\rm sonic}=0.7$ and 35 kpc.
The full 0.01--300 kpc survey and forbidden-zone bisection have not been run.

## Changes

- GOTHAM profile generation uses adaptive RK45 without a compulsory radial step
  cap. Both absolute and relative tolerances are $10^{-9}$ for the logarithmic
  dependent variables. The generic solver's defaults remain unchanged.
- A 2049-node cubic Hermite spline in $\ln r$ stores the angularly averaged
  potential and its analytic radial derivative. Circular speed and its slope
  come from that same spline, so the force and potential remain consistent.
  The original 64-point angular quadrature is used only during setup. The full
  cylindrical potential remains available for disk and cluster generation.
- Dense ODE output samples profiles at $\Delta\ln r\leq0.002$, retaining the
  adaptive knots around sharp transitions. Output resolution no longer controls
  the integration step size.
- One cooling table and potential table serve an entire survey invocation.
  Spawned workers receive these objects through their initializer; they do not
  reload cooling files for each radius.
- Profiles are cached as non-pickled NPZ files. Cache fingerprints include the
  input contents, numerical settings, relevant source files, cooling arrays,
  and NumPy/SciPy versions. A changed fingerprint requires a new output directory.
- Independent points can use a process pool. Exact cache hits bypass the solver
  and do not start workers. Failed searches are recorded and retried on later
  invocations, rather than cached as evidence of physical nonexistence.

The existing complete-run exporter also uses the tabulated spherical potential
and adaptive integration. Its historical gas-disk and cluster defaults remain
unchanged; use the new library command for CGM-only output.

## Warm starts and the outer condition

`shoot_from_R_sonic(..., x_guess=...)` accepts a previous shooting parameter.
A failed guess updates the original full search bracket, rather than restricting
the search to an assumed nearby interval. The library's `--warm-start` option
uses the nearest completed profile's parameter in subsequent batches.

**Nearby-radius warm starts are opt-in.** The current outer prescription accepts
the first trial that reaches the requested outer radius while satisfying the
Bernoulli filter. It can accept a finite interval of shooting parameters. For
the small fiducial, changing a successful $x$ by just 0.01% changed the outer
temperature by about 5%, even though the mass flux changed by only about 0.06%.
Consequently, a warm start is not guaranteed to reproduce the original search's
selected profile. The default retains the original bisection sequence and uses
only exact cached profiles. Reusing the same accepted $x$ at the same radius
also speeds up convergence checks without changing which flow is being tested.

The Bernoulli acceptance rule is unchanged. Its potential-zero dependence and
the need for a physically specified outer condition remain documented in
`TODO.md`. No rejected or failed shoot is automatically labeled a forbidden
physical radius.

## Fiducial results

Both fiducials reach the paper-sized domain, 0.1--244.8 kpc. The reference used
direct angular quadrature, a maximum logarithmic step of 0.002, and the previous
$10^{-8}$ integration tolerances. The optimized run used the settings above.
Both used the same $10^{-4}$ sonic offset and original shooting sequence.
These are local wall-clock measurements, not universal performance guarantees.

| Sonic radius | Mass flux [$M_\odot/{\rm yr}$] | Reference solve | Optimized solve | Largest profile difference |
| --- | --- | --- | --- | --- |
| 0.7 kpc | 0.14599430 | 5.21 s | 0.46 s | 0.011% |
| 35 kpc | 1.00753458 | 2.11 s | 0.21 s | 0.00088% |

The profile differences cover temperature, density, and inward speed, compared
on a common logarithmic radius coordinate. The selected shooting parameters
are identical, and mass fluxes agree to better than $3\times10^{-9}$ relative.
Serial and two-worker outputs are identical. The two-point serial job took
about 0.71 s including setup and saving, versus 0.85 s with two workers. Startup
overhead outweighs parallelism for this small batch. Reloading both cached
profiles took about 0.008 s within the already-started Python process.

Additional checks:

- Over 0.001--400 kpc, tabulated circular speeds agree with the direct angular
  quadrature to $4.5\times10^{-9}$ relative; the absolute error in
  $d\ln v_c/d\ln r$ is below $7.4\times10^{-6}$.
- Increasing the table to 4097 nodes or tightening integration tolerances to
  $10^{-10}$ changes the accepted profiles by less than 0.011%. Reducing the
  sonic offset to $10^{-5}$ changes them by less than 0.023%.
- These refinement checks hold the accepted $x$ fixed. A cold search at
  $10^{-10}$ encountered a failed integration on a rejected, nearly sonic
  trial in the large-radius case; this is a numerical failure, not evidence
  against the accepted profile.
- Mass-flux conservation is within $2.3\times10^{-16}$ relative. Finite-difference
  momentum and energy residuals, normalized by $v_c^2$, are below $3.4\times10^{-4}$
  and $2.5\times10^{-3}$ respectively, excluding the immediate sonic offset and
  endpoint stencils. Cooling-table slope changes limit differentiated samples.
- All 38 tests pass, including the existing notebook execution and complete
  GOTHAM export tests, direct-versus-tabulated profiles, cache invalidation,
  warm-start fallback, multiprocessing, and incomplete-profile rejection.

**Outer coverage remains a separate limitation.** The current GOTHAM mesh needs
a profile extending to 354.529 kpc. The 0.7 kpc fiducial reaches it. The 35 kpc
fiducial does not produce an accepted shoot under the present prescription,
including a reference retry with a $10^{-8}$ shooting-bracket tolerance.
The 244.8 kpc profiles are explicitly marked `covers_mesh: false`; they must not
be treated as full-domain ICs for that mesh. No three-dimensional simulation or
long relaxation was run for these performance changes.

## Output and use

```bash
python generate_cgm_library.py /path/to/athinput.gotham /path/to/library \
    --r-sonic 0.7 35 --r-max 244.8 --workers 2
```

Each accepted radius produces `rsonic_<radius>.npz`, containing physical-unit
radius, temperature, density, inward speed, Mach number, cooling time, flow time,
and JSON metadata. Cooling time can be negative in net-heating gas. When the
profile lies within the supported cooling-table range, a matching `.txt` file
contains `r T rho v_inflow` in GOTHAM code units. The temperature column is
$P/\rho$ in those units. `parameters.json` records the model and numerical
settings, and `summary.json` records the most recent requested batch's outcomes.

No profile file is written for a numerical failure, unsuccessful shoot, or
incomplete radial coverage. `covers_mesh` and `within_cooling_table` are separate
metadata flags. The command does not create a runnable athinput or configure
the simulation's gas-disk switch. The full IC library still requires the outer
condition decision before its excluded intervals can be interpreted physically.
