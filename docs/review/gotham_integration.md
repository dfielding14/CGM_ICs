# GOTHAM integration review

Review date: 2026-10-23. Generator baseline: `04ddec89`; AthenaK baseline:
`179a7038` on `gotham_dev`. The new command is `generate_gotham_ics.py`.
It produces a complete run bundle from the actual GOTHAM athinput, without editing
the source input, historical tables, or notebooks. The current example notebooks
did not reproduce the former production export.

## Corrections

1. **Retain net heating.** The Wiersma reader clipped negative net rates before
   taking their logarithm. Photoheating was therefore erased, with spurious large
   logarithmic slopes around the clipping boundary. Interpolation now preserves
   the signed net rate, bilinearly in log temperature and log hydrogen density,
   matching GOTHAM. Slopes differentiate that same interpolant and are defined
   only where net cooling is positive.
2. **Match the composition and gas metallicity.** Use the $Y=0.248$ H/He table
   shipped in GOTHAM. The previous helium prescription selected $Y=0.238$ for
   subsolar metallicity. For the CGM, gas metallicity is
   `problem/metallicity * (1 - dust/D_floor)`; `dust/d_z_init` applies to the disk.
   High-density lookup clamps to the final table row. The exporter rejects flows
   outside GOTHAM's temperature range or below its minimum lookup density.
3. **Use the simulation's gravity and units.** Read the NFW, stellar
   Miyamoto-Nagai, and outer-halo parameters directly from the input. Match the
   constants in AthenaK's unit system. Cluster circular speeds and the CGM
   rotation normalization now use this potential instead of a spherical proxy
   or a separately hard-coded speed.
4. **Normalize the disk in its actual vertical potential.** The former export
   normalized using only stellar gravity over a restricted vertical interval.
   GOTHAM reconstructs the disk in the full potential. The new normalization uses
   that same potential over the mesh's vertical extent.
5. **Distinguish integration failure from reaching the boundary.** Failed ODE
   integration raises an error; nonterminal unbound crossings no longer override
   a later successful boundary arrival. The exporter requires both radial ends
   and enough outer coverage for corner ghost cells. A failed shoot cannot
   quietly produce a truncated run bundle.
6. **Correct rotating-flow bookkeeping.** The generic rotating solver conserves
   $j=v_c(R_c)R_c$, so its centrifugal term is $j^2/r^2$, and its Bernoulli
   parameter includes $j^2/(2r^2)$. The old expression implicitly assumed a flat
   circular-speed curve. Mass flux at the offset shooting radius now uses that
   actual radius. Unbound trial temperatures reduce the bracket instead of
   aborting the search. The GOTHAM export itself uses the nonrotating radial flow.
7. **Improve root solving and diagnostics.** Solve the sonic density in log space
   so the absolute root tolerance does not swamp dilute densities. Gas mass is
   the integral over the resolved interval, with zero mass assigned below its
   inner edge. Cooling-radius diagnostics select the outermost bracketed positive
   cooling-time crossing and return NaN if none exists. Cooling-driven rotation
   and turbulence diagnostics are marked undefined in net-heating gas.

## Model and export conventions

The one-dimensional CGM uses the angular mean of the full cylindrical potential,

$$
\bar\Phi(r)=\int_0^1\Phi(r\sqrt{1-\mu^2},r\mu)\,d\mu,
\qquad v_c^2=r\,\frac{d\bar\Phi}{dr}.
$$

The angular quadrature has 64 nodes. The radial force is differentiated
analytically before angular averaging. The default sonic radius is 0.7 kpc;
the inner radial endpoint is 0.1 kpc. The outer temperature condition is

$$
T(R_b)=f_T\frac{\mu m_p}{\gamma k_B}v_c^2(R_b),
\qquad R_b=1.5R_{200}\ \text{by default}.
$$

Both $R_b$ (`--r-boundary`, kpc) and $f_T$ (`--f-t`, default 1) are parameters.
The default $f_T=1$ equates adiabatic sound speed and circular speed at $R_b$.
The sonic temperature is adjusted to match this boundary; density and mass
inflow rate are outputs. The endpoint includes both the physical matching radius
and corner ghost cells, whichever is farther out. The solver requires
$0.1\,{\rm kpc}<R_{\rm sonic}<R_b$ and permits sonic points beyond the mesh.

An alternative `--cooling-to-flow [q_b]` imposes $t_{\rm cool}/t_{\rm flow}=q_b$
at the same $R_b$, with $q_b=1$ when the argument is omitted. The cooling time
uses the net cooling rate from the flow equations, and $t_{\rm flow}=r/|v_r|$.
This sets the steady entropy slope $d\ln(P/\rho^\gamma)/d\ln r=1/q_b$ and
lets the boundary temperature emerge. Heating-dominated boundary states are
rejected. `--f-t` and `--cooling-to-flow` are mutually exclusive; metadata
records the selected condition and target.

These boundaries replace the exporter's former bound/unbound trial selection.
The outer-halo potential is confining and has no finite zero at infinity, so
the sign of its Bernoulli parameter cannot define escape. Both modes
therefore disable that rejection and are invariant under additive potential
offsets. The older `shoot_from_R_sonic` retains its behavior for notebook users.
The former baseline export's approximately $0.146\,M_\odot\,{\rm yr}^{-1}$ was
selected by that older criterion and is not a prediction of the new boundary.

The disk retains $M_d=2.5\times10^9\,M_\odot$, $T_d=10^4$ K, and a 20 kpc
outer cutoff. Its surface density is exponential with the input's stellar
radial scale. At each radius,

$$
\rho(R,z)=\rho_0(R)\exp\left[-\frac{\Phi(R,z)-\Phi(R,0)}{c_T^2}\right],
\qquad c_T^2=\frac{k_B T_d}{\mu m_p},
$$

and the vertical integral sets $\rho_0$ to recover that surface density.
The analytic mass normalization includes the radial cutoff. GOTHAM computes
pressure-supported disk rotation when reading the profile. This continuum
normalization does not guarantee the same mass on an underresolved mesh.

Use `--no-disk` to omit this initial gas component. The exporter writes
`problem/initialize_disk = false` and omits the disk table and its file path;
GOTHAM skips both reading and adding the gas disk. The same setting in an input
file is preserved on regeneration. Stellar-disk gravity, CGM rotation and
turbulence, cluster generation, and feedback timing are retained. Existing
inputs without this setting continue to initialize a gas disk.

The historical cluster prescription is retained: 1500 clusters inside 0.5 kpc,
uniform in disk area, heights within 10 pc, masses between $10^4$ and $10^6$
$M_\odot$ drawn from a slope-two power law, formation times within 10 Myr,
and approximately 10 km/s vertical velocity dispersion. Planar speeds are
equatorial circular speeds in the full potential. The random seed is recorded.

| File | Columns |
| --- | --- |
| `tables/cgm_profiles.txt` | radius, $P/\rho$, density, positive inward speed; all code units |
| `tables/disk_profiles.txt` | cylindrical radius, $P/\rho$, midplane density, unused zero velocity; all code units |
| `tables/clusters.txt` | x, y, z, vx, vy, vz, formation time in code units; mass in solar masses |

The profile temperature column is **isothermal** $P/\rho$, without a factor of
$\gamma$. GOTHAM applies the inward sign to CGM radial velocity. The exporter
preserves case-sensitive input names, including `T_max`, and updates the three
file paths and `v_circ`. It does not change the feedback schedule: the reviewed
template has `time/tlim=10` and `SN/delay=10`, intentionally ending before the
first supernova.

## Validation

- Python tests cover execution of all example notebooks, raw-HDF5
  cooling comparisons, signed heating, dilute sonic roots, rotating forces,
  solver failure handling, disk mass, cluster units, and exported mass flux.
- Outer-temperature regressions recover a known transonic flow, check invariance
  under potential offsets, separate matching radius from profile extent, and
  reject a matched branch that encounters a second sonic point during extension.
- Independent momentum and energy residuals are checked along an exported flow.
  Halving its maximum log-radius step changes temperature, density, and velocity
  by less than $2\times10^{-4}$ in the regression case.
- AthenaK's `tst/scripts/hydro/gotham_ics.py` compiles the actual C++ potential and
  cooling implementation. Across 2048 samples the maximum relative potential
  discrepancy is $6.1\times10^{-15}$. The maximum cooling discrepancy divided by
  $|\Lambda|+10^{-23}\,{\rm erg\,cm^3\,s^{-1}}$ is $7.0\times10^{-6}$; the small
  difference comes from rounded log axes in GOTHAM's historical table.
- A serial CPU GOTHAM build loads the generated tables and 1500 clusters and
  advances two steps on a $32^3$ mesh. Output fields remain finite, with positive
  density and internal energy. This startup check disables AMR and reduces the
  turbulence mode count; it is not a production convergence test.

## Remaining approximations

- **Three-dimensional equilibrium:** the CGM is a spherical cooling flow in
  averaged gravity. GOTHAM adds rotation and turbulence afterwards. The inner
  rotating flow and disk interface will relax; this is not an exact steady 3D
  solution. The table is clamped below 0.1 kpc.
- **Thermal equilibrium:** the CGM solver uses optically thin, redshift-zero
  photoionization-equilibrium cooling. GOTHAM additionally uses local heating and
  cell-width shielding. The imposed $10^4$ K disk is in vertical hydrostatic
  balance, not guaranteed thermal equilibrium under those terms.
- **Composition:** constant $\mu=0.62$, hydrogen fraction 0.75, and scaled solar
  metals remain approximations. This does not add a nonequilibrium ionization or
  variable-composition model.
- **Disk support and gravity:** no gas self-gravity, turbulent pressure, or
  magnetic pressure is included in the vertical construction. The gas disk is
  25% of the template's stellar mass, so self-gravity is a material omission
  shared with the simulation, not an IC-only inconsistency that this fixes.
- **Clusters:** the spatial, mass, formation-time, and vertical-kick distributions
  remain prescribed; they are not a self-consistent stellar distribution function.
- **Scope:** historical exported tables, pickles, and unused alternative potential
  models were not recalibrated. GPU, multi-rank MPI, production AMR, feedback, and
  long-time relaxation were not validated in this review.
