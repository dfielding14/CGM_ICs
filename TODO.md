# TODO

- [ ] Resolve the Bernoulli rejection criterion before interpreting limits on sonic radius.

  The default sonic-point shooting in [solve_ode.py](pysrc/solve_ode.py)
  (`terminal_unbound=True`) rejects an outward solution that starts with positive
  Bernoulli parameter or reaches zero:

  $$
  B = \frac{v_r^2}{2} + \frac{c_s^2}{\gamma-1} + \bar\Phi(r).
  $$

  Here $c_s$ is the adiabatic sound speed and the exported CGM flow is nonrotating.
  The [GOTHAM outer-halo potential](pysrc/HaloPotential_new.py)
  (`OuterHaloPotential.phi_kms2`) contains positive $\sqrt{r}$ and $r^2$ terms,
  so it has no finite zero at infinity. Consequently, $B=0$ is not a physical
  escape threshold for this potential. Adding a constant to $\Phi$ changes
  the accepted solutions without changing the gravitational acceleration or
  the flow equations. This can impose an artificial upper limit on sonic radius.

  Define the intended finite-radius outer condition and replace or justify this
  filter before a large-sonic-radius survey. The current shooter accepts the
  first trial that reaches `R_max`; it does not match a prescribed outer density,
  temperature, or pressure. Verify that the revised selection is invariant under
  an additive shift of potential (with any reference energy shifted consistently),
  then reassess which sonic radii are possible. The criterion remains unchanged
  until this item is resolved.
