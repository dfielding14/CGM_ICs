"""Cooling/flow boundary tests against an independently chosen sonic state."""

import numpy as np
import pytest

import HaloPotential_new as Halo
import WiersmaCooling as Cool
import solve_ode as CF


@pytest.fixture(scope='module')
def timescale_problem():
    potential = Halo.PowerLaw(0, 100, 100)
    cooling = Cool.Constant_Cooling(1e-22)
    conditions = CF._calculate_sonic_point_conditions(.3, 1., potential, cooling, -1.)
    reference = CF._integrate_from_sonic_point(
        conditions, potential, cooling, 1., 10., .8, 1e-4,
        False, True, 2e4, .03, False, 'RK45', check_unbound=False)
    ratio = reference.cooling_time_Gyr[-1] / reference.flow_time_Gyr[-1]
    return potential, cooling, reference, ratio


def solve(potential, cooling, ratio, outer_radius=10., **kwargs):
    return CF.shoot_from_R_sonic_outer_boundary(
        potential, cooling, R_sonic_kpc=1., R_boundary_kpc=10.,
        boundary_value=ratio, R_max_kpc=outer_radius, R_min_kpc=.8,
        epsilon=1e-4, max_step=.03, boundary='timescale', **kwargs)


def test_recovers_known_flow_and_timescale_ratio(timescale_problem):
    potential, cooling, reference, ratio = timescale_problem
    solution = solve(potential, cooling, ratio)
    assert solution.cooling_time_Gyr[-1] > 0.
    assert solution.cooling_time_Gyr[-1] / solution.flow_time_Gyr[-1] == pytest.approx(
        ratio, rel=2e-5)
    assert solution.mass_flow_rate_Msun_per_yr == pytest.approx(
        reference.mass_flow_rate_Msun_per_yr, rel=2e-5)


def test_extent_and_potential_zero_do_not_change_selection(timescale_problem):
    _, cooling, reference, ratio = timescale_problem
    potential = Halo.PowerLaw(0, 100, 100, R_phi0_kpc=.001)
    solution = solve(potential, cooling, ratio, outer_radius=20.)
    assert solution.radius_kpc[-1] == pytest.approx(20.)
    assert solution.mass_flow_rate_Msun_per_yr == pytest.approx(
        reference.mass_flow_rate_Msun_per_yr, rel=2e-5)
    # The entropy equation predicts d ln(P/rho^gamma)/d ln r = t_flow/t_cool.
    offset = 1e-4
    states = solution.result.sol(np.log(10.) + np.array([-offset, offset]))
    ln_entropy = states[0] - (CF.GC.GAMMA - 1.) * states[1]
    slope = np.diff(ln_entropy)[0] / (2 * offset)
    assert slope == pytest.approx(1. / ratio, rel=2e-4)
    solution.sample_profile()
    assert solution.radius_kpc[-1] == pytest.approx(20.)
    assert np.all(np.isfinite(solution.temperature_K))


@pytest.mark.parametrize('ratio', [0., -1., np.nan, np.inf])
def test_invalid_ratio_is_rejected(timescale_problem, ratio):
    potential, cooling, _, _ = timescale_problem
    with pytest.raises(ValueError):
        solve(potential, cooling, ratio)


def test_unknown_boundary_is_rejected(timescale_problem):
    potential, cooling, _, _ = timescale_problem
    with pytest.raises(ValueError):
        CF.shoot_from_R_sonic_outer_boundary(
            potential, cooling, 1., 10., 1., 20., .8, boundary='unknown')


@pytest.mark.parametrize('cooling_time', [-1., 0., np.inf])
def test_non_cooling_boundary_cannot_be_a_match(timescale_problem, monkeypatch, cooling_time):
    potential, cooling, _, _ = timescale_problem
    # Isolate the boundary selector: the trial ODE remains a valid cooling flow.
    monkeypatch.setattr(CF.CGMSolution, 'cooling_time_Gyr', property(
        lambda self: np.full_like(self.radius_kpc, cooling_time)))
    with pytest.raises(RuntimeError, match='bracket'):
        solve(potential, cooling, 1.)


def test_two_narrow_hot_branch_matches_are_rejected(monkeypatch):
    from types import SimpleNamespace

    # Both crossings lie below the first nonzero node of the linear scan.
    def conditions(x, *args):
        return {'dlnM_dlnR': -1., 'x': x}

    def integrate(state, *args, **kwargs):
        x = state['x']
        ratio = .5 + 1e5 * (x - .0003) * (x - .0008)
        return SimpleNamespace(stop_reason=CF.StopReason.MAX_RADIUS,
                               radius_kpc=np.array([10.]),
                               cooling_time_Gyr=np.array([ratio]),
                               flow_time_Gyr=np.array([1.]))

    monkeypatch.setattr(CF, '_calculate_sonic_point_conditions', conditions)
    monkeypatch.setattr(CF, '_integrate_from_sonic_point', integrate)
    with pytest.raises(RuntimeError, match='Multiple .* brackets'):
        solve(Halo.PowerLaw(0, 100, 100), Cool.Constant_Cooling(1e-22), .5)
