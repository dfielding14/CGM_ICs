"""Physical boundary checks independent of the GOTHAM exporter."""

import numpy as np
import pytest

import HaloPotential_new as Halo
import WiersmaCooling as Cool
import solve_ode as CF


@pytest.fixture(scope='module')
def boundary_problem():
    potential = Halo.PowerLaw(0, 100, 100)
    cooling = Cool.Constant_Cooling(1e-22)
    conditions = CF._calculate_sonic_point_conditions(.3, 1., potential, cooling, -1.)
    reference = CF._integrate_from_sonic_point(
        conditions, potential, cooling, 1., 10., .8, 1e-4,
        False, True, 2e4, .03, False, 'RK45')
    return potential, cooling, reference


def solve(potential, cooling, temperature, outer_radius=10.):
    return CF.shoot_from_R_sonic_outer_temperature(
        potential, cooling, R_sonic_kpc=1., R_boundary_kpc=10.,
        T_boundary_K=temperature, R_max_kpc=outer_radius, R_min_kpc=.8,
        epsilon=1e-4, max_step=.03)


def test_recovers_known_flow_and_matches_temperature(boundary_problem):
    potential, cooling, reference = boundary_problem
    solution = solve(potential, cooling, reference.temperature_K[-1])
    assert solution.stop_reason == CF.StopReason.MAX_RADIUS
    assert solution.radius_kpc[0] == pytest.approx(.8)
    assert solution.radius_kpc[-1] == pytest.approx(10.)
    assert solution.temperature_K[-1] == pytest.approx(reference.temperature_K[-1], rel=2e-5)
    assert solution.mass_flow_rate_Msun_per_yr == pytest.approx(
        reference.mass_flow_rate_Msun_per_yr, rel=2e-5)
    radii = solution.radius_kpc
    assert np.all(solution.mach_number[radii < .999] > 1.)
    assert np.all(solution.mach_number[radii > 1.001] < 1.)


def test_boundary_is_independent_of_profile_extent(boundary_problem):
    potential, cooling, reference = boundary_problem
    solution = solve(potential, cooling, reference.temperature_K[-1], outer_radius=20.)
    assert solution.radius_kpc[-1] == pytest.approx(20.)
    actual = np.exp(np.interp(np.log(10.), np.log(solution.radius_kpc),
                              np.log(solution.temperature_K)))
    assert actual == pytest.approx(reference.temperature_K[-1], rel=2e-4)
    assert solution.mass_flow_rate_Msun_per_yr == pytest.approx(
        reference.mass_flow_rate_Msun_per_yr, rel=2e-5)


def test_resampling_preserves_extension_beyond_temperature_boundary(boundary_problem):
    potential, cooling, reference = boundary_problem
    solution = solve(potential, cooling, reference.temperature_K[-1], outer_radius=20.)
    knots = solution.result.t.copy()
    values = solution.result.y.copy()
    assert solution.result.sol.ts[-1] == pytest.approx(np.log(20.))
    assert np.any(knots > np.log(10.))
    solution.sample_profile()
    indices = np.searchsorted(solution.result.t, knots)
    np.testing.assert_array_equal(solution.result.t[indices], knots)
    np.testing.assert_allclose(solution.result.y[:, indices], values, rtol=0, atol=1e-10)
    assert solution.radius_kpc[-1] == pytest.approx(20.)


def test_potential_zero_cannot_select_the_flow(boundary_problem):
    _, cooling, reference = boundary_problem
    positive_potential = Halo.PowerLaw(0, 100, 100, R_phi0_kpc=.001)
    assert positive_potential.phi_kms2(.8) > 0.
    solution = solve(positive_potential, cooling, reference.temperature_K[-1])
    assert solution.mass_flow_rate_Msun_per_yr == pytest.approx(
        reference.mass_flow_rate_Msun_per_yr, rel=2e-5)
    assert solution.temperature_K[-1] == pytest.approx(reference.temperature_K[-1], rel=2e-5)


@pytest.mark.parametrize('changes', [
    {'R_sonic_kpc': 0.}, {'R_boundary_kpc': .5},
    {'R_min_kpc': 2.}, {'T_boundary_K': -1.},
    {'T_boundary_K': np.nan}, {'R_boundary_kpc': np.inf},
])
def test_invalid_boundary_inputs_are_rejected(changes):
    parameters = dict(R_sonic_kpc=1., R_boundary_kpc=10.,
                      T_boundary_K=1e6, R_max_kpc=20., R_min_kpc=.8)
    parameters.update(changes)
    with pytest.raises(ValueError):
        CF.shoot_from_R_sonic_outer_temperature(
            Halo.PowerLaw(0, 100, 100), Cool.Constant_Cooling(1e-22), **parameters)


def test_impossible_target_is_not_reported_as_a_solution(boundary_problem):
    potential, cooling, _ = boundary_problem
    with pytest.raises((ValueError, RuntimeError)):
        solve(potential, cooling, 1e30)


def test_matching_temperature_does_not_guarantee_extension(boundary_problem):
    potential, cooling, _ = boundary_problem
    conditions = CF._calculate_sonic_point_conditions(.62, 1., potential, cooling, None)
    reference = CF._integrate_from_sonic_point(
        conditions, potential, cooling, 1., 10., .8, 1e-4,
        False, False, 2e4, .03, False, 'RK45', check_unbound=False)
    assert reference.stop_reason == CF.StopReason.MAX_RADIUS
    # This physical branch encounters a second sonic point near 15.6 kpc.
    with pytest.raises(RuntimeError, match='extension stops'):
        solve(potential, cooling, reference.temperature_K[-1], outer_radius=20.)


def test_exact_scan_node_is_one_solution(boundary_problem):
    potential, cooling, _ = boundary_problem
    x = float(np.linspace(1e-5, 1 - 1e-5, 33)[10])
    conditions = CF._calculate_sonic_point_conditions(x, 1., potential, cooling, None)
    reference = CF._integrate_from_sonic_point(
        conditions, potential, cooling, 1., 10., .8, 1e-4,
        False, False, 2e4, .03, False, 'RK45', check_unbound=False)
    solution = solve(potential, cooling, reference.temperature_K[-1])
    assert solution.mass_flow_rate_Msun_per_yr == pytest.approx(
        conditions['mass_flow_rate_Msun_per_yr'], rel=2e-5)
    assert solution.temperature_K[-1] == pytest.approx(reference.temperature_K[-1], rel=2e-5)
