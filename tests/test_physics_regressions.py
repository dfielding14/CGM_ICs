import h5py
import numpy as np
import pytest
from scipy.integrate import cumulative_trapezoid
from types import SimpleNamespace

import HaloPotential_new as Halo
import WiersmaCooling as Cool
import solve_ode as CF
import unit_system as us


def test_wiersma_keeps_signed_net_rates_and_matches_gotham_table():
    cooling = Cool.Wiersma_Cooling(1.0 / 3.0, 0.0)
    with h5py.File(Cool.DATA_DIR / 'z_0.000.hdf5') as table:
        temperatures = table['Metal_free/Temperature_bins'][...]
        densities = table['Metal_free/Hydrogen_density_bins'][...]
        expected = (table['Metal_free/Net_Cooling'][1, ...].astype(float)
                    + table['Total_Metals/Net_cooling'][...].astype(float) / 3.0)
    actual = cooling.lambda_cgs(temperatures[:, None], densities[None, :])
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=1e-35)
    assert cooling.lambda_cgs(1e4, 1e-4) < 0.0


def test_cooling_slopes_differentiate_the_evaluated_rate():
    cooling = Cool.Wiersma_Cooling(1.0 / 3.0, 0.0)
    temperature, density = 234567.0, 0.0007654
    step = 1e-6
    rate = cooling.lambda_cgs(temperature, density)
    for axis, slope in enumerate((cooling.dln_lambda_dln_T,
                                  cooling.dln_lambda_dln_rho)):
        lo, hi = np.array([temperature, density]), np.array([temperature, density])
        lo[axis] *= np.exp(-step)
        hi[axis] *= np.exp(step)
        expected = (cooling.lambda_cgs(*hi) - cooling.lambda_cgs(*lo)) / (2 * step * rate)
        assert slope(temperature, density) == pytest.approx(expected, rel=1e-6)


def test_cooling_clamps_high_density_and_rejects_invalid_inputs():
    cooling = Cool.Wiersma_Cooling(1.0 / 3.0, 0.0)
    assert cooling.lambda_cgs(1e6, 1e5) == cooling.lambda_cgs(1e6, 1.0)
    assert cooling.dln_lambda_dln_rho(1e6, 1e5) == 0.0
    with pytest.raises(ValueError):
        cooling.lambda_cgs(0.0, 1e-4)
    with pytest.raises(ValueError):
        cooling.dln_lambda_dln_T(1e4, 1e-4)


def test_failed_integrator_is_not_a_success(monkeypatch):
    original = CF.scipy.integrate.solve_ivp

    def failed(*args, **kwargs):
        result = original(*args, **kwargs)
        result.success, result.status = False, -1
        result.message = 'Required step size is less than spacing between numbers.'
        return result

    monkeypatch.setattr(CF.scipy.integrate, 'solve_ivp', failed)
    with pytest.raises(RuntimeError, match='Required step size'):
        CF.IntegrateFlowEquations(.01, 1e6, 1e-27,
                                 Halo.PowerLaw(0, 100, 100),
                                 Cool.Constant_Cooling(1e-22), 1, 10, 20,
                                 check_unbound=False)


def test_nonterminal_unbound_flow_can_reach_outer_boundary():
    solution = CF.IntegrateFlowEquations(.01, 1e6, 1e-27,
                                        Halo.PowerLaw(0, 10, 100),
                                        Cool.Constant_Cooling(1e-22), 1, 10, 20,
                                        terminal_unbound=False)
    assert solution.stop_reason == CF.StopReason.MAX_RADIUS
    assert solution.radius_kpc[-1] == pytest.approx(20)


def test_gas_mass_is_integral_over_the_resolved_radial_interval():
    solution = CF.IntegrateFlowEquations(.01, 1e6, 1e-27,
                                        Halo.PowerLaw(0, 100, 100),
                                        Cool.Constant_Cooling(1e-22), 1, 10, 20)
    r = solution.radius_kpc * us.KPC_TO_CM
    expected = cumulative_trapezoid(4 * np.pi * r**2 * solution.density_cgs,
                                    r, initial=0) / us.MSUN_TO_G
    np.testing.assert_allclose(solution.gas_mass_Msun, expected, rtol=1e-12)


def test_rotating_flow_conserves_specific_angular_momentum():
    potential = Halo.PowerLaw(-0.2, 100, 100)
    radius, rcirc, temperature, density = 20.0, 5.0, 1e6, 1e-27
    y = np.log([temperature, density])
    cooling = Cool.Constant_Cooling(1e-22)
    static = CF._create_ode_system(.01, potential, cooling, 0, 1)(np.log(radius), y)
    rotating = CF._create_ode_system(.01, potential, cooling, rcirc, 1)(np.log(radius), y)
    cs2 = us.sound_speed_squared_kms2(temperature, CF.GC.GAMMA, CF.GC.MU)
    speed = us.mass_flow_to_velocity_kms(.01, radius, density)
    j = potential.vc_kms(rcirc) * rcirc
    expected = (j / radius)**2 / cs2 / (1 - speed**2 / cs2)
    assert rotating[1] - static[1] == pytest.approx(expected, rel=1e-12)


def test_sonic_conditions_satisfy_regularized_momentum_equation():
    potential = Halo.PowerLaw(-0.1, 100, 100)
    for x in (0.2, 0.5):
        conditions = CF._calculate_sonic_point_conditions(x, 10, potential,
                                                         Cool.Constant_Cooling(1e-22), -1)
        a, q = conditions['dlnT_dlnR'], conditions['dlnrho_dlnR']
        t = conditions['tflow_to_tcool']
        # Differentiate the momentum equation at M=1 and use continuity + entropy.
        lhs = (q + 2) * (-4 - 2*q - a)
        rhs = (t / CF.GC.GAMMA) * (3 + 2*q - a) + 2*x*(-0.2 - a)
        assert lhs == pytest.approx(rhs, abs=1e-12)


def test_sonic_density_search_resolves_dilute_flows():
    cooling = Cool.Constant_Cooling(1e-22)
    nH = 1e-10
    density = nH * us.M_P_G / CF.GC.X
    temperature, radius, velocity = 1e6, 10.0, 100.0
    cs2 = us.sound_speed_squared_kms2(temperature, CF.GC.GAMMA, CF.GC.MU)
    ratio = us.flow_time_Gyr(radius, velocity) / us.cooling_time_Gyr(
        density, cs2, nH, 1e-22, CF.GC.GAMMA, CF.GC.GAMMA_M1)
    actual = CF._calc_rho_from_tflow_to_tcool(velocity, ratio, temperature, radius, cooling)
    assert actual == pytest.approx(density, rel=1e-10, abs=0)


def test_cooling_radius_ignores_heating_and_does_not_extrapolate():
    solution = SimpleNamespace(radius_kpc=np.array([1., 2., 4., 8., 16.]),
                               cooling_time_Gyr=np.array([-1., .1, 10., .1, 10.]))
    assert CF.CGMSolution.cooling_radius_kpc(solution, 1.) == pytest.approx(np.sqrt(8 * 16))
    assert np.isnan(CF.CGMSolution.cooling_radius_kpc(solution, 100.))
