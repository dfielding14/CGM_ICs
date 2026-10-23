import json

import numpy as np
import pytest

import generate_cgm_library as library
import generate_gotham_ics as IC


@pytest.fixture(scope="module")
def fiducial(tmp_path_factory):
    source = tmp_path_factory.mktemp("fiducial") / "athinput"
    source.write_text("""<units>
length_cgs = 3.0856775809623245e21
mass_cgs = 3.036951775493658e39
time_cgs = 3.15576e15
mu = 0.62
<potential>
r_scale = 11.605294024177795
rho_scale = 2.1743800601559324
mass_gal = 6547.387043624801
scale_gal = 1.4691189466064514
z_gal = 0.15
r_200 = 144.76736317842207
rho_mean = 2.5759499304527406e-05
<mesh>
nx1 = 384
nx2 = 384
nx3 = 384
x1min = -200
x1max = 200
x2min = -200
x2max = 200
x3min = -200
x3max = 200
nghost = 4
<hydro>
gamma = 1.6666666666666667
<problem>
metallicity = 0.3333333333333333
<dust>
D_floor = 1e-7
""")
    config = IC.read_input(source)
    units = IC.CodeUnits(config)
    return source, IC.GalaxyPotential(config, units), IC.cooling_from_input(config)


def test_potential_table_force_and_slope(fiducial):
    direct = fiducial[1]
    table = IC.TabulatedPotential(direct, 0.001, 400)
    radius = np.geomspace(0.001, 400, 10003)
    np.testing.assert_allclose(table.phi_kms2(radius), direct.phi_kms2(radius), rtol=0, atol=1e-6)
    np.testing.assert_allclose(table.vc_kms(radius), direct.vc_kms(radius), rtol=1e-7)
    np.testing.assert_allclose(table.dln_vc_dln_r(radius), direct.dln_vc_dln_r(radius), rtol=0, atol=2e-5)
    step = 1e-4
    middle = radius[10:-10]
    force = (table.phi_kms2(middle * np.exp(step))
             - table.phi_kms2(middle * np.exp(-step))) / (2 * step)
    np.testing.assert_allclose(force, table.vc_kms(middle)**2, rtol=3e-7)
    assert np.isnan(table.phi_kms2(500))


@pytest.mark.parametrize("radius", [0.7, 35.0])
def test_adaptive_profiles_match_direct_small_steps(fiducial, radius):
    _, direct, cooling = fiducial
    table = IC.TabulatedPotential(direct, 0.1, 244.8)
    reference = IC.CF.shoot_from_R_sonic(direct, cooling, radius, 244.8, 0.1,
                                         max_step=0.002, epsilon=1e-4)
    fast = IC.CF.shoot_from_R_sonic(table, cooling, radius, 244.8, 0.1,
                                    max_step=np.inf, epsilon=1e-4, atol=1e-9, rtol=1e-9)
    assert fast.shooting_x == reference.shooting_x
    assert fast.mass_flow_rate_Msun_per_yr == pytest.approx(reference.mass_flow_rate_Msun_per_yr, rel=1e-7)
    assert fast.result.nfev < reference.result.nfev / 4
    fast.sample_profile()
    assert max(np.diff(np.log(fast.radius_kpc))) <= 0.002 * (1 + 1e-10)
    for name in ("temperature_K", "density_cgs", "velocity_kms"):
        expected = np.interp(np.log(fast.radius_kpc), np.log(reference.radius_kpc), np.log(getattr(reference, name)))
        assert max(abs(np.expm1(np.log(getattr(fast, name)) - expected))) < 2e-4
    flux = IC.us.mass_flow_rate_from_density_velocity_radius(fast.density_cgs, fast.velocity_kms, fast.radius_kpc)
    np.testing.assert_allclose(flux, fast.mass_flow_rate_Msun_per_yr, rtol=1e-12)
    warm = IC.CF.shoot_from_R_sonic(table, cooling, radius, 244.8, 0.1,
        max_step=np.inf, epsilon=1e-4, atol=1e-9, rtol=1e-9, x_guess=fast.shooting_x)
    assert warm.shooting_trials == 1
    assert warm.mass_flow_rate_Msun_per_yr == fast.mass_flow_rate_Msun_per_yr


def test_serial_parallel_and_exact_cache_agree(fiducial, tmp_path, monkeypatch):
    source = fiducial[0]
    serial = library.run_survey(source, tmp_path / "serial", r_max_kpc=244.8)
    parallel = library.run_survey(source, tmp_path / "parallel", r_max_kpc=244.8, workers=2)
    for a, b in zip(serial["results"], parallel["results"]):
        assert a["status"] == b["status"] == "accepted"
        assert not a["covers_mesh"]
        assert a["x_guess"] is b["x_guess"] is None
        with np.load(a["profile"]) as left, np.load(b["profile"]) as right:
            for field in ("radius_kpc", "temperature_K", "density_cgs", "velocity_kms"):
                np.testing.assert_array_equal(left[field], right[field])

    def forbidden_solve(*args, **kwargs):
        raise AssertionError("a cached profile must not be solved again")
    monkeypatch.setattr(IC.CF, "shoot_from_R_sonic", forbidden_solve)
    reused = library.run_survey(source, tmp_path / "serial", r_max_kpc=244.8, workers=2)
    assert all(entry["cache_hit"] for entry in reused["results"])
    with pytest.raises(ValueError, match="settings changed"):
        library.run_survey(source, tmp_path / "serial", r_max_kpc=245)
    assert json.loads((tmp_path / "serial/parameters.json").read_text())["gas_disk"] is False


def test_incomplete_high_fiducial_is_not_exported(fiducial, tmp_path):
    result = library.run_survey(fiducial[0], tmp_path, sonic_radii=[35])
    assert result["results"][0]["status"] == "no_accepted_shoot"
    assert not list(tmp_path.glob("rsonic_*.npz"))
    assert not list(tmp_path.glob("rsonic_*.txt"))


def test_warm_guess_retains_full_search_and_is_explicit(fiducial, tmp_path):
    _, direct, cooling = fiducial
    table = IC.TabulatedPotential(direct, 0.1, 244.8)
    solution = IC.CF.shoot_from_R_sonic(table, cooling, 0.7, 244.8, 0.1,
        max_step=np.inf, epsilon=1e-4, x_guess=0.1)
    assert solution is not None and solution.shooting_trials > 1
    result = library.run_survey(fiducial[0], tmp_path, sonic_radii=[0.7, 0.70007],
                                r_max_kpc=244.8, warm_start=True)
    first, second = result["results"]
    assert first["x_guess"] is None
    assert second["x_guess"] == first["shooting_x"]
    assert second["status"] == "accepted"


def test_numerically_failed_warm_guess_retries_original_search(fiducial, monkeypatch):
    _, direct, cooling = fiducial
    table = IC.TabulatedPotential(direct, 0.1, 244.8)
    original = IC.CF._integrate_from_sonic_point
    calls = 0

    def fail_once(*args, **kwargs):
        nonlocal calls
        calls += 1
        if calls == 1:
            raise RuntimeError("failed warm trial")
        return original(*args, **kwargs)

    monkeypatch.setattr(IC.CF, "_integrate_from_sonic_point", fail_once)
    result = IC.CF.shoot_from_R_sonic(table, cooling, .7, 244.8, .1,
        max_step=np.inf, epsilon=1e-4, atol=1e-9, rtol=1e-9, x_guess=.438)
    assert result.shooting_x == pytest.approx(0.4379939013671875, abs=1e-15)
    assert result.shooting_trials == 12
