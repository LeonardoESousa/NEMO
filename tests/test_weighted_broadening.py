import numpy as np
import pytest
from scipy.integrate import quad, trapezoid
from scipy.optimize import minimize_scalar
from scipy.stats import norm

from nemo import analysis


@pytest.mark.parametrize("solvent_sd", [0.0, 0.05, 0.2])
def test_weighted_bandwidth_matches_integrated_mise(solvent_sd):
    # Independent amplitudes: both weight groups sample the same reference density.
    e = np.repeat(norm.ppf((np.arange(30)+0.5)/30)*0.2, 2)
    w = np.tile([0.2, 1.8], 30)
    n = len(e)
    mean = np.average(e, weights=w)
    total_var = np.average((e-mean)**2, weights=w) + solvent_sd**2
    mean_w, second_w = np.mean(w), np.mean(w**2)

    def risk(h):
        individual = lambda x: norm.pdf(x, scale=np.sqrt(solvent_sd**2+h*h))
        expected = lambda x: mean_w*norm.pdf(x, scale=np.sqrt(total_var+h*h))
        target = lambda x: mean_w*norm.pdf(x, scale=np.sqrt(total_var))
        bound = 12*np.sqrt(total_var+h*h)
        self_norm = quad(lambda x: individual(x)**2, -bound, bound)[0]
        mean_norm = quad(lambda x: expected(x)**2, -bound, bound)[0]
        bias = quad(lambda x: (expected(x)-target(x))**2, -bound, bound)[0]
        return (second_w*self_norm-mean_norm)/n + bias

    oracle = minimize_scalar(risk, bounds=(1e-4, 1.0), method="bounded",
                             options={"xatol": 1e-10})
    widths = analysis.total_reorganization_energy(
        np.full(n, solvent_sd**2/(2*0.026)), 0.026, e, w
    )
    h = np.sqrt(max(0, widths[0]**2-solvent_sd**2))
    np.testing.assert_allclose(h, oracle.x, rtol=2e-5, atol=2e-7)


def test_dark_geometries_do_not_enter_bandwidth_moments():
    e = np.array([1., 1.2, 1.4, 5., 8.])
    w = np.array([1., 2., 1., 0., 0.])
    lam = np.array([0., 0.02, 0.1, 0., 0.1])
    before = analysis.total_reorganization_energy(lam, 0.026, e, w)
    e[3:] = [-1000., 2000.]
    lam[3:] = [10., 20.]
    after = analysis.total_reorganization_energy(lam, 0.026, e, w)
    np.testing.assert_array_equal(after[:3], before[:3])


def test_constant_weights_and_amplitude_scaling():
    e = np.linspace(1., 2., 20)
    lam = np.linspace(0., 0.2, len(e))
    uniform = analysis.total_reorganization_energy(lam, 0.026, e)
    np.testing.assert_array_equal(uniform, analysis.total_reorganization_energy(
        lam, 0.026, e, np.ones_like(e)))
    w = np.linspace(0.1, 1., len(e))
    a = analysis.total_reorganization_energy(lam, 0.026, e, w)
    b = analysis.total_reorganization_energy(lam, 0.026, e, 100*w)
    np.testing.assert_allclose(a, b, rtol=1e-6)


def test_multiple_bands_keep_independent_weights_and_common_added_variance():
    e = np.column_stack((np.linspace(1., 2., 20), np.linspace(2., 4., 20)))
    lam = np.column_stack((np.zeros(20), np.linspace(0., 0.2, 20)))
    w = np.column_stack((np.arange(20.), np.arange(20.)[::-1]))
    widths = analysis.total_reorganization_energy(lam, 0.026, e, w)
    assert np.isfinite(widths).all() and np.all(widths > 0)
    for j in range(2):
        np.testing.assert_allclose(widths[:, j], analysis.total_reorganization_energy(
            lam[:, j], 0.026, e[:, j], w[:, j]))
        np.testing.assert_allclose(np.ptp(widths[:, j]**2-2*lam[:, j]*0.026), 0, atol=1e-15)


@pytest.mark.parametrize("fallback", [0.28, 0.026])
def test_single_geometry_and_dark_channels_are_finite(fallback):
    lam = np.array([[0., 0.1]])
    expected = np.sqrt(2*lam*0.026+fallback**2)
    actual = analysis.total_reorganization_energy(
        lam, 0.026, np.array([[2., 3.]]), np.array([[0., 1.]]), fallback=fallback)
    np.testing.assert_allclose(actual, expected)
    lam = np.repeat(lam, 5, axis=0)
    actual = analysis.total_reorganization_energy(
        lam, 0.026, np.ones((5, 2))*2, np.zeros((5, 2)), fallback=fallback)
    np.testing.assert_allclose(actual, np.sqrt(2*lam*0.026+fallback**2))


@pytest.mark.parametrize("initial", ["S1", "T1"])
@pytest.mark.parametrize("dark", ["emission", "isc", "all"])
def test_dark_channels_keep_zero_rates_and_finite_reporting(initial, dark):
    from test_rates_ordering import crossed_states
    data = crossed_states()
    for col in data:
        if (dark in ("emission", "all") and col.startswith("osce_")) or (
                dark in ("isc", "all") and col.startswith("soc_")):
            data[col] = 0.0
    result, spectrum, details = analysis.rates(initial, (1., 1.), data=data, detailed=True)
    assert np.isfinite(result.select_dtypes(include=[np.number]).to_numpy()).all()
    assert np.isfinite(spectrum.to_numpy()).all()
    if dark in ("emission", "all"):
        assert spectrum.rate == 0
        assert result.iloc[0]["AvgConc(%)"] == 0
        np.testing.assert_array_equal(details[f"{initial}->S0"], 0)
    if dark in ("isc", "all"):
        np.testing.assert_array_equal(result.iloc[1:]["Rate(s^-1)"], 0)
        np.testing.assert_array_equal(result.iloc[1:]["AvgConc(%)"], 0)


@pytest.mark.parametrize("initial", ["S1", "T1"])
def test_process_amplitude_scaling_preserves_widths_and_ensemble_normalization(initial):
    from test_rates_ordering import crossed_states
    data = crossed_states()
    # Add dark geometries to exercise the original N in the final rate average.
    for col in data:
        if col.startswith(("osce_", "soc_")):
            data.loc[0, col] = 0.0
    before, emi_before, details_before = analysis.rates(initial, (3., 1.4), data=data.copy(), detailed=True)
    for col in data:
        if col.startswith("osce_"):
            data[col] *= 4
        elif col.startswith("soc_"):
            data[col] *= 2
    after, emi_after, details_after = analysis.rates(initial, (3., 1.4), data=data, detailed=True)
    np.testing.assert_allclose(after["Rate(s^-1)"], 4*before["Rate(s^-1)"], rtol=1e-8)
    np.testing.assert_allclose(after["AvgSigma(eV)"], before["AvgSigma(eV)"], rtol=1e-10)
    np.testing.assert_allclose(emi_after.rate, 4*emi_before.rate)
    np.testing.assert_allclose(details_after["sigma"], details_before["sigma"])
    np.testing.assert_allclose(before["Rate(s^-1)"], details_before.iloc[:, -len(before):].mean().to_numpy()
                               / np.r_[analysis.HBAR_EV, np.ones(len(before)-1)])


def test_absorption_uses_oscillator_weights_without_renormalizing_dark_geometries():
    from test_rates_ordering import crossed_states
    data = crossed_states()
    data["osc_s1"] = [0., 0.2]
    data["osc_s2"] = [0., 0.4]
    before, detail_before = analysis.absorption("S0", (3., 1.4), data=data.copy(), detailed=True)
    data[["osc_s1", "osc_s2"]] *= 4
    after, detail_after = analysis.absorption("S0", (3., 1.4), data=data.copy(), detailed=True)
    np.testing.assert_allclose(after.iloc[:, 1:], 4*before.iloc[:, 1:], rtol=1e-10)
    np.testing.assert_allclose(after.Energy, before.Energy)
    np.testing.assert_allclose(detail_after.filter(regex="^sigma_"), detail_before.filter(regex="^sigma_"))
    area = trapezoid(before.Total, before.Energy)
    expected = detail_before.filter(regex="^S0->").sum(axis=1).mean()
    np.testing.assert_allclose(area, expected, rtol=1e-3)
    data[["osc_s1", "osc_s2"]] = 0.
    dark = analysis.absorption("S0", (1., 1.), data=data)
    np.testing.assert_array_equal(dark.iloc[:, 1:], 0.)
