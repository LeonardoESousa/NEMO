from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import nemo.analysis as analysis


def crossed_states():
    data = pd.DataFrame({"kbT": [0.1, 0.1], "gamma_s0": [0.0, 0.0]})
    energies = {"s": [[2.0, 2.35], [2.0, 2.5]],
                "t": [[2.3, 2.6, 2.2], [2.3, 2.1, 2.4]]}
    chis = {"s": [[0.1, 1.0], [0.1, 1.0]],
            "t": [[0.1, 1.0, 0.4], [0.1, 0.5, 0.3]]}
    for spin in ("s", "t"):
        for i in range(len(energies[spin][0])):
            data[f"e_{spin}{i+1}"] = np.array(energies[spin])[:, i]
            data[f"chi_{spin}{i+1}"] = np.array(chis[spin])[:, i]
            data[f"gamma_{spin}{i+1}"] = 0.0
            data[f"osce_{spin}{i+1}"] = 0.01 * (i + 1)
    for i in range(2):
        for j in range(3):
            data[f"soc_s{i+1}_t{j+1}"] = 0.001 * (3 * i + j + 1)
            data[f"soc_t{j+1}_s{i+1}"] = 0.001 * (3 * i + j + 1)
    for i in range(3):
        data[f"soc_t{i+1}_s0"] = 0.002 * (i + 1)
    return data


def test_reorder_rectangular_axes_and_input_preservation():
    initial = np.array([[2., 1.], [1., 2.]])
    final = np.array([[3., 1., 2.], [2., 3., 1.]])
    socs = np.arange(12.).reshape(2, 6)
    original = socs.copy()
    actual = analysis.reorder(initial, final, initial, final, initial, final, socs)[-1]
    expected = np.array([[[4., 5., 3.], [1., 2., 0.]],
                         [[8., 6., 7.], [11., 9., 10.]]])
    np.testing.assert_array_equal(actual, expected)
    np.testing.assert_array_equal(socs, original)


@pytest.mark.parametrize("initial", ["S1", "S2", "T1", "T2"])
def test_rates_keep_one_state_identity_for_every_process(initial):
    data = crossed_states()
    spin = initial[0].lower()
    other = "t" if spin == "s" else "s"
    rank = int(initial[1:]) - 1
    eps, nr = 3., np.sqrt(2.)
    ast, aopt = 0.5, 1. / 3.
    arrays = {}
    orders = {}
    for label in (spin, other):
        energies = analysis.fetch(data, [f"^e_{label}"])
        chi = analysis.fetch(data, [f"^chi_{label}"])
        arrays[label] = energies, chi
        # The starting state is equilibrated; the endpoint uses final_state.
        alpha = ast if label == spin else aopt
        orders[label] = np.argsort(energies - chi * alpha, axis=1, kind="stable")
    assert np.any(orders[other] != np.argsort(
        arrays[other][0] - arrays[other][1]*ast, axis=1, kind="stable"
    ))
    if spin == "s":
        data = data.drop(columns=[c for c in data if c.startswith(("osce_t", "soc_t"))])
    else:
        data = data.drop(columns=[c for c in data if c.startswith(("osce_s", "soc_s"))])
    _, emission, details = analysis.rates(initial, (eps, nr), data=data, detailed=True)
    constant = nr**2 * analysis.E_CHARGE**2 / (
        2*np.pi*analysis.HBAR_EV*analysis.MASS_E*analysis.LIGHT_SPEED**3*analysis.EPSILON_0
    )
    if spin == "t": constant /= 3
    expected_radiative = []
    for row in range(len(data)):
        source = orders[spin][row, rank]
        ei, ci = arrays[spin][0][row, source], arrays[spin][1][row, source]
        photon = ei - ci*ast - ci*(ast-aopt)
        width = np.sqrt(2*ci*(ast-aopt)*data.kbT[row] + data.kbT[row]**2)
        np.testing.assert_allclose(details.eng[row], photon)
        np.testing.assert_allclose(details.sigma[row], width)
        expected_radiative.append(constant*photon**3/ei*data[f"osce_{spin}{source+1}"][row]/analysis.HBAR_EV)
        for target_rank, target in enumerate(orders[other][row]):
            ef, cf = arrays[other][0][row, target], arrays[other][1][row, target]
            gap = ef - cf*aopt - (ei-ci*ast)
            sigma = np.sqrt(2*cf*(ast-aopt)*data.kbT[row] + data.kbT[row]**2)
            h = data[f"soc_{spin}{source+1}_{other}{target+1}"][row]
            expected = 2*np.pi/analysis.HBAR_EV*h**2*np.exp(-gap**2/(2*sigma**2))/(np.sqrt(2*np.pi)*sigma)
            np.testing.assert_allclose(details[f"{initial}~>{other.upper()}{target_rank+1}"][row], expected, rtol=1e-12, atol=0)
        if spin == "t":
            h = data[f"soc_t{source+1}_s0"][row]
            expected = 2*np.pi/analysis.HBAR_EV*h**2*np.exp(-photon**2/(2*width**2))/(np.sqrt(2*np.pi)*width)
            np.testing.assert_allclose(details[f"{initial}~>S0"][row], expected, rtol=1e-12, atol=0)
    np.testing.assert_allclose(emission.rate, np.mean(expected_radiative))
    shuffled = data.sample(frac=1, axis=1, random_state=41)
    result, spectrum, breakdown = analysis.rates(initial, (eps, nr), data=shuffled, detailed=True)
    expected_result, expected_spectrum, expected_breakdown = analysis.rates(initial, (eps, nr), data=data, detailed=True)
    pd.testing.assert_frame_equal(result, expected_result)
    pd.testing.assert_frame_equal(spectrum, expected_spectrum)
    pd.testing.assert_frame_equal(breakdown, expected_breakdown)


def test_fetch_orders_multidigit_states_numerically():
    data = pd.DataFrame({f"soc_t{i}_s{j}": [100*i+j]
                         for i in (10, 2, 1) for j in (10, 2, 1)})
    np.testing.assert_array_equal(analysis.fetch(data, [r"^soc_t\d+_s[1-9]\d*$"]),
                                 [[101, 102, 110, 201, 202, 210, 1001, 1002, 1010]])


@pytest.mark.parametrize("initial,spin", [("S2", "s"), ("T2", "t")])
def test_gather_names_full_emission_arrays_from_root_one(initial, spin, monkeypatch):
    monkeypatch.chdir(Path(__file__).resolve().parent / "tddft")
    frame, _ = analysis._build_gather_dataframe(
        initial, ["Geometry-1-.log"], 5, "tddft", 2.38, 1.49
    )
    assert [c for c in frame if c.startswith("osce_")] == [f"osce_{spin}{i}" for i in range(1, 6)]
