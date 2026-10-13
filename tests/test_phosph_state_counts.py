from pathlib import Path

import numpy as np
import pytest

from nemo import parser


@pytest.mark.parametrize("n_singlets,n_triplets", [(5, 5), (4, 5), (5, 4)])
def test_moment_uses_retained_states(monkeypatch, n_singlets, n_triplets):
    # The vacuum job prints six states of each spin; fewer survive gathering.
    ground_soc = (np.arange(1, 7) + 0.5j)[:, None] * 1e-4
    excited_soc = (np.arange(2, 8) - 0.3j)[:, None] * 1e-4
    monkeypatch.setattr(parser, "soc_s0", lambda *args: ground_soc)
    monkeypatch.setattr(parser, "soc_t1", lambda *args: excited_soc)
    singlets = np.arange(4., 10.)[:n_singlets]
    triplets = np.arange(1., 7.)[:n_triplets] * 0.4
    ground_dipoles = np.arange(1., 22.).reshape(7, 3) / 10
    triplet_dipoles = np.arange(1., 19.).reshape(6, 3) / 7
    expected = np.zeros(3, dtype=complex)
    for t in range(n_triplets):
        expected += ground_soc[t, 0] / (-triplets[t]) * triplet_dipoles[t]
    expected += ground_soc[0, 0] / triplets[0] * ground_dipoles[0]
    for s in range(n_singlets):
        expected += (excited_soc[s, 0].conjugate()
                     / (triplets[0] - singlets[s]) * ground_dipoles[s + 1])
    # Identical mocked SOCs for each of the three magnetic sublevels.
    expected = 3 * np.sum(np.abs(expected)**2) * (8.4783533e-30)**2
    actual = parser.moment("test.log", singlets, triplets, 0,
                           ground_dipoles, triplet_dipoles, 0)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)


def test_phosph_osc_ignores_later_jobs(tmp_path, monkeypatch):
    source = Path(__file__).parent / "tddft/Geometries/Geometry-1-.log"
    first_job = source.read_text().split("Have a nice day")[0] + "Have a nice day\n"
    geometries = tmp_path / "Geometries"
    geometries.mkdir()
    (geometries / "one.log").write_text(first_job)
    # A later job with another copy of the same tables must not add SOCs/dipoles.
    (geometries / "two.log").write_text(first_job + first_job)
    monkeypatch.chdir(tmp_path)
    data = parser.parse_block(first_job)
    singlets, triplets = data["singlets"], data["triplets"]
    expected = parser.phosph_osc("one.log", len(triplets), singlets, triplets)
    actual = parser.phosph_osc("two.log", len(triplets), singlets, triplets)
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0)
    # Exercise the real Q-Chem tables with both spin manifolds truncated.
    reduced = parser.phosph_osc("two.log", 4, singlets[:4], triplets[:4])
    assert reduced.shape == (1, 4)
    assert np.all(np.isfinite(reduced))
    assert np.all(reduced > 0)
