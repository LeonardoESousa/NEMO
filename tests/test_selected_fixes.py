import numpy as np
import pytest

from nemo import analysis, eom
from test_rates_ordering import crossed_states


def write_log(tmp_path, monkeypatch, text):
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'Geometries').mkdir()
    (tmp_path / 'Geometries' / 'test.log').write_text(text)
    return 'test.log'


def test_ground_dipole_is_in_atomic_units(tmp_path, monkeypatch):
    file = write_log(tmp_path, monkeypatch,
                     'Dipole Moment (Debye)\nX 2.541746 Y -5.083492 Z 7.625238\n')
    np.testing.assert_allclose(eom.pega_dipole_ground(file), [[1., -2., 3.]])


@pytest.mark.parametrize('root', [1, 2, 11, 12])
def test_triplet_soc_selects_exact_root(tmp_path, monkeypatch, root):
    log = ''
    for printed_root in [11, 12, 1, 2]:
        log += (f'State A: eomee_ccsd/rhfref/singlets: 1/A\n'
                f'State B: eomee_ccsd/rhfref/triplets: {printed_root}/A\n'
                'Arithmetically averaged transition SO matrices\n'
                f'Hso(L+) = ({printed_root},0)\n'
                f'SOCC = {printed_root}\n'
                '--------------------------------------\n')
    file = write_log(tmp_path, monkeypatch, log)
    expected = [[root * 0.12398 / 1000]]
    np.testing.assert_allclose(eom.pega_soc_triplet(file, root - 1), expected)
    np.testing.assert_allclose(eom.soc_t1(file, '1', root - 1), expected)


@pytest.mark.parametrize('root', [1, 2, 11, 12])
def test_ground_soc_selects_exact_root(tmp_path, monkeypatch, root):
    log = ''
    for printed_root in [11, 12, 1, 2]:
        log += ('State A: ccsd: 0/A\n'
                f'State B: eomee_ccsd/rhfref/triplets: {printed_root}/A\n'
                'Arithmetically averaged transition SO matrices\n'
                f'SOCC = {printed_root}\n')
    file = write_log(tmp_path, monkeypatch, log)
    np.testing.assert_allclose(eom.pega_soc_ground(file, root - 1),
                               [[root * 0.12398 / 1000]])


@pytest.mark.parametrize('initial', ['S1', 'T1'])
def test_reported_emission_coupling(initial):
    data = crossed_states().iloc[[0, 0]].reset_index(drop=True)
    results, emission = analysis.rates(initial, (3., np.sqrt(2.)), data=data)
    expected = 1000 * np.sqrt(analysis.HBAR_EV * emission.rate / (2 * np.pi))
    np.testing.assert_allclose(results.iloc[0]['AvgCoupling(meV)'], expected)
