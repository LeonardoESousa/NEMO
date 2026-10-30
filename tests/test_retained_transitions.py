from pathlib import Path

import numpy as np
import pytest

from nemo import analysis, eom, parser


def write_transitions(path, module, spin, count, missing=None):
    label = 'Singlet' if spin == 'S' else 'Triplet'
    # Interleaved combined IDs, as in TDDFT singlet/triplet calculations.
    ids = np.arange(1, 2 * count, 2)
    lines = []
    if module is parser:
        lines.append('TDDFT/TDA Excitation Energies')
        for root, state_id in enumerate(ids, 1):
            lines.extend([f'Excited state {state_id}: excitation energy (eV) = {root}',
                          f'Multiplicity: {label}', 'Strength : 0.1'])
        lines.extend(['---------------------------------------------------',
                      f'Transition Moments Between {label} Excited States',
                      'States X Y Z Strength(a.u.)', '-----------------------'])
    # Deliberately reverse row and pair order to verify target identities.
    for target in range(count, 1, -1):
        for source in range(target - 1, 0, -1):
            if (source, target) == missing:
                continue
            strength = source * 100 + target
            if module is parser:
                lines.append(f'{ids[target-1]} {ids[source-1]} 0 0 0 {strength}')
            else:
                kind = 'singlets' if spin == 'S' else 'triplets'
                lines.extend([f'State A: eomee_ccsd/rhfref/{kind}: {target}/A',
                              f'State B: eomee_ccsd/rhfref/{kind}: {source}/A',
                              f'Oscillator strength (a.u.): {strength}',
                              '-----------------------'])
    lines.extend(['-----------------------', 'Have a nice day'])
    if module is parser:
        # Later PCM tables use the same IDs but different strengths.
        lines.extend([f'Transition Moments Between {label} Excited States',
                      '1 3 0 0 0 99999', '-----------------------', 'Have a nice day'])
    path.write_text('\n'.join(lines))


@pytest.mark.parametrize('module', [parser, eom])
@pytest.mark.parametrize('spin', ['S', 'T'])
@pytest.mark.parametrize('root', [1, 2, 5])
@pytest.mark.parametrize('counts', [(5, 6), (6, 5)])
def test_retained_transition_values(tmp_path, monkeypatch, module, spin, root, counts):
    monkeypatch.chdir(tmp_path)
    Path('Geometries').mkdir()
    files = ['first.log', 'second.log']
    for file, count in zip(files, counts):
        write_transitions(Path('Geometries') / file, module, spin, count)
    actual = module.pega_oscs(files, f'{spin}{root}', nstates=5)
    expected = [[100 * root + t for t in range(root + 1, 6)]] * 2
    np.testing.assert_array_equal(actual, np.array(expected).reshape(2, 5-root))


@pytest.mark.parametrize('module', [parser, eom])
def test_missing_transition_is_not_zero_filled(tmp_path, monkeypatch, module):
    monkeypatch.chdir(tmp_path)
    Path('Geometries').mkdir()
    write_transitions(Path('Geometries/bad.log'), module, 'T', 6, missing=(1, 3))
    with pytest.raises(ValueError, match=r'bad.log: missing oscillator strengths from T1'):
        module.pega_oscs(['bad.log'], 'T1', nstates=5)


@pytest.mark.parametrize('module', [parser, eom])
@pytest.mark.parametrize('counts', [(5, 6), (6, 5)])
def test_soc_stacking_uses_retained_targets(monkeypatch, module, counts):
    values = {file: np.arange(1., count + 1)[None, :]
              for file, count in zip(['first.log', 'second.log'], counts)}
    monkeypatch.setattr(module, 'pega_soc_triplet', lambda file, root: values[file])
    np.testing.assert_array_equal(
        module.avg_socs(list(values), 'triplet', 0, nstates=5),
        [np.arange(1., 6)] * 2,
    )
    values['second.log'] = np.ones((1, 4))
    with pytest.raises(ValueError, match='second.log: missing SOC targets'):
        module.avg_socs(list(values), 'triplet', 0, nstates=5)


@pytest.mark.parametrize('module,calculation', [(parser, 'tddft'), (eom, 'eom-ccsd')])
@pytest.mark.parametrize('spin', ['S', 'T'])
def test_gather_columns_match_retained_energies(tmp_path, monkeypatch, module, calculation, spin):
    monkeypatch.chdir(tmp_path)
    Path('Geometries').mkdir()
    files = ['first.log', 'second.log']
    for file, count in zip(files, [6, 5]):
        write_transitions(Path('Geometries') / file, module, spin, count)
    energies = np.arange(1., 6)
    monkeypatch.setitem(analysis.get_energies, calculation,
                        lambda file: (energies, energies, energies, energies,
                                      energies, 0., energies, energies))
    monkeypatch.setitem(analysis.get_phosph_osc, calculation,
                        lambda *args: np.ones((1, 5)))
    for name in ['pega_soc_singlet', 'pega_soc_triplet']:
        monkeypatch.setattr(module, name,
                            lambda file, root: np.ones((1, 6 if file == files[0] else 5)))
    monkeypatch.setattr(module, 'pega_soc_ground', lambda file, root: np.ones((1, 1)))
    data, _ = analysis._build_gather_dataframe(
        spin + '1', files, 5, calculation, 3., 1.5)
    osc_columns = [c for c in data if c.startswith('osc_' + spin.lower())]
    assert osc_columns == [f'osc_{spin.lower()}{t}' for t in range(2, 6)]
    np.testing.assert_array_equal(data[osc_columns].to_numpy(), [[102, 103, 104, 105]] * 2)
