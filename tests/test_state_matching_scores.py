from pathlib import Path

import numpy as np
import pytest

from nemo import parser, tools


@pytest.mark.parametrize('vacuum,pcm,expected', [
    ([{(1, 1): 1.}, {(2, 1): 1.}],
     [{(2, 1): -1.}, {(1, 1): 1.}], [1, 0]),
    ([{(1, 1): 1.}, {(2, 1): 1.}],
     [{(1, 1): 1.}, {(3, 1): 1.}], [0, 1]),
    ([{(1, 1): 1.}], [{(2, 1): 1.}, {(1, 1): 1.}], [1]),
    ([], [{(1, 1): 1.}], []),
    ([{(1, 1): 1.}], [], []),
])
def test_scores_describe_actual_assignment_without_changing_default(vacuum, pcm, expected):
    old_api = parser.match_compositions(vacuum, pcm)
    indices, scores = parser.match_compositions(vacuum, pcm, return_scores=True)
    np.testing.assert_array_equal(indices, expected)
    np.testing.assert_array_equal(indices, old_api)
    assert scores.shape == (len(vacuum), len(pcm))
    for i in range(len(vacuum)):
        for j in range(len(pcm)):
            assert scores[i, j] == pytest.approx(parser.compute_similarity(vacuum[i], pcm[j]))


def write_jobs(tmp_path, vacuum, pcm, pcm_energies=None):
    def block(compositions, energies, corrections=False):
        lines = ['Total energy = -10.0', 'Dielectric 3.0',
                 'OpticalDielectric 2.0', 'TDDFT/TDA Excitation Energies']
        for root, (comp, energy) in enumerate(zip(compositions, energies), 1):
            lines.extend([f'Excited state {root}: excitation energy (eV) = {energy}',
                          'Multiplicity: Singlet', 'Strength: 0.1'])
            for (occupied, virtual), amplitude in comp.items():
                lines.append(f'D( {occupied}) --> V( {virtual}) amplitude = {amplitude}')
        lines.append('---------------------------------------------------')
        if corrections:
            lines.append('Excited-state properties with   relaxed density')
            for _ in compositions:
                lines.extend(['SS-PCM correction = -0.1 eV',
                              'LR-PCM correction = -0.02 eV'])
            lines.append('------------------------ END OF SUMMARY -----------------------')
        return '\n'.join(lines) + '\nHave a nice day\n'
    vac_energies = np.arange(len(vacuum)) + 2.
    if pcm_energies is None:
        pcm_energies = np.arange(len(pcm)) + 1.9
    content = block(vacuum, vac_energies)
    content += block(pcm, pcm_energies, corrections=True)
    path = tmp_path / 'matching.log'
    path.write_text(content)
    return path


def test_zero_score_for_missing_counterpart_is_reported_without_extra_output(tmp_path, capsys):
    log = write_jobs(tmp_path, [{(1, 1): 1.}, {(2, 1): 1.}],
                     [{(1, 1): 1.}, {(3, 1): 1.}])
    parsed = parser.pega_energias(log, return_matching=True)
    np.testing.assert_array_equal(parsed[-1]['s']['scores'], [1., 0.])
    tools.susceptibility_check(log)
    output = capsys.readouterr().out
    assert 'E_vac(eV)' in output and 'MatchScore' in output
    rows = [line.split() for line in output.splitlines()]
    assert len(rows) == 4  # Header, S0, S1 and the poorly paired S2.
    assert rows[-1] == ['S2', '3.000', '0.360', '0.200', '0.5795']
    assert 'PCM' not in output and 'WARNING' not in output and 'NOTE' not in output


def test_matching_uses_all_calculated_pcm_candidates(tmp_path, capsys):
    log = write_jobs(tmp_path, [{(1, 1): 1.}, {(2, 1): 1.}],
                     [{(1, 1): 1.}, {(3, 1): 1.}, {(2, 1): 1.}])
    default = parser.pega_energias(log)
    diagnostic = parser.pega_energias(log, return_matching=True)
    assert len(default) == 8 and len(diagnostic) == 9
    for old, new in zip(default, diagnostic):
        np.testing.assert_array_equal(old, new)
    info = diagnostic[-1]['s']
    assert info['matrix'].shape == (2, 3)
    assert info['matrix'][1, 2] == 1.
    np.testing.assert_array_equal(info['pcm_roots'], [1, 3])
    assert info['scores'][1] == 1.
    tools.susceptibility_check(log)
    assert 'Composition similarities' not in capsys.readouterr().out


def test_raw_negative_shift_is_available_without_extra_console_output(tmp_path, capsys):
    log = write_jobs(tmp_path, [{(1, 1): 1.}], [{(1, 1): 1.}], [2.5])
    diagnostic = parser.pega_energias(log, return_matching=True)
    np.testing.assert_allclose(diagnostic[6], [0.])
    np.testing.assert_allclose(diagnostic[-1]['s']['raw_shift'], [-0.5])
    tools.susceptibility_check(log)
    assert len(capsys.readouterr().out.splitlines()) == 3


def test_single_job_does_not_claim_independent_match(tmp_path, capsys):
    log = write_jobs(tmp_path, [{(1, 1): 1.}], [{(1, 1): 1.}])
    # The second job supplies the PCM correction format also needed for a single job.
    log.write_text(log.read_text().split('Have a nice day')[1] + 'Have a nice day\n')
    tools.susceptibility_check(log)
    assert len(capsys.readouterr().out.splitlines()) == 3


def test_real_log_report_preserves_parser_outputs_and_selected_pcm_roots(capsys):
    log = Path(__file__).parent / 'tddft/Geometries/Geometry-2-.log'
    default = parser.pega_energias(log)
    diagnostic = parser.pega_energias(log, return_matching=True)
    for old, new in zip(default, diagnostic):
        np.testing.assert_array_equal(old, new)
    info = diagnostic[-1]['s']
    np.testing.assert_array_equal(info['pcm_roots'], [1, 2, 3, 5, 4])
    np.testing.assert_allclose(info['scores'],
                               info['matrix'][np.arange(5), info['pcm_roots']-1])
    tools.susceptibility_check(log)
    table = capsys.readouterr().out
    rows = {line.split()[0]: line.split() for line in table.splitlines()
            if len(line.split()) == 5 and line.split()[0] in ('S4', 'S5')}
    assert float(rows['S4'][4]) == pytest.approx(info['assignment_scores'][3], abs=5e-5)
    assert float(rows['S5'][4]) == pytest.approx(info['assignment_scores'][4], abs=5e-5)
