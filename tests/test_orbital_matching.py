from copy import deepcopy

import numpy as np
import pytest

from nemo import parser


def reference_data():
    return {
        'atoms': ['C', 'N', 'H', 'O'],
        'geometry': np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]),
        'tdm_singlets': np.array([[1., 0., 0.], [0., 2., 0.]]),
        'comp_singlets': [{(114, 1): 1.}, {(118, 2): 1.}],
        'comp_triplets': [{(114, 3): 0.9, (118, 1): 0.2}],
    }


@pytest.mark.parametrize('rotate', [False, True])
def test_singlet_inference_transfers_to_triplets_and_is_phase_and_frame_invariant(rotate):
    vacuum = reference_data()
    pcm = deepcopy(vacuum)
    pcm['comp_singlets'][0] = {(113, 1): 1.}
    pcm['comp_triplets'][0] = {(113, 3): 0.9, (118, 1): 0.2}
    if rotate:
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        pcm['geometry'] = pcm['geometry'] @ rotation + [2., -1., 4.]
        pcm['tdm_singlets'] = pcm['tdm_singlets'] @ rotation
    pcm['tdm_singlets'][0] *= -1  # The overall excited-state sign is arbitrary.
    occupied, virtual = parser._orbital_relabeling(vacuum, pcm)
    assert occupied == {113: 114, 114: 113}
    assert virtual == {}
    mapped = parser._relabel_compositions(pcm['comp_triplets'], occupied, virtual)
    assert parser.compute_similarity(vacuum['comp_triplets'][0], pcm['comp_triplets'][0]) < 0.01
    assert parser.compute_similarity(vacuum['comp_triplets'][0], mapped[0]) == pytest.approx(1.)


def test_virtual_permutation_uses_singlets_only():
    vacuum = reference_data()
    pcm = deepcopy(vacuum)
    vacuum['comp_singlets'][0] = {(118, 7): 1.}
    pcm['comp_singlets'][0] = {(118, 6): 1.}
    # Altering triplets must have no influence on the inferred orbital map.
    pcm['comp_triplets'] = [{(1, 99): 1.}]
    occupied, virtual = parser._orbital_relabeling(vacuum, pcm)
    assert occupied == {}
    assert virtual == {6: 7, 7: 6}


def test_occupied_and_virtual_labels_can_change_together():
    vacuum = reference_data()
    pcm = deepcopy(vacuum)
    pcm['comp_singlets'][0] = {(113, 3): 1.}
    occupied, virtual = parser._orbital_relabeling(vacuum, pcm)
    assert occupied == {113: 114, 114: 113}
    assert virtual == {1: 3, 3: 1}


@pytest.mark.parametrize('missing', ['dark', 'weak', 'moments', 'geometry', 'linear', 'changed_geometry', 'composition'])
def test_unavailable_or_unreliable_anchors_leave_original_labels(missing):
    vacuum = reference_data()
    pcm = deepcopy(vacuum)
    pcm['comp_singlets'][0] = {(113, 1): 1.}
    if missing == 'dark':
        pcm['tdm_singlets'][:] = 0.
    elif missing == 'weak':
        pcm['tdm_singlets'] *= 20.
    elif missing == 'moments':
        pcm['tdm_singlets'][:] = np.nan
    elif missing == 'geometry':
        pcm['geometry'] = np.empty((0, 3))
    elif missing == 'linear':
        vacuum['geometry'] = pcm['geometry'] = np.array([[0., 0., 0.], [1., 0., 0.],
                                                       [2., 0., 0.], [3., 0., 0.]])
    elif missing == 'changed_geometry':
        pcm['geometry'][1] += [0.5, 0., 0.]
    else:
        vacuum['comp_singlets'] = [{(114, 1): 0.}, {}]
    assert parser._orbital_relabeling(vacuum, pcm) == ({}, {})


def test_small_amplitude_changes_do_not_trigger_extra_orbital_swaps():
    vacuum = reference_data()
    pcm = deepcopy(vacuum)
    vacuum['comp_singlets'] = [{(112, 1): 0.70, (113, 1): 0.71}, {(118, 2): 1.}]
    pcm['comp_singlets'] = [{(112, 1): 0.71, (113, 1): 0.70}, {(118, 2): 1.}]
    assert parser._orbital_relabeling(vacuum, pcm) == ({}, {})


def test_one_consistent_frontier_cycle_and_virtual_swap_recover_triplet_order():
    vacuum = reference_data()
    pcm = deepcopy(vacuum)
    vacuum['tdm_singlets'] = pcm['tdm_singlets'] = np.diag([1., 2., 3.])
    vacuum['comp_singlets'] = [{(114, 1): 1.}, {(113, 1): 0.8, (118, 7): 0.6}, {(112, 1): 1.}]
    pcm['comp_singlets'] = [{(112, 1): 1.}, {(114, 1): 0.8, (118, 6): 0.6}, {(113, 1): 1.}]
    vacuum['comp_triplets'] = [{(114, 2): 0.9, (118, 7): 0.2}, {(113, 1): 1.}, {(112, 1): 1.}]
    pcm['comp_triplets'] = [{(114, 1): 1.}, {(113, 1): 1.}, {(112, 2): 0.9, (118, 6): 0.2}]
    occupied, virtual = parser._orbital_relabeling(vacuum, pcm)
    assert occupied == {112: 114, 113: 112, 114: 113}
    assert virtual == {6: 7, 7: 6}
    mapped = parser._relabel_compositions(pcm['comp_triplets'], occupied, virtual)
    roots, scores = parser.match_compositions(vacuum['comp_triplets'], mapped, return_scores=True)
    np.testing.assert_array_equal(roots, [2, 0, 1])
    np.testing.assert_allclose(scores[np.arange(3), roots], 1.)


@pytest.mark.parametrize('prefix, header', [
    ('', 'TDDFT/TDA Excitation Energies'),
    ('X: ', 'TDDFT Excitation Energies'),
])
def test_log_parsing_and_matching_apply_the_same_map_to_both_spins(tmp_path, prefix, header):
    def job(occupied, pcm=False):
        text = '''Total energy = -10.0
Standard Nuclear Orientation (Angstroms)
  1 C 0.0 0.0 0.0
  2 N 1.0 0.0 0.0
  3 H 0.0 1.0 0.0
  4 O 0.0 0.0 1.0
Nuclear Repulsion Energy = 1.0
TDDFT/TDA Excitation Energies
'''
        text = text.replace('TDDFT/TDA Excitation Energies', header)
        for root, (spin, virtual) in enumerate([('Triplet', 2), ('Singlet', 1)], 1):
            text += (f'Excited state {root}: excitation energy (eV) = {root}\n'
                     f'Multiplicity: {spin}\nTrans. Mom.: {int(spin == "Singlet")} X 0.0 Y 0.0 Z\n'
                     f'Strength: 0.1\n{prefix}D( {occupied}) --> V( {virtual}) amplitude = 1.0\n')
        text += '---------------------------------------------------\n'
        if pcm:
            text += 'Excited-state properties with   relaxed density\n'
            for _ in range(2):
                text += 'SS-PCM correction = -0.1 eV\nLR-PCM correction = -0.01 eV\n'
            text += '------------------------ END OF SUMMARY -----------------------\n'
        return text + 'Have a nice day\n'
    path = tmp_path / 'orbitals.log'
    path.write_text(job(114) + job(113, pcm=True))
    parsed = parser.pega_energias(path, return_matching=True)
    assert len(parser.pega_energias(path)) == 8
    for spin in ('s', 't'):
        np.testing.assert_array_equal(parsed[-1][spin]['pcm_roots'], [1])
        np.testing.assert_allclose(parsed[-1][spin]['scores'], [1.])


@pytest.mark.parametrize('spin', ['Singlet', 'Triplet'])
def test_tddft_y_amplitudes_do_not_overwrite_x_composition(spin):
    block = f'''TDDFT Excitation Energies
Excited state 1: excitation energy (eV) = 2.0
Multiplicity: {spin}
Trans. Mom.: 1.0 X 0.0 Y 0.0 Z
Strength: 0.1
    X: D( 89) --> V( 1) amplitude = -9.0E-1
    X: D( 88) --> V( 2) amplitude = 0.3
    Y: D( 89) --> V( 1) amplitude = 0.1
    Y: D( 87) --> V( 3) amplitude = 0.2
---------------------------------------------------
'''
    data = parser.parse_block(block)
    key = 'comp_singlets' if spin == 'Singlet' else 'comp_triplets'
    assert data[key][0] == {(89, 1): -0.9, (88, 2): 0.3}
