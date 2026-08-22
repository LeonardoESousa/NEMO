import numpy as np

from nemo import tools


def _two_fragment_molecule():
    geom = np.array([
        [0.0, 0.0, 0.0],
        [1.5, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [1.5, -1.0, 0.0],
    ])
    atoms = ["C", "C", "H", "H"]
    return geom, atoms, tools.adjacency(geom, atoms)


def test_rotate_fragment_preserves_internal_distances():
    geom, _, _ = _two_fragment_molecule()
    before = tools.distance_matrix(geom)

    tools._rotate_fragment(geom, 0, 1, np.array([1, 3]), np.pi / 2)
    after = tools.distance_matrix(geom)

    assert np.allclose(before[np.ix_([1, 3], [1, 3])],
                       after[np.ix_([1, 3], [1, 3])])
    assert np.allclose(geom[0], [0.0, 0.0, 0.0])
    assert np.allclose(geom[1], [1.5, 0.0, 0.0])


def test_ring_bond_is_not_a_rotor():
    adj = np.array([
        [0, 1, 1],
        [1, 0, 1],
        [1, 1, 0],
    ])

    assert tools._fragment_after_cut(adj, 0, 1) is None


def test_canonical_torsion_matches_finite_counterrotation():
    geom, atoms, adj = _two_fragment_molecule()
    bond = tools._find_rotatable_bonds(geom, atoms, adj)[0]
    field, rotor = tools._canonical_torsion(geom, *bond)
    atom_a, atom_b, fragment, other, weight_f, weight_o = rotor
    epsilon = 1e-6

    rotated = geom.copy()
    tools._rotate_fragment(
        rotated, atom_a, atom_b, fragment, weight_f * epsilon
    )
    tools._rotate_fragment(
        rotated, atom_a, atom_b, other, -weight_o * epsilon
    )

    assert np.allclose((rotated - geom) / epsilon, field, atol=1e-6)
    assert np.isclose(weight_f + weight_o, 1.0)


def test_subspace_combines_torsion_distributed_over_modes():
    geom, atoms, adj = _two_fragment_molecule()
    bond = tools._find_rotatable_bonds(geom, atoms, adj)[0]
    field, _ = tools._canonical_torsion(geom, *bond)
    modes = np.zeros((len(geom), 3, 3))
    modes[:, :, 0] = 0.6 * field
    modes[:, :, 1] = -0.8 * field
    modes[2, 0, 1] += 0.2
    modes[3, 1, 2] = 1.0
    frequencies = np.array([25.0, 80.0, 150.0])
    frequencies *= tools.LIGHT_SPEED * 100 * 2 * np.pi

    model = tools._build_torsional_subspace(
        geom, atoms, adj, frequencies, modes, cutoff_cm=100.0
    )

    assert np.array_equal(model["low_modes"], [0, 1])
    assert len(model["rotors"]) == 1
    assert model["angle_from_q"].shape == (1, 2)
    reconstructed = (
        model["residual_modes"].reshape(-1, 2)
        + field.reshape(-1, 1).dot(model["angle_from_q"])
    )
    assert np.allclose(reconstructed, modes[:, :, :2].reshape(-1, 2))


def test_sampled_subspace_keeps_q_and_connectivity():
    geom, atoms, adj = _two_fragment_molecule()
    frequencies = np.array([50.0])
    frequencies *= tools.LIGHT_SPEED * 100 * 2 * np.pi
    bond = tools._find_rotatable_bonds(geom, atoms, adj)[0]
    field, _ = tools._canonical_torsion(geom, *bond)
    modes = field[:, :, None]
    model = tools._build_torsional_subspace(
        geom, atoms, adj, frequencies, modes
    )
    scale = np.array([1.0])
    seed = 1234

    sampled, magnitudes, rejected = tools.sample_single_geometry(
        (geom, atoms, adj, scale, modes, model, seed, 10)
    )

    q = np.random.RandomState(seed).normal(scale=scale)[0]
    assert rejected == 0
    assert np.isclose(magnitudes[0, 0], q)
    assert np.array_equal(adj, tools.adjacency(sampled, atoms))
    assert np.isclose(np.linalg.norm(sampled[2] - sampled[0]), 1.0)


def test_empty_torsional_subspace_is_plain_linear_sampling():
    geom = np.array([[0.0, 0.0, 0.0]])
    atoms = ["He"]
    adj = np.zeros((1, 1))
    modes = np.zeros((1, 3, 1))
    modes[0, 0, 0] = 1.0
    frequencies = np.array([150.0])
    frequencies *= tools.LIGHT_SPEED * 100 * 2 * np.pi
    model = tools._build_torsional_subspace(
        geom, atoms, adj, frequencies, modes
    )
    seed = 1234

    sampled, magnitudes, rejected = tools.sample_single_geometry(
        (geom, atoms, adj, np.array([0.1]), modes, model, seed, 10)
    )

    q = np.random.RandomState(seed).normal(scale=[0.1])[0]
    assert rejected == 0
    assert np.allclose(sampled[0], [q, 0.0, 0.0])
    assert np.isclose(magnitudes[0, 0], q)
