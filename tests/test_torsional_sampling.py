"""Physical and numerical regression checks; also runnable with unittest."""
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from joblib import Parallel

from nemo import tools


class TorsionalSamplingTests(unittest.TestCase):
    def setUp(self):
        self.geom = np.array([[0., 0., 0.], [1.5, 0., 0.],
                              [0., 1., 0.], [1.5, -1., 0.]])
        self.atoms = ["C", "C", "H", "H"]
        self.mass = np.array([12., 12., 1., 2.])
        self.adj = tools.adjacency(self.geom, self.atoms)
        self.freq = np.array([50.]) * tools.LIGHT_SPEED * 200 * np.pi
        bond = tools._find_rotatable_bonds(self.geom, self.atoms, self.adj)[0]
        field, self.rotor = tools._canonical_torsion(self.geom, *bond, self.mass)
        root = np.repeat(np.sqrt(self.mass), 3)
        rigid = tools._rigid_basis(self.geom, self.mass)
        weighted = root * field.ravel()
        self.mode = ((weighted - rigid @ (rigid.T @ weighted)) / root).reshape(4, 3, 1)

    def model(self, **kwargs):
        return tools._build_torsional_subspace(
            self.geom, self.atoms, self.adj, self.freq, self.mode,
            atomic_masses=self.mass, **kwargs)

    def test_mass_weighted_pure_torsion_has_no_linear_residual(self):
        model = self.model()
        np.testing.assert_allclose(model["angle_from_q"], [[1.]], atol=1e-12)
        np.testing.assert_allclose(model["residual_modes"], 0., atol=1e-12)
        self.assertAlmostEqual(self.rotor[-2], 1 / 3)
        self.assertAlmostEqual(self.rotor[-1], 2 / 3)

    def test_small_amplitude_limit_includes_residual_motion(self):
        self.mode[2, 1, 0] += 0.15
        model = self.model()
        for scale in (1e-4, 1e-5):
            sampled, q, rejected = tools.sample_single_geometry(
                (self.geom, self.atoms, self.adj, np.array([scale]),
                 self.mode, model, 1234, 10))
            self.assertEqual(rejected, 0)
            np.testing.assert_allclose((sampled - self.geom) / q[0, 0],
                                       self.mode[:, :, 0], atol=scale)

    def test_large_pure_torsion_preserves_bonds_and_sampled_amplitudes(self):
        model = self.model()
        i, j = np.where(np.triu(self.adj, 1))
        reference = np.linalg.norm(self.geom[i] - self.geom[j], axis=1)
        for seed in range(100):
            sampled, q, rejected = tools.sample_single_geometry(
                (self.geom, self.atoms, self.adj, np.array([2.]),
                 self.mode, model, seed, 10))
            self.assertEqual(rejected, 0)
            self.assertEqual(q[0, 0], np.random.RandomState(seed).normal(scale=2.))
            np.testing.assert_allclose(np.linalg.norm(sampled[i] - sampled[j], axis=1),
                                       reference, atol=1e-12)

    def test_above_cutoff_requires_character_and_angular_amplitude(self):
        self.freq *= 4
        self.assertEqual(len(self.model(scales=np.array([0.3]))["low_modes"]), 1)
        self.assertEqual(len(self.model(scales=np.array([0.01]))["low_modes"]), 0)
        self.assertEqual(len(self.model(scales=np.array([0.3]), min_angle=None)["low_modes"]), 0)
        self.mode[:] = 0
        self.mode[2, 1, 0] = 1.
        self.assertEqual(len(self.model(scales=np.array([0.3]))["low_modes"]), 0)

    def test_compressed_bond_and_nonfinite_geometry_are_rejected(self):
        bad = self.geom.copy()
        bad[2, 1] = 0.2
        self.assertTrue(np.array_equal(tools.adjacency(bad, self.atoms), self.adj))
        self.assertFalse(tools._sampling_geometry_valid(bad, self.geom, self.atoms, self.adj))
        bad[2, 1] = np.nan
        self.assertFalse(tools._sampling_geometry_valid(bad, self.geom, self.atoms, self.adj))

    def test_rigid_motion_is_not_mistaken_for_internal_rotation(self):
        rigid_mode = np.cross(np.array([0., 0., 1.]), self.geom)
        self.mode = rigid_mode[:, :, None]
        np.testing.assert_allclose(self.model()["angle_from_q"], 0., atol=1e-12)

    def test_invalid_inputs_fail_before_sampling(self):
        for kwargs in ({"temp": -1}, {"num_geoms": 0}, {"max_attempts": 0},
                       {"rotor_svd_cutoff": 0}, {"bond_tolerance": 1}):
            args = dict(freqlog="missing.log", num_geoms=1, temp=300)
            args.update(kwargs)
            with self.assertRaises(ValueError):
                tools.sample_geometries(**args)
        with patch.object(tools.nemo.parser, "pega_geom", return_value=(self.geom, self.atoms)), \
                patch.object(tools.nemo.parser, "pega_modos", return_value=self.mode):
            for frequency in (-1., 0., np.nan):
                with patch.object(tools.nemo.parser, "pega_freq",
                                  return_value=(np.array([frequency]), np.array([1.]))):
                    with self.assertRaises(ValueError):
                        tools.sample_geometries("unused", 1, 300)

    def test_real_log_public_api_is_reproducible(self):
        log = str(Path(__file__).parent / "tddft" / "freqs1.log")
        with patch.object(tools, "Parallel", side_effect=lambda **kwargs: Parallel(n_jobs=1)):
            first = tools.sample_geometries(log, 8, 300, seed=42)
            second = tools.sample_geometries(log, 8, 300, seed=42)
            zero = tools.sample_geometries(log, 2, 0, seed=42)
        np.testing.assert_array_equal(first[0], second[0])
        np.testing.assert_array_equal(first[2], second[2])
        self.assertEqual(first[2].shape, (18, 3, 8))
        self.assertTrue(np.all(np.isfinite(zero[2])))
        geom, atoms = tools.nemo.parser.pega_geom(log)
        adj = tools.adjacency(geom, atoms)
        for geometry in first[2].transpose(2, 0, 1):
            self.assertTrue(tools._sampling_geometry_valid(geometry, geom, atoms, adj))


if __name__ == "__main__":
    unittest.main()
