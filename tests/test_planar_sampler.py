# Copyright 2022 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS F ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest

import dimod
import numpy as np

from dwave.samplers.planar import PlanarGraphSolver
from dwave.samplers.planar import PlanarGraphSampler


class TestGroundStateBQM(unittest.TestCase):
    def test_energy(self):
        bqm = dimod.BinaryQuadraticModel({}, {(0, 1): 1, (1, 2): 1}, 0, "SPIN")
        ss = PlanarGraphSolver().sample(bqm)
        np.testing.assert_array_equal(ss.record.energy, bqm.energies(ss))

    def test_noPosProvided_threeVariables(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction('a', 'b', +1.0)
        bqm.add_interaction('b', 'c', +1.0)
        bqm.add_interaction('c', 'a', +1.0)

        sample = PlanarGraphSolver().sample(bqm)
        self.assertDictEqual(sample.first.sample, {'a': 1, 'b': -1, 'c': -1})

    def test_noPosProvided_tenVariables(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        for i in range(0, 10):
            bqm.add_interaction(i, i+1, 1)

        sample = PlanarGraphSolver().sample(bqm)
        self.assertDictEqual(sample.first.sample,
                             {0: 1, 1: -1, 2: 1, 3: -1, 4: 1, 5: -1, 6: 1, 7: -1, 8: 1, 9: -1, 10: 1}
                             )

    def test_NAE3SAT_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction('a', 'b', +1.0)
        bqm.add_interaction('b', 'c', +1.0)
        bqm.add_interaction('c', 'a', +1.0)

        pos = {'a': (0, 0), 'b': (1, 0), 'c': (0, 1)}

        sample = PlanarGraphSolver().sample(bqm, pos)
        self.assertDictEqual(sample.first.sample, {'a': 1, 'b': -1, 'c': -1})

    def test_grid_15x15(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        for x in range(15):
            for y in range(15):
                bqm.add_interaction((x, y), (x + 1, y), 1)
                bqm.add_interaction((x, y), (x + 1, y + 1), 1)
                bqm.add_interaction((x, y), (x, y + 1), 1)

        sample = PlanarGraphSolver().sample(bqm)

        self.assertEqual(set(sample.first.sample.values()), {-1, +1})
        self.assertDictEqual(
            sample.first.sample,
            {(0, 0): 1, (0, 1): -1, (0, 2): 1, (0, 3): -1, (0, 4): 1, (0, 5): -1, (0, 6): 1, (0, 7): -1, (0, 8): 1,
             (0, 9): -1, (0, 10): 1, (0, 11): -1, (0, 12): 1, (0, 13): -1, (0, 14): 1, (0, 15): -1, (1, 0): -1,
             (1, 1): -1, (1, 2): 1, (1, 3): -1, (1, 4): 1, (1, 5): -1, (1, 6): 1, (1, 7): -1, (1, 8): 1, (1, 9): -1,
             (1, 10): 1, (1, 11): -1, (1, 12): 1, (1, 13): -1, (1, 14): 1, (1, 15): -1, (2, 0): 1, (2, 1): 1, (2, 2): 1,
             (2, 3): -1, (2, 4): 1, (2, 5): -1, (2, 6): 1, (2, 7): -1, (2, 8): 1, (2, 9): -1, (2, 10): 1, (2, 11): -1,
             (2, 12): 1, (2, 13): -1, (2, 14): 1, (2, 15): -1, (3, 0): -1, (3, 1): -1, (3, 2): -1, (3, 3): -1,
             (3, 4): 1, (3, 5): -1, (3, 6): 1, (3, 7): -1, (3, 8): 1, (3, 9): -1, (3, 10): 1, (3, 11): -1, (3, 12): 1,
             (3, 13): -1, (3, 14): 1, (3, 15): -1, (4, 0): 1, (4, 1): 1, (4, 2): 1, (4, 3): 1, (4, 4): 1, (4, 5): -1,
             (4, 6): 1, (4, 7): -1, (4, 8): 1, (4, 9): -1, (4, 10): 1, (4, 11): -1, (4, 12): 1, (4, 13): -1, (4, 14): 1,
             (4, 15): -1, (5, 0): -1, (5, 1): -1, (5, 2): -1, (5, 3): -1, (5, 4): -1, (5, 5): -1, (5, 6): 1, (5, 7): -1,
             (5, 8): 1, (5, 9): -1, (5, 10): 1, (5, 11): -1, (5, 12): 1, (5, 13): -1, (5, 14): 1, (5, 15): -1,
             (6, 0): 1, (6, 1): 1, (6, 2): 1, (6, 3): 1, (6, 4): 1, (6, 5): 1, (6, 6): 1, (6, 7): -1, (6, 8): 1,
             (6, 9): -1, (6, 10): 1, (6, 11): -1, (6, 12): 1, (6, 13): -1, (6, 14): 1, (6, 15): -1, (7, 0): -1,
             (7, 1): -1, (7, 2): -1, (7, 3): -1, (7, 4): -1, (7, 5): -1, (7, 6): -1, (7, 7): -1, (7, 8): 1, (7, 9): -1,
             (7, 10): 1, (7, 11): -1, (7, 12): 1, (7, 13): -1, (7, 14): 1, (7, 15): -1, (8, 0): 1, (8, 1): 1, (8, 2): 1,
             (8, 3): 1, (8, 4): 1, (8, 5): 1, (8, 6): 1, (8, 7): 1, (8, 8): 1, (8, 9): -1, (8, 10): 1, (8, 11): -1,
             (8, 12): 1, (8, 13): -1, (8, 14): 1, (8, 15): -1, (9, 0): -1, (9, 1): -1, (9, 2): -1, (9, 3): -1,
             (9, 4): -1, (9, 5): -1, (9, 6): -1, (9, 7): -1, (9, 8): -1, (9, 9): -1, (9, 10): 1, (9, 11): -1,
             (9, 12): 1, (9, 13): -1, (9, 14): 1, (9, 15): -1, (10, 0): 1, (10, 1): 1, (10, 2): 1, (10, 3): 1,
             (10, 4): 1, (10, 5): 1, (10, 6): 1, (10, 7): 1, (10, 8): 1, (10, 9): 1, (10, 10): 1, (10, 11): -1,
             (10, 12): 1, (10, 13): -1, (10, 14): 1, (10, 15): -1, (11, 0): -1, (11, 1): -1, (11, 2): -1, (11, 3): -1,
             (11, 4): -1, (11, 5): -1, (11, 6): -1, (11, 7): -1, (11, 8): -1, (11, 9): -1, (11, 10): -1, (11, 11): -1,
             (11, 12): 1, (11, 13): -1, (11, 14): 1, (11, 15): -1, (12, 0): 1, (12, 1): 1, (12, 2): 1, (12, 3): 1,
             (12, 4): 1, (12, 5): 1, (12, 6): 1, (12, 7): 1, (12, 8): 1, (12, 9): 1, (12, 10): 1, (12, 11): 1,
             (12, 12): 1, (12, 13): -1, (12, 14): 1, (12, 15): -1, (13, 0): -1, (13, 1): -1, (13, 2): -1, (13, 3): -1,
             (13, 4): -1, (13, 5): -1, (13, 6): -1, (13, 7): -1, (13, 8): -1, (13, 9): -1, (13, 10): -1, (13, 11): -1,
             (13, 12): -1, (13, 13): -1, (13, 14): 1, (13, 15): -1, (14, 0): 1, (14, 1): 1, (14, 2): 1, (14, 3): 1,
             (14, 4): 1, (14, 5): 1, (14, 6): 1, (14, 7): 1, (14, 8): 1, (14, 9): 1, (14, 10): 1, (14, 11): 1,
             (14, 12): 1, (14, 13): 1, (14, 14): 1, (14, 15): -1, (15, 0): -1, (15, 1): -1, (15, 2): -1, (15, 3): -1,
             (15, 4): -1, (15, 5): -1, (15, 6): -1, (15, 7): -1, (15, 8): -1, (15, 9): -1, (15, 10): -1, (15, 11): -1,
             (15, 12): -1, (15, 13): -1, (15, 14): -1, (15, 15): -1}
        )

    def test_grid_15x15_ferromagnet(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        for x in range(15):
            for y in range(15):
                bqm.add_interaction((x, y), (x + 1, y), -1)
                bqm.add_interaction((x, y), (x + 1, y + 1), -1)
                bqm.add_interaction((x, y), (x, y + 1), -1)

        sample = PlanarGraphSolver().sample(bqm)

        # should all be the same
        self.assertEqual(set(sample.first.sample.values()), {-1})
        self.assertDictEqual(
            sample.first.sample,
            {(0, 0): -1, (0, 1): -1, (0, 2): -1, (0, 3): -1, (0, 4): -1, (0, 5): -1, (0, 6): -1, (0, 7): -1, (0, 8): -1,
             (0, 9): -1, (0, 10): -1, (0, 11): -1, (0, 12): -1, (0, 13): -1, (0, 14): -1, (0, 15): -1, (1, 0): -1,
             (1, 1): -1, (1, 2): -1, (1, 3): -1, (1, 4): -1, (1, 5): -1, (1, 6): -1, (1, 7): -1, (1, 8): -1, (1, 9): -1,
             (1, 10): -1, (1, 11): -1, (1, 12): -1, (1, 13): -1, (1, 14): -1, (1, 15): -1, (2, 0): -1, (2, 1): -1,
             (2, 2): -1, (2, 3): -1, (2, 4): -1, (2, 5): -1, (2, 6): -1, (2, 7): -1, (2, 8): -1, (2, 9): -1,
             (2, 10): -1, (2, 11): -1, (2, 12): -1, (2, 13): -1, (2, 14): -1, (2, 15): -1, (3, 0): -1, (3, 1): -1,
             (3, 2): -1, (3, 3): -1, (3, 4): -1, (3, 5): -1, (3, 6): -1, (3, 7): -1, (3, 8): -1, (3, 9): -1,
             (3, 10): -1, (3, 11): -1, (3, 12): -1, (3, 13): -1, (3, 14): -1, (3, 15): -1, (4, 0): -1, (4, 1): -1,
             (4, 2): -1, (4, 3): -1, (4, 4): -1, (4, 5): -1, (4, 6): -1, (4, 7): -1, (4, 8): -1, (4, 9): -1,
             (4, 10): -1, (4, 11): -1, (4, 12): -1, (4, 13): -1, (4, 14): -1, (4, 15): -1, (5, 0): -1, (5, 1): -1,
             (5, 2): -1, (5, 3): -1, (5, 4): -1, (5, 5): -1, (5, 6): -1, (5, 7): -1, (5, 8): -1, (5, 9): -1,
             (5, 10): -1, (5, 11): -1, (5, 12): -1, (5, 13): -1, (5, 14): -1, (5, 15): -1, (6, 0): -1, (6, 1): -1,
             (6, 2): -1, (6, 3): -1, (6, 4): -1, (6, 5): -1, (6, 6): -1, (6, 7): -1, (6, 8): -1, (6, 9): -1,
             (6, 10): -1, (6, 11): -1, (6, 12): -1, (6, 13): -1, (6, 14): -1, (6, 15): -1, (7, 0): -1, (7, 1): -1,
             (7, 2): -1, (7, 3): -1, (7, 4): -1, (7, 5): -1, (7, 6): -1, (7, 7): -1, (7, 8): -1, (7, 9): -1,
             (7, 10): -1, (7, 11): -1, (7, 12): -1, (7, 13): -1, (7, 14): -1, (7, 15): -1, (8, 0): -1, (8, 1): -1,
             (8, 2): -1, (8, 3): -1, (8, 4): -1, (8, 5): -1, (8, 6): -1, (8, 7): -1, (8, 8): -1, (8, 9): -1,
             (8, 10): -1, (8, 11): -1, (8, 12): -1, (8, 13): -1, (8, 14): -1, (8, 15): -1, (9, 0): -1, (9, 1): -1,
             (9, 2): -1, (9, 3): -1, (9, 4): -1, (9, 5): -1, (9, 6): -1, (9, 7): -1, (9, 8): -1, (9, 9): -1,
             (9, 10): -1, (9, 11): -1, (9, 12): -1, (9, 13): -1, (9, 14): -1, (9, 15): -1, (10, 0): -1, (10, 1): -1,
             (10, 2): -1, (10, 3): -1, (10, 4): -1, (10, 5): -1, (10, 6): -1, (10, 7): -1, (10, 8): -1, (10, 9): -1,
             (10, 10): -1, (10, 11): -1, (10, 12): -1, (10, 13): -1, (10, 14): -1, (10, 15): -1, (11, 0): -1,
             (11, 1): -1, (11, 2): -1, (11, 3): -1, (11, 4): -1, (11, 5): -1, (11, 6): -1, (11, 7): -1, (11, 8): -1,
             (11, 9): -1, (11, 10): -1, (11, 11): -1, (11, 12): -1, (11, 13): -1, (11, 14): -1, (11, 15): -1,
             (12, 0): -1, (12, 1): -1, (12, 2): -1, (12, 3): -1, (12, 4): -1, (12, 5): -1, (12, 6): -1, (12, 7): -1,
             (12, 8): -1, (12, 9): -1, (12, 10): -1, (12, 11): -1, (12, 12): -1, (12, 13): -1, (12, 14): -1,
             (12, 15): -1, (13, 0): -1, (13, 1): -1, (13, 2): -1, (13, 3): -1, (13, 4): -1, (13, 5): -1, (13, 6): -1,
             (13, 7): -1, (13, 8): -1, (13, 9): -1, (13, 10): -1, (13, 11): -1, (13, 12): -1, (13, 13): -1,
             (13, 14): -1, (13, 15): -1, (14, 0): -1, (14, 1): -1, (14, 2): -1, (14, 3): -1, (14, 4): -1, (14, 5): -1,
             (14, 6): -1, (14, 7): -1, (14, 8): -1, (14, 9): -1, (14, 10): -1, (14, 11): -1, (14, 12): -1, (14, 13): -1,
             (14, 14): -1, (14, 15): -1, (15, 0): -1, (15, 1): -1, (15, 2): -1, (15, 3): -1, (15, 4): -1, (15, 5): -1,
             (15, 6): -1, (15, 7): -1, (15, 8): -1, (15, 9): -1, (15, 10): -1, (15, 11): -1, (15, 12): -1, (15, 13): -1,
             (15, 14): -1, (15, 15): -1}
        )


class TestLogPartitionFunction(unittest.TestCase):
    def _brute_force_logZ(self, bqm, beta=1.0):
        import itertools
        en = [-beta * bqm.energy(dict(zip(bqm.variables, config)))
              for config in itertools.product((-1, 1), repeat=len(bqm))]
        return np.log(np.sum(np.exp(en)))

    def test_three_path_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
        bqm.add_interaction(0, 1, .69)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 0, .5)
        bqm.add_offset(0)

        sampleset = PlanarGraphSampler().sample(bqm)
        logZ = sampleset.info['log_partition_function']

        self.assertAlmostEqual(logZ, self._brute_force_logZ(bqm))

    def test_four_path_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
        bqm.add_interaction(0, 1, +1.0)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 3, +1.0)
        bqm.add_interaction(0, 3, +1.0)
        bqm.add_offset(1.8)

        sampleset = PlanarGraphSampler().sample(bqm)
        logZ = sampleset.info['log_partition_function']

        self.assertAlmostEqual(logZ, self._brute_force_logZ(bqm))

    def test_square_with_chord(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {0: -0.0, 1: -0.0, 2: -0.0, 3: -0.0},
            {(1, 2): 2, (0, 1): 1, (1, 3): 1, (2, 3): 1, (0, 2): 1})

        sampleset = PlanarGraphSampler().sample(bqm)
        logZ = sampleset.info['log_partition_function']

        self.assertAlmostEqual(logZ, self._brute_force_logZ(bqm))

    def test_beta(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
        bqm.add_interaction(0, 1, .69)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 0, .5)

        for beta in (0.1, 0.5, 2.0):
            sampleset = PlanarGraphSampler().sample(bqm, beta=beta)
            logZ = sampleset.info['log_partition_function']
            self.assertAlmostEqual(logZ, self._brute_force_logZ(bqm, beta=beta))

    def test_log_partition_in_info(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)
        bqm.add_interaction('a', 'b', +1.0)
        bqm.add_interaction('b', 'c', +1.0)
        bqm.add_interaction('c', 'a', +1.0)

        sampleset = PlanarGraphSampler().sample(bqm)

        self.assertIn('log_partition_function', sampleset.info)

    def test_disconnected_log_partition(self):
        # two independent triangles; each component has its own flip symmetry
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): -1.0, (0, 2): 0.5,
                 (3, 4): 0.8, (4, 5): 1.2, (3, 5): -0.6})

        sampleset = PlanarGraphSampler().sample(bqm, num_reads=1, seed=0)

        self.assertAlmostEqual(sampleset.info['log_partition_function'],
                               self._brute_force_logZ(bqm))

    def test_square_lattice_sampling(self):
        # contracting a square lattice produces intermediate graphs that are
        # planar but not biconnected; sampling must handle these via the
        # spanning-tree contraction and block decomposition
        L = 4
        edges = {((i, j), (i, j + 1)) for i in range(L) for j in range(L - 1)}
        edges |= {((i, j), (i + 1, j)) for i in range(L - 1) for j in range(L)}
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {e: -1.0 for e in edges})

        sampleset = PlanarGraphSampler().sample(bqm, num_reads=3, seed=0)

        self.assertEqual(sampleset.record.sample.shape, (3, L * L))
        self.assertAlmostEqual(sampleset.info['log_partition_function'],
                               self._brute_force_logZ(bqm))

    def test_too_few_variables(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {(0, 1): 1})

        with self.assertRaises(ValueError):
            PlanarGraphSampler().sample(bqm)


class TestSampling(unittest.TestCase):
    def _exact_probabilities(self, bqm, beta=1.0):
        import itertools
        variables = list(bqm.variables)
        spin = bqm.change_vartype(dimod.SPIN, inplace=False)
        weights = {config: np.exp(-beta * spin.energy(dict(zip(variables, config))))
                   for config in itertools.product((-1, 1), repeat=len(variables))}
        Z = sum(weights.values())
        return variables, {c: w / Z for c, w in weights.items()}

    def _max_abs_error(self, bqm, beta=1.0, num_reads=2000, seed=42):
        from collections import Counter
        variables, exact = self._exact_probabilities(bqm, beta=beta)
        sampleset = PlanarGraphSampler().sample(bqm, beta=beta,
                                                num_reads=num_reads, seed=seed)
        order = list(sampleset.variables)
        idx = [order.index(v) for v in variables]

        def to_spin(x):
            return x if bqm.vartype is dimod.SPIN else 2 * x - 1

        counts = Counter(tuple(to_spin(int(row[i])) for i in idx)
                         for row in sampleset.record.sample)
        return max(abs(counts.get(c, 0) / num_reads - p) for c, p in exact.items())

    def test_triangle_distribution(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0})
        self.assertLess(self._max_abs_error(bqm), 0.05)

    def test_square_with_chord_distribution(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(1, 2): 2, (0, 1): 1, (1, 3): 1, (2, 3): 1, (0, 2): 1})
        self.assertLess(self._max_abs_error(bqm), 0.05)

    def test_asymmetric_couplings_distribution(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 0.7, (1, 2): -1.3, (0, 2): 0.4, (2, 3): 0.9, (0, 3): -0.5})
        self.assertLess(self._max_abs_error(bqm), 0.05)

    def test_beta_distribution(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 0.7, (1, 2): -1.3, (0, 2): 0.4, (2, 3): 0.9, (0, 3): -0.5})
        self.assertLess(self._max_abs_error(bqm, beta=2.0), 0.05)

    def test_disconnected_distribution(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): -1.0, (0, 2): 0.5,
                 (3, 4): 0.8, (4, 5): 1.2, (3, 5): -0.6})
        self.assertLess(self._max_abs_error(bqm), 0.05)

    def test_binary_vartype(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0}
        ).change_vartype(dimod.BINARY, inplace=False)
        sampleset = PlanarGraphSampler().sample(bqm, num_reads=10, seed=0)
        self.assertEqual(sampleset.vartype, dimod.BINARY)
        self.assertTrue(set(np.unique(sampleset.record.sample)).issubset({0, 1}))
        self.assertLess(self._max_abs_error(bqm), 0.05)

    def test_num_reads_and_energies(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0})
        sampleset = PlanarGraphSampler().sample(bqm, num_reads=7, seed=0)
        self.assertEqual(sampleset.record.sample.shape, (7, 3))
        np.testing.assert_allclose(
            sampleset.record.energy,
            bqm.energies((sampleset.record.sample, list(sampleset.variables))))


class TestMarginals(unittest.TestCase):
    def _brute_correlations(self, bqm, beta=1.0):
        import itertools
        variables = list(bqm.variables)
        spin = bqm.change_vartype(dimod.SPIN, inplace=False)
        weights = {config: np.exp(-beta * spin.energy(dict(zip(variables, config))))
                   for config in itertools.product((-1, 1), repeat=len(variables))}
        Z = sum(weights.values())

        def to_native(s):
            return s if bqm.vartype is dimod.SPIN else (s + 1) // 2

        correlations = {}
        for u, v in bqm.quadratic:
            c = 0.0
            for su in (-1, 1):
                for sv in (-1, 1):
                    p = sum(w / Z for config, w in weights.items()
                            if dict(zip(variables, config))[u] == su
                            and dict(zip(variables, config))[v] == sv)
                    c += to_native(su) * to_native(sv) * p
            correlations[(u, v)] = c
        return correlations

    def test_correlations_spin(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(1, 2): 2, (0, 1): 1, (1, 3): 1, (2, 3): 1, (0, 2): 1})
        sampleset = PlanarGraphSampler().sample(bqm, beta=0.8, num_reads=1, seed=0)
        exact = self._brute_correlations(bqm, beta=0.8)
        for edge, c in sampleset.info['correlations'].items():
            self.assertAlmostEqual(c, exact[edge])

    def test_correlations_binary(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): -0.7, (0, 2): 0.4}
        ).change_vartype(dimod.BINARY, inplace=False)
        sampleset = PlanarGraphSampler().sample(bqm, num_reads=1, seed=0)
        exact = self._brute_correlations(bqm)
        for edge, c in sampleset.info['correlations'].items():
            self.assertAlmostEqual(c, exact[edge])

    def test_marginals_false_omits_them(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {(0, 1): 1.0, (1, 2): 1.0, (0, 2): 1.0})
        sampleset = PlanarGraphSampler().sample(bqm, num_reads=1, seed=0,
                                                marginals=False)
        self.assertNotIn('correlations', sampleset.info)


