import itertools
import unittest

import dimod

from orang import OrangSolver


class TestOrangSolver(unittest.TestCase):
    def test_empty(self):
        sampler = OrangSolver()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        samples = sampler.sample(bqm)

        self.assertEqual(len(samples), 0)
        dimod.testing.assert_response_energies(samples, bqm)

    def test_single_variable(self):
        sampler = OrangSolver()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {})

        samples = sampler.sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        self.assertEqual(list(samples), [{'a': 1}])
        dimod.testing.assert_response_energies(samples, bqm)

    def test_single_interaction(self):
        sampler = OrangSolver()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {'ab': 1})

        samples = sampler.sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        self.assertEqual(list(samples), [{'a': 1, 'b': -1}])
        dimod.testing.assert_response_energies(samples, bqm)

    def test_chain(self):

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {(v, v+1): -1 for v in range(99)})

        samples = OrangSolver().sample(bqm, num_reads=3)
        dimod.testing.assert_response_energies(samples, bqm)

        self.assertEqual(len(samples), 3)
        ground0, ground1, excited = samples.samples()

        # the two ground states should be all -1 or all 1
        self.assertEqual(len(set(ground0.values())), 1)
        self.assertEqual(len(set(ground1.values())), 1)
        self.assertEqual(set(ground1.values()).union(ground0.values()), {-1, 1})

        # first excited should have one frustrated edge
        self.assertEqual(sum(excited[u] != excited[v] for u, v in bqm.quadratic), 1)

    def test_clique(self):
        bqm = dimod.BinaryQuadraticModel.from_qubo({pair: -1 for pair in itertools.combinations(range(20), 2)})

        samples = OrangSolver().sample(bqm, num_reads=2)
        dimod.testing.assert_response_energies(samples, bqm)

        self.assertEqual(len(samples), 2)
        ground, excited = samples.samples()

        self.assertEqual(set(ground.values()), {1})
        self.assertEqual(sum(excited.values()), len(excited) - 1)
