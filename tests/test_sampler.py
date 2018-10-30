import itertools
import unittest

import dimod
import numpy as np

from orang import OrangSampler


class TestOrangSampler(unittest.TestCase):
    def test_empty(self):
        sampler = OrangSampler()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        samples = sampler.sample(bqm)

        self.assertEqual(len(samples), 0)
        dimod.testing.assert_response_energies(samples, bqm)

    def test_single_variable(self):
        sampler = OrangSampler()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {})

        samples = sampler.sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        dimod.testing.assert_response_energies(samples, bqm)

    def test_single_interaction(self):
        sampler = OrangSampler()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {'ab': 1})

        samples = sampler.sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        dimod.testing.assert_response_energies(samples, bqm)

    def test_larger_problem(self):
        sampler = OrangSampler()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {'ab': 1, 'bc': -1, 'cd': +1})

        samples = sampler.sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        dimod.testing.assert_response_energies(samples, bqm)

    # def test_logz_three_path_bqm(self):
    #     bqm = dimod.BinaryQuadraticModel.empty(dimod.BINARY)

    #     bqm.add_interaction(0, 1, .69)
    #     bqm.add_interaction(1, 2, +1.0)
    #     bqm.add_interaction(2, 0, .5)
    #     bqm.add_offset(0)

    #     b = 1

    #     samples = OrangSampler().sample(bqm, num_reads=1, beta=b)

    #     logZ = samples.info['log_partition_function']

    #     en = []
    #     for config in itertools.product((-1, 1), repeat=len(bqm)):
    #         sample = dict(zip(range(len(bqm)), config))
    #         en.append(bqm.energy(sample))

    #     print(en)

    #     print(bqm.binary.offset + np.log(2))

    #     print(np.log(np.sum(np.exp(-b*np.asarray(en)))), logZ)

    #     self.assertAlmostEqual(np.log(np.sum(np.exp(-b*np.asarray(en)))), logZ)
