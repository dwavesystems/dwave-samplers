import itertools
import unittest
import random

import dimod
import numpy as np

from dwave.samplers.planar import log_partition_bqm


class TestLogPartitionBQM(unittest.TestCase):
    def test_three_path_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction(0, 1, .69)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 0, .5)

        pos = {0: (0, 0), 1: (1, 0), 2: (0, 1)}

        logZ = log_partition_bqm(bqm, pos)

        en = list(-bqm.energy(dict(zip(range(len(bqm)), config)))
                  for config in itertools.product((-1, 1), repeat=len(bqm)))

        self.assertAlmostEqual(np.log(np.sum(np.exp(en))), logZ)

    def test_four_path_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction(0, 1, +1.0)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 3, +1.0)
        bqm.add_interaction(0, 3, +1.0)
        bqm.offset += 1.8

        pos = {0: (+1, +1),
               1: (-1, +1),
               2: (-1, -1),
               3: (+1, -1)}

        logZ = log_partition_bqm(bqm, pos)

        en = list(-bqm.energy(dict(zip(range(len(bqm)), config)))
                  for config in itertools.product((-1, 1), repeat=len(bqm)))

        self.assertAlmostEqual(np.log(np.sum(np.exp(en))), logZ)

    # def test_FrustTriangleL39(self):
    #     from tests.data import bqm_L39, pos_L39
    #
    #     _ = log_partition_bqm(bqm_L39, pos_L39)

    def test_square_with_chord(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({0: -0.0, 1: -0.0, 2: -0.0, 3: -0.0},
                                                    {(1, 2): 2, (0, 1): 1, (1, 3): 1, (2, 3): 1, (0, 2): 1})
        pos = {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)}

        logZ = log_partition_bqm(bqm, pos)

        en = []
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(range(len(bqm)), config))
            en.append(bqm.energy(sample))

        self.assertAlmostEqual(np.log(np.sum(np.exp(-1*np.asarray(en)))), logZ)

    def test_square_with_chord_2(self):
        bqm = dimod.BinaryQuadraticModel({0: -0.0, 1: -0.0, 2: -0.0, 3: -0.0},
                                         {(1, 2): 200, (0, 1): 100, (1, 3): 100, (2, 3): 100, (0, 2): 100},
                                         -0.0, dimod.SPIN)

        pos = {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)}

        logZ = log_partition_bqm(bqm, pos)

        en = []
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(range(len(bqm)), config))
            en.append(bqm.energy(sample))

        self.assertAlmostEqual(np.log(np.sum(np.exp(-1*np.asarray(en)))), logZ)

    def test_functional_square_with_chord_random(self):
        nodes = list(range(4))
        edges = [[1, 2], [0, 1], [1, 3], [2, 3], [0, 2]]
        pos = {0: (0, 0), 1: (0, 1), 2: (1, 0), 3: (1, 1)}
        for __ in range(10):
            bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

            random.shuffle(nodes)
            for v in nodes:
                bqm.add_variable(v, 0)

            random.shuffle(edges)
            for edge in edges:
                random.shuffle(edge)
            for u, v in edges:
                bqm.add_interaction(u, v, random.uniform(-1, 1))

            logZ = log_partition_bqm(bqm, pos)

            en = list(-bqm.energy(dict(zip(range(len(bqm)), config)))
                      for config in itertools.product((-1, 1), repeat=len(bqm)))

            self.assertAlmostEqual(np.log(np.sum(np.exp(en))), logZ)
