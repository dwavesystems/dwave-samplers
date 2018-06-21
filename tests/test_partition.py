import itertools
import unittest

import dimod
import numpy as np

import savanna


class TestLogPartitionBQM(unittest.TestCase):
    def test_three_path_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction(0, 1, .69)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 0, .5)
        bqm.add_offset(0)

        pos = {0: (0, 0), 1: (1, 0), 2: (0, 1)}

        logZ = savanna.log_partition_bqm(bqm, pos)

        en = list(-bqm.energy(dict(zip(range(len(bqm)), config)))
                  for config in itertools.product((-1, 1), repeat=len(bqm)))

        self.assertAlmostEqual(np.log(np.sum(np.exp(en))), logZ)

    def test_four_path_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction(0, 1, +1.0)
        bqm.add_interaction(1, 2, +1.0)
        bqm.add_interaction(2, 3, +1.0)
        bqm.add_interaction(0, 3, +1.0)
        bqm.add_offset(1.8)

        pos = {0: (+1, +1),
               1: (-1, +1),
               2: (-1, -1),
               3: (+1, -1)}

        logZ = savanna.log_partition_bqm(bqm, pos)

        en = list(-bqm.energy(dict(zip(range(len(bqm)), config)))
                  for config in itertools.product((-1, 1), repeat=len(bqm)))

        self.assertAlmostEqual(np.log(np.sum(np.exp(en))), logZ)

    def test_FrustTriangleL39(self):
        from tests.data import bqm_L39, pos_L39

        logZ = savanna.log_partition_bqm(bqm_L39, pos_L39)
