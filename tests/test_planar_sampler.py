import unittest

import dimod

from dwave.samplers.planar import PlanarGraphSampler


# from dwave.samplers.planar.sampler_old import ground_state_bqm


class TestGroundStateBQM(unittest.TestCase):
    def test_NAE3SAT_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction('a', 'b', +1.0)
        bqm.add_interaction('b', 'c', +1.0)
        bqm.add_interaction('c', 'a', +1.0)

        pos = {'a': (0, 0), 'b': (1, 0), 'c': (0, 1)}

        sample = PlanarGraphSampler().sample(bqm, pos)

        self.assertEqual(set(sample.values()), {-1, +1})

    def test_grid_15x15(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        for x in range(15):
            for y in range(15):
                bqm.add_interaction((x, y), (x + 1, y), 1)
                bqm.add_interaction((x, y), (x + 1, y + 1), 1)
                bqm.add_interaction((x, y), (x, y + 1), 1)

        def pos(v): return v

        sample = PlanarGraphSampler().sample(bqm, pos)

        self.assertEqual(set(sample.values()), {-1, +1})

    def test_grid_15x15_ferromagnet(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        for x in range(15):
            for y in range(15):
                bqm.add_interaction((x, y), (x + 1, y), -1)
                bqm.add_interaction((x, y), (x + 1, y + 1), -1)
                bqm.add_interaction((x, y), (x, y + 1), -1)

        def pos(v): return v

        sample = PlanarGraphSampler().sample(bqm, pos)

        # should all be the same
        self.assertEqual(len(set(sample.values())), 1)
        self.assertTrue(set(sample.values()).issubset({-1, 1}))
