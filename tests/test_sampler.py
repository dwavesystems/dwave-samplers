import unittest

import dimod

import savanna


class TestGroundStateBQM(unittest.TestCase):
    def test_NAE3SAT_bqm(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        bqm.add_interaction('a', 'b', +1.0)
        bqm.add_interaction('b', 'c', +1.0)
        bqm.add_interaction('c', 'a', +1.0)

        pos = {'a': (0, 0), 'b': (1, 0), 'c': (0, 1)}

        sample = savanna.ground_state_bqm(bqm, pos)

        self.assertEqual(set(sample.values()), {-1, +1})
