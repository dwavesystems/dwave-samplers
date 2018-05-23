import unittest

import dimod

import savanna


class TestGroundStateSample(unittest.TestCase):
    def test_one_interaction_equality(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): -1})

        rotation_system = {'a': ['b'], 'b': ['a']}

        sample = savanna.ground_state_sample(bqm, rotation_system)

        self.assertIn(sample, [{'a': +1, 'b': +1}, {'a': -1, 'b': -1}])

    def test_one_interaction_not_equality(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): 1})

        rotation_system = {'a': ['b'], 'b': ['a']}

        sample = savanna.ground_state_sample(bqm, rotation_system)

        self.assertIn(sample, [{'a': -1, 'b': +1}, {'a': +1, 'b': -1}])