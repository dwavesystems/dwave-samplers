import itertools
import unittest

import dimod
# import numpy as np

from orang import solve_polynomial


class Test_solve_polynomial(unittest.TestCase):

    def check_energies(self, poly, samples):
        for sample, energy in samples.data(['sample', 'energy']):
            self.assertAlmostEqual(poly_energy(poly, sample), energy)

    def test_empty(self):
        samples = solve_polynomial({}, 'SPIN')
        self.assertEqual(len(samples), 0)

    def test_cubed(self):
        samples3 = solve_polynomial({'aaa': -1}, 'SPIN')
        samples1 = solve_polynomial({'a': -1}, 'SPIN')
        samples13 = solve_polynomial({'aaa': -.5, 'a': -.5}, 'SPIN')

        self.check_energies({'a': -1}, samples3)
        self.check_energies({'a': -1}, samples1)
        self.check_energies({'a': -1}, samples13)

    def test_one_interaction(self):
        poly = {('a', 'b'): -1, ('a',): .2}

        samples = solve_polynomial(poly, dimod.SPIN, num_reads=5)
        self.check_energies(poly, samples)
        self.assertEqual(len(samples), 4)

    def test_3interaction(self):
        poly = {('a', 'b', 'c'): 1}
        samples = solve_polynomial(poly, dimod.SPIN, num_reads=100)
        self.check_energies(poly, samples)
        self.assertEqual(len(samples), 8)

    def test_offset(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({0: .5}, {(v, v+1): -1 for v in range(99)}, offset=.15)

        poly = bqm_to_polynomial(bqm)

        samples = solve_polynomial(poly, dimod.SPIN, num_reads=100)
        self.check_energies(poly, samples)
        dimod.testing.assert_response_energies(samples, bqm)

    def test_chain_binary(self):

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {(v, v+1): -1 for v in range(99)})

        poly = bqm_to_polynomial(bqm.binary)

        samples = solve_polynomial(poly, 'BINARY', 3)
        self.check_energies(poly, samples)

        self.assertEqual(len(samples), 3)
        ground0, ground1, excited = samples.samples()

        # the two ground states should be all 0 or all 1
        self.assertEqual(len(set(ground0.values())), 1)
        self.assertEqual(len(set(ground1.values())), 1)
        self.assertEqual(set(ground1.values()).union(ground0.values()), {0, 1})

        # first excited should have one frustrated edge
        self.assertEqual(sum(excited[u] != excited[v] for u, v in bqm.quadratic), 1)

    def test_chain_spin(self):

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {(v, v+1): -1 for v in range(9)})

        poly = bqm_to_polynomial(bqm)

        samples = solve_polynomial(poly, 'SPIN', 3)
        self.check_energies(poly, samples)

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

        poly = bqm_to_polynomial(bqm)

        samples = solve_polynomial(poly, bqm.vartype, 2)
        self.check_energies(poly, samples)

        self.assertEqual(len(samples), 2)
        ground, excited = samples.samples()

        self.assertEqual(set(ground.values()), {1})
        self.assertEqual(sum(excited.values()), len(excited) - 1)


def poly_energy(poly, sample):
    energy = 0
    for interaction, bias in poly.items():
        mul = 1
        for v in interaction:
            mul *= sample[v]
        energy += bias*mul
    return energy


def bqm_to_polynomial(bqm):
    poly = {(v,): bias for v, bias in bqm.linear.items() if bias}
    poly.update(bqm.quadratic)

    if bqm.offset:
        poly[tuple()] = bqm.offset

    return poly
