import unittest
import itertools

import dimod
import networkx as nx

import savanna


class TestAgreementIsingEnergy(unittest.TestCase):
    def test_empty(self):
        G = nx.Graph()

        self.assertEqual(savanna.agreement_energy({}, G), 0.0)

    def test_path_2(self):
        G = nx.Graph()
        G.add_edge(0, 1, weight=1.0)

        self.assertEqual(savanna.agreement_energy({0: 0, 1: 0}, G), 0.0)
        self.assertEqual(savanna.agreement_energy({0: 1, 1: 1}, G), 0.0)
        self.assertEqual(savanna.agreement_energy({0: 0, 1: 1}, G), 1.0)
        self.assertEqual(savanna.agreement_energy({0: 1, 1: 0}, G), 1.0)


class TestBQMToAgreementIsing(unittest.TestCase):
    def test_functional_empty(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.BINARY)
        G, offset = savanna.bqm_to_agreement_graph(bqm)
        self.assertEqual(bqm.energy({}),
                         savanna.agreement_energy({}, G, offset))

    def test_functional_one_eq_interaction(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): -1})
        G, offset = savanna.bqm_to_agreement_graph(bqm)

        variables = list(bqm)
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(variables, config))

            bqm_en = bqm.energy(sample)
            graph_en = savanna.agreement_energy(sample, G, offset)

            self.assertEqual(bqm_en, graph_en, '{}, bqm: {}, graph: {}'.format(sample, bqm_en, graph_en))

    def test_functional_one_ne_interaction(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): +1})
        G, offset = savanna.bqm_to_agreement_graph(bqm)

        variables = list(bqm)
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(variables, config))

            bqm_en = bqm.energy(sample)
            graph_en = savanna.agreement_energy(sample, G, offset)

            self.assertEqual(bqm_en, graph_en, '{}, bqm: {}, graph: {}'.format(sample, bqm_en, graph_en))

    def test_functional_larger(self):

        J = {(u, v): u * v for u, v in itertools.combinations(range(-3, 3), 2)}

        bqm = dimod.BinaryQuadraticModel.from_ising({}, J, offset=1.3)
        G, offset = savanna.bqm_to_agreement_graph(bqm)

        variables = list(bqm)
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(variables, config))

            bqm_en = bqm.energy(sample)
            graph_en = savanna.agreement_energy(sample, G, offset)

            self.assertAlmostEqual(bqm_en, graph_en, msg='{}, bqm: {}, graph: {}'.format(sample, bqm_en, graph_en))
