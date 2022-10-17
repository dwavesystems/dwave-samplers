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

import itertools
import unittest

import dimod
import networkx as nx

from dwave.samplers.planar.util import bqm_to_multigraph


def _agreement_energy(state, G, offset=0.0):
    en = offset
    for u, v, data in G.edges(data=True):
        if state[u] != state[v]:
            en += data['weight']
    return en


class TestAgreementIsingEnergy(unittest.TestCase):
    def test_empty(self):
        G = nx.Graph()

        self.assertEqual(_agreement_energy({}, G), 0.0)

    def test_path_2(self):
        G = nx.Graph()
        G.add_edge(0, 1, weight=1.0)

        self.assertEqual(_agreement_energy({0: 0, 1: 0}, G), 0.0)
        self.assertEqual(_agreement_energy({0: 1, 1: 1}, G), 0.0)
        self.assertEqual(_agreement_energy({0: 0, 1: 1}, G), 1.0)
        self.assertEqual(_agreement_energy({0: 1, 1: 0}, G), 1.0)


class TestBQMToAgreementIsing(unittest.TestCase):
    def test_functional_empty(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.BINARY)
        G, offset = bqm_to_multigraph(bqm)
        self.assertEqual(bqm.energy({}),
                         _agreement_energy({}, G, offset))

    def test_functional_one_eq_interaction(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): -1})
        G, offset = bqm_to_multigraph(bqm)

        variables = bqm.variables
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(variables, config))

            bqm_en = bqm.energy(sample)
            graph_en = _agreement_energy(sample, G, offset)

            self.assertEqual(bqm_en, graph_en, '{}, bqm: {}, graph: {}'.format(sample, bqm_en, graph_en))

    def test_functional_one_ne_interaction(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {('a', 'b'): +1})
        G, offset = bqm_to_multigraph(bqm)

        variables = bqm.variables
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(variables, config))

            bqm_en = bqm.energy(sample)
            graph_en = _agreement_energy(sample, G, offset)

            self.assertEqual(bqm_en, graph_en, '{}, bqm: {}, graph: {}'.format(sample, bqm_en, graph_en))

    def test_functional_larger(self):

        J = {(u, v): u * v for u, v in itertools.combinations(range(-3, 3), 2)}

        bqm = dimod.BinaryQuadraticModel.from_ising({}, J, offset=1.3)
        G, offset = bqm_to_multigraph(bqm)

        variables = bqm.variables
        for config in itertools.product((-1, 1), repeat=len(bqm)):
            sample = dict(zip(variables, config))

            bqm_en = bqm.energy(sample)
            graph_en = _agreement_energy(sample, G, offset)

            self.assertAlmostEqual(bqm_en, graph_en, msg='{}, bqm: {}, graph: {}'.format(sample, bqm_en, graph_en))
