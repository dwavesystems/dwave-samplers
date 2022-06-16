# Copyright 2019 D-Wave Systems Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import itertools
import unittest

import dimod
import networkx as nx
import numpy as np

from dwave.samplers.tree.solve import solve_bqm_wrapper
from dwave.samplers.tree.sample import sample_bqm_wrapper
from dwave.samplers.tree.utilities import elimination_order_width, min_fill_heuristic


class TestWrappers(unittest.TestCase):
    def test_single_interaction_spin(self):
        bqm = dimod.BQM({0: .1, 1: 0}, {(0,1): -1}, 0, 'SPIN')
        _, elimination_order = min_fill_heuristic(bqm)

        # Test solver
        sample, energies = solve_bqm_wrapper(bqm=bqm,
                                              order=elimination_order,
                                              max_complexity=5)
        np.testing.assert_array_equal(sample, [[-1, -1]])
        np.testing.assert_array_equal(energies, [-1.1])

        # Test sampler
        samples, data = sample_bqm_wrapper(bqm=bqm,
                                           beta=3,
                                           max_complexity=5,
                                           order=elimination_order,
                                           marginals=True,
                                           num_reads=3,
                                           seed=222)

        # Check that at least one sample has the ground state
        self.assertIn(np.array([-1,-1]), samples)
        
        # Check that all info is returned
        keys = ['log_partition_function', 'variable_marginals', 'interaction_marginals', 'interactions']
        for key in keys:
            self.assertIn(key, data)

    def test_empty(self):
        bqm = dimod.BQM({},{}, 0, 'SPIN')
        _, elimination_order = min_fill_heuristic(bqm)

        with self.assertRaises(ValueError):
            solve_bqm_wrapper(bqm=bqm,
                              order=elimination_order,
                              max_complexity=5)

        with self.assertRaises(ValueError):
            sample_bqm_wrapper(bqm=bqm,
                               beta=1,
                               max_complexity=2,
                               order=elimination_order)

    def test_invalid_order(self):
        bqm = dimod.BQM({0: 0, 1: 1, 2: -1}, {}, 0, 'BINARY')

        with self.assertRaises(ValueError):
            elimination_order = [0,1]
            solve_bqm_wrapper(bqm=bqm, max_complexity=2, order=elimination_order)

        with self.assertRaises(ValueError):
            elimination_order = ['a','b','c']
            sample_bqm_wrapper(bqm=bqm, beta=1, max_complexity=2, order=elimination_order)


class TestTreewidth(unittest.TestCase):

    def check_order(self, bqm, treewidth, order):
        self.assertEqual(set(bqm.variables), set(order))
        self.assertEqual(treewidth, elimination_order_width(bqm, order))

    def test_clique(self):
        bqm = dimod.generators.gnp_random_bqm(10, 1, 'BINARY')

        tw, order = min_fill_heuristic(bqm)
        self.assertEqual(tw, 9)
        self.check_order(bqm, tw, order)

    def test_cycle(self):
        bqm = dimod.BQM('BINARY')
        for v in range(43):
            bqm.add_quadratic(v, (v + 1) % 43, 1)

        tw, order = min_fill_heuristic(bqm)
        self.assertEqual(tw, 2)
        self.check_order(bqm, tw, order)

    def test_empty(self):
        bqm = dimod.BQM('BINARY')
        tw, order = min_fill_heuristic(bqm)
        self.assertEqual(tw, 0)
        self.assertEqual(order, [])

    def test_exceptions(self):
        bqm = dimod.from_networkx_graph(nx.complete_graph(6), vartype='BINARY')
        order = range(4)

        with self.assertRaises(ValueError):
            elimination_order_width(bqm, order)

        order = range(7)
        with self.assertRaises(ValueError):
            elimination_order_width(bqm, order)

    def test_graphs(self):
        H = nx.complete_graph(2)
        H.add_edge(2, 3)

        graphs = [nx.complete_graph(7),
                  nx.balanced_tree(5, 3),
                  nx.barbell_graph(8, 11),
                  nx.cycle_graph(5),
                  H]

        for G in graphs:
            bqm = dimod.from_networkx_graph(G, vartype='BINARY')
            tw, order = min_fill_heuristic(bqm)
            self.check_order(bqm, tw, order)

    def test_path(self):
        with self.subTest(n=1):
            bqm = dimod.Binary('x')
            tw, order = min_fill_heuristic(bqm)
            self.assertEqual(tw, 1)
            self.assertEqual(order, ['x'])
            self.check_order(bqm, tw, order)

        with self.subTest(n=2):
            bqm = dimod.Binary('x')*dimod.Binary('y')
            tw, order = min_fill_heuristic(bqm)
            self.assertEqual(tw, 1)
            self.assertIn(order, [['x', 'y'], ['y', 'x']])
            self.check_order(bqm, tw, order)

        with self.subTest(n=5):
            bqm = dimod.BQM('BINARY')
            for v in range(4):
                bqm.add_quadratic(v, v+1, 1)
            tw, order = min_fill_heuristic(bqm)
            self.assertEqual(tw, 1)
            self.assertIn(order[0], [0, 4])  # starts at one of the ends
            self.check_order(bqm, tw, order)

        # all possible variable orderings for a 5path
        for combo in itertools.permutations(range(5)):
            with self.subTest(path=combo):
                bqm = dimod.BQM('BINARY')
                bqm.add_variables_from((v, 0) for v in combo)
                for v in range(4):
                    bqm.add_quadratic(v, v+1, 0)
                tw, order = min_fill_heuristic(bqm)
                self.assertEqual(tw, 1)
                self.assertIn(order[0], [0, 4])  # starts at one of the ends
                self.check_order(bqm, tw, order)

        with self.subTest(path='edbde'):
            bqm = dimod.BQM('BINARY')
            bqm.add_variables_from((v, 0) for v in 'edbca')
            bqm.add_quadratic_from({'ab': 1, 'bc': 1, 'cd': 1, 'de': 1})
            tw, order = min_fill_heuristic(bqm)
            self.assertEqual(tw, 1)

            self.check_order(bqm, tw, order)
