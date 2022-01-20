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

import unittest

import dimod
import numpy as np
import dwave_networkx as dnx

from orang.solve import solve_bqm_wrapper
from orang.sample import sample_bqm_wrapper

class TestWrappers(unittest.TestCase):
    def test_single_interaction_spin(self):
        bqm = dimod.BQM({0: .1, 1: 0}, {(0,1): -1}, 0, 'SPIN')
        _, elimination_order = dnx.min_fill_heuristic(bqm.adj)

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
        _, elimination_order = dnx.min_fill_heuristic(bqm.adj)

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
