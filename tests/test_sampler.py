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

import numpy as np
import dimod

from greedy.sampler import SteepestDescentSampler


class TestSteepestDescentSampler(unittest.TestCase):

    def setUp(self):
        self.convex_bqm = dimod.BQM.from_ising({0: 2, 1: 2}, {(0, 1): -1})

    def test_instantiation(self):
        sampler = SteepestDescentSampler()
        dimod.testing.assert_sampler_api(sampler)

    def test_small_convex_ising(self):
        # a convex section of hyperbolic paraboloid in the Ising space,
        # with global minimum at (-1,-1)
        h = {0: 2, 1: 2}
        J = {(0, 1): -1}

        ss = SteepestDescentSampler().sample_ising(h, J)

        self.assertEqual(len(ss), 1)
        self.assertEqual(len(ss.variables), 2)
        self.assertEqual(ss.record.sample.shape, (1, 2))
        np.testing.assert_array_equal(ss.record.sample[0], [-1, -1])

    def test_small_convex_qubo(self):
        # a convex section of hyperbolic paraboloid in the QUBO space,
        # with global minimum at (0,0)
        Q = {(0, 0): 2, (1, 1): 2, (0, 1): -1}

        ss = SteepestDescentSampler().sample_qubo(Q)

        self.assertEqual(len(ss), 1)
        self.assertEqual(len(ss.variables), 2)
        self.assertEqual(ss.record.sample.shape, (1, 2))
        np.testing.assert_array_equal(ss.record.sample[0], [0, 0])
