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

from dwave.samplers.greedy.descent import steepest_gradient_descent


class SteepestGradientDescentCython(unittest.TestCase):
    large_sparse_opt = False

    def test_steepest_gradient_descent_on_small_convex_ising(self):
        """The sampler must converge to a global minimum of a convex Ising problem."""

        # a convex section of hyperbolic paraboloid in the Ising space,
        # with global minimum at (-1,-1):
        #
        #   h = {0: 2, 1: 2}
        #   J = {(0, 1): -1}

        num_samples = 1
        linear_biases = [2, 2]
        coupler_starts, coupler_ends, coupler_weights = [0], [1], [-1]
        initial_states = np.array([[1, 1]], dtype=np.int8)

        samples, energies, num_steps = steepest_gradient_descent(
            num_samples, linear_biases, coupler_starts, coupler_ends,
            coupler_weights, initial_states, self.large_sparse_opt)

        self.assertEqual(samples.shape, (1, 2))
        np.testing.assert_array_equal(samples, [[-1, -1]])
        np.testing.assert_array_equal(energies, [-5])


class SteepestGradientDescentLargeSparseCython(SteepestGradientDescentCython):
    large_sparse_opt = True
