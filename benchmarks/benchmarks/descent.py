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

import numpy as np
import greedy


class SteepestGradientDescentCython(object):
    params = [1, 1000, 1000000]
    param_names = ['num_samples']

    def setup(self, num_samples):
        self.linear_biases = [2, 2]
        self.coupler_starts = [0]
        self.coupler_ends = [1]
        self.coupler_weights = [-1]
        self.initial_states = np.tile(np.array([1, 1], dtype=np.int8), (num_samples, 1))

    def time_single_flip(self, num_samples):
        samples, energies, info = greedy.descent.steepest_gradient_descent(
            num_samples, self.linear_biases, self.coupler_starts,
            self.coupler_ends, self.coupler_weights, self.initial_states)
