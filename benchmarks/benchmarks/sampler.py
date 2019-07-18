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

import dimod
import greedy


class SteepestDescentSampler(object):

    def setup(self):
        self.sampler = greedy.SteepestDescentSampler()

    def time_single_flip_1_read(self):
        self.sampler.sample_ising({0: 2, 1: 2}, {(0, 1): -1})

    def time_single_flip_1k_reads(self):
        self.sampler.sample_ising({0: 2, 1: 2}, {(0, 1): -1}, num_reads=1000)

    def time_single_flip_1M_reads(self):
        self.sampler.sample_ising({0: 2, 1: 2}, {(0, 1): -1}, num_reads=1000000)
