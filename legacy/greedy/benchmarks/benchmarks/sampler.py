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

import networkx as nx

import dimod
import greedy


class SteepestDescentSimple(object):
    params = [1, 1000, 1000000]
    param_names = ['num_reads']
    repeat = (3, 10, 60)

    def setup(self, num_reads):
        self.sampler = greedy.SteepestDescentSampler()
        self.h = {0: 2, 1: 2}
        self.J = {(0, 1): -1}

    def time_single_flip(self, num_reads):
        self.sampler.sample_ising(self.h, self.J, num_reads=num_reads, seed=0)


class SteepestDescentComplete(object):
    """Test steepest descent on complete graphs."""

    params = ([100, 1000, 2000], [1, 10])
    param_names = ['graph_size', 'num_reads']
    repeat = (3, 10, 60)
    timeout = 300

    def setup(self, graph_size, num_reads):
        self.sampler = greedy.SteepestDescentSampler()
        self.bqm = dimod.generators.random.ran_r(r=1, graph=graph_size, seed=0)

    def time_ran1(self, graph_size, num_reads):
        self.sampler.sample(self.bqm, num_reads=num_reads, seed=0)


class SteepestDescentSparse(object):
    """Test steepest descent on Erdős-Rényi sparse graphs with varying density."""

    params = ([100, 1000, 2000], [0.05, 0.1, 0.25, 0.5], [1, 10])
    param_names = ['graph_size', 'graph_density', 'num_reads']
    repeat = (3, 10, 60)
    timeout = 300

    def setup(self, graph_size, graph_density, num_reads):
        self.sampler = greedy.SteepestDescentSampler()
        self.graph = nx.fast_gnp_random_graph(n=graph_size, p=graph_density, seed=0)
        self.bqm = dimod.generators.random.ran_r(r=1, graph=self.graph, seed=0)

    def time_ran1(self, graph_size, graph_density, num_reads):
        self.sampler.sample(self.bqm, num_reads=num_reads, seed=0)


class SteepestDescentLargeSparse(object):
    """Test steepest descent on large and sparse Erdős-Rényi graphs.

    Note:
        Because of practical limitations on BQM size in `dimod`, we limit
        the number of edges in test problems below to ~5M.
    """

    params = ([2000, 10000, 15000], [0.01, 0.05], [1, 10])
    param_names = ['graph_size', 'graph_density', 'num_reads']
    repeat = (3, 10, 60)
    timeout = 300

    def setup(self, graph_size, graph_density, num_reads):
        self.sampler = greedy.SteepestDescentSampler()
        self.graph = nx.fast_gnp_random_graph(n=graph_size, p=graph_density, seed=0)
        self.bqm = dimod.generators.random.ran_r(r=1, graph=self.graph, seed=0)

    def time_ran1(self, graph_size, graph_density, num_reads):
        self.sampler.sample(self.bqm, num_reads=num_reads, seed=0, large_sparse_opt=True)


class SteepestDescentCompleteOverBQMTypes(object):
    """Test steepest descent on complete graphs."""

    params = ([dimod.AdjArrayBQM,
               dimod.AdjDictBQM,
               dimod.AdjMapBQM,
               dimod.AdjVectorBQM], [1, 10])
    param_names = ['bqm_type', 'num_reads']
    repeat = (3, 10, 60)
    timeout = 300
    graph_size = 1000

    def setup(self, bqm_type, num_reads):
        self.sampler = greedy.SteepestDescentSampler()
        self.bqm = dimod.generators.random.ran_r(
            r=1, graph=self.graph_size, cls=bqm_type, seed=0)

    def time_ran1(self, bqm_type, num_reads):
        self.sampler.sample(self.bqm, num_reads=num_reads, seed=0)
