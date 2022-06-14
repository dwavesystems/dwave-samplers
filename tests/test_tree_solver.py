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

import inspect
import itertools
import unittest

import dimod

from dwave.samplers.tree import TreeDecompositionSolver


class TestConstruction(unittest.TestCase):
    def test_construction(self):
        sampler = TreeDecompositionSolver()
        dimod.testing.assert_sampler_api(sampler)

        # check that the args exposed by parameters is consistent with the
        # sampler inputs
        args = {arg for arg in inspect.getfullargspec(sampler.sample).args
                if arg != 'self' and arg != 'bqm'}
        self.assertEqual(set(sampler.parameters), args)

        self.assertEqual(sampler.properties, {'max_treewidth': 25})


@dimod.testing.load_sampler_bqm_tests(TreeDecompositionSolver())
class TestSample(unittest.TestCase):
    def test_empty(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        sampleset = TreeDecompositionSolver().sample(bqm)
        dimod.testing.assert_response_energies(sampleset, bqm)

    def test_empty_num_reads(self):
        bqm = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        sampleset = TreeDecompositionSolver().sample(bqm, num_reads=10)
        self.assertEqual(len(sampleset), 10)
        dimod.testing.assert_response_energies(sampleset, bqm)

    def test_consistent_dtype(self):
        bqm_empty = dimod.BinaryQuadraticModel.empty(dimod.BINARY)
        bqm = dimod.BinaryQuadraticModel.from_qubo({(0, 0): -1, (0, 1): 1})

        sampleset_empty = TreeDecompositionSolver().sample(bqm_empty)
        sampleset = TreeDecompositionSolver().sample(bqm)

        self.assertEqual(sampleset_empty.record.sample.dtype,
                         sampleset.record.sample.dtype)
        self.assertEqual(sampleset_empty.record.energy.dtype,
                         sampleset.record.energy.dtype)

    def test_single_variable_spin(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {})

        samples = TreeDecompositionSolver().sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        self.assertEqual(list(samples), [{'a': 1}])
        dimod.testing.assert_response_energies(samples, bqm)

    def test_single_variable_binary(self):
        sampleset = TreeDecompositionSolver().sample_qubo({(0, 0): 1}, num_reads=1)
        self.assertEqual(sampleset.first.sample, {0: 0})

    def test_single_interaction(self):
        sampler = TreeDecompositionSolver()

        dimod.testing.assert_sampler_api(sampler)

        bqm = dimod.BinaryQuadraticModel.from_ising({'a': -1}, {'ab': 1})

        samples = sampler.sample(bqm, num_reads=1)

        self.assertEqual(len(samples), 1)
        self.assertEqual(list(samples), [{'a': 1, 'b': -1}])
        dimod.testing.assert_response_energies(samples, bqm)

    def test_chain(self):

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {(v, v+1): -1 for v in range(99)})

        samples = TreeDecompositionSolver().sample(bqm, num_reads=3)
        dimod.testing.assert_response_energies(samples, bqm)

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

        samples = TreeDecompositionSolver().sample(bqm, num_reads=2)
        dimod.testing.assert_response_energies(samples, bqm)

        self.assertEqual(len(samples), 2)
        ground, excited = samples.samples()

        self.assertEqual(set(ground.values()), {1})
        self.assertEqual(sum(excited.values()), len(excited) - 1)

    def test_offset(self):
        bqm = dimod.BinaryQuadraticModel.from_qubo({pair: -1 for pair in itertools.combinations(range(4), 2)})
        bqm.offset += -5

        samples = TreeDecompositionSolver().sample(bqm, num_reads=2)
        dimod.testing.assert_response_energies(samples, bqm)

    def test_num_reads_gt_max_samples(self):
        bqm = dimod.BinaryQuadraticModel.from_qubo({(0, 0): -1, (0, 1): 1})

        # there are only 4 possible samples for the bqm, say we want 101 reads
        sampleset = TreeDecompositionSolver().sample(bqm, num_reads=101)

        self.assertEqual(sum(sampleset.record.num_occurrences), 101)
