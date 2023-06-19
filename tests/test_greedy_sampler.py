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
from parameterized import parameterized, parameterized_class

import dimod
from dimod.testing.sampler import BQM_SUBCLASSES

from dwave.samplers.greedy.sampler import SteepestDescentSampler


BQM_CLASSES = [(cls, ) for cls in BQM_SUBCLASSES]


@parameterized_class([
    # vary large & sparse problem optimization via sampling params
    {"params": dict(large_sparse_opt=False)},
    {"params": dict(large_sparse_opt=True)},
])
@dimod.testing.load_sampler_bqm_tests(SteepestDescentSampler)
class TestSteepestDescentSampler(unittest.TestCase):

    def test_instantiation(self):
        """The sampler must conform to `dimod.Sampler` interface."""

        sampler = SteepestDescentSampler()
        dimod.testing.assert_sampler_api(sampler)

    @parameterized.expand(BQM_CLASSES)
    def test_edges(self, BQM):
        """The sampler correctly handles edge cases."""

        sampler = SteepestDescentSampler()

        # empty bqm
        ss = sampler.sample(BQM.empty('SPIN'), **self.params)
        self.assertEqual(ss.first.sample, {})

        # single-variable problem
        ss = sampler.sample(BQM.from_ising({'x': 1}, {}), **self.params)
        self.assertEqual(ss.first.sample, {'x': -1})

    @parameterized.expand(BQM_CLASSES)
    def test_validation(self, BQM):
        """Inputs are validated."""

        sampler = SteepestDescentSampler()
        empty = BQM.from_ising({}, {})
        bqm = BQM.from_ising({'x': 1}, {})

        with self.assertRaises(TypeError):
            sampler.sample(bqm, num_reads=2.3, **self.params)

        with self.assertRaises(ValueError):
            sampler.sample(bqm, num_reads=0, **self.params)

        with self.assertRaises(ValueError):
            sampler.sample(bqm, initial_states=(), **self.params)

        with self.assertRaises(TypeError):
            sampler.sample(bqm, seed=2.3, **self.params)

        with self.assertRaises(ValueError):
            sampler.sample(bqm, seed=-1, **self.params)

        with self.assertRaises(ValueError):
            sampler.sample(bqm, initial_states_generator='invalid', **self.params)

        init = dimod.SampleSet.from_samples({'y': 1}, vartype='SPIN', energy=0)
        with self.assertRaises(ValueError):
            sampler.sample(bqm, initial_states=init, **self.params)

        init = dimod.SampleSet.from_samples({0: 1}, vartype='SPIN', energy=0)
        with self.assertRaises(ValueError):
            sampler.sample(empty, initial_states=init, **self.params)

    def test_small_convex_ising(self):
        """The sampler must converge to a global minimum of a convex Ising problem."""

        # a convex section of hyperbolic paraboloid in the Ising space,
        # with global minimum at (-1,-1)
        h = {0: 2, 1: 2}
        J = {(0, 1): -1}

        ss = SteepestDescentSampler().sample_ising(h, J, **self.params)

        self.assertEqual(len(ss), 1)
        self.assertEqual(len(ss.variables), 2)
        self.assertEqual(ss.record.sample.shape, (1, 2))
        np.testing.assert_array_equal(ss.record.sample[0], [-1, -1])

    def test_small_convex_qubo(self):
        """The sampler must converge to a global minimum of a convex QUBO problem."""

        # a convex section of hyperbolic paraboloid in the QUBO space,
        # with global minimum at (0,0)
        Q = {(0, 0): 2, (1, 1): 2, (0, 1): -1}

        ss = SteepestDescentSampler().sample_qubo(Q, **self.params)

        self.assertEqual(len(ss), 1)
        self.assertEqual(len(ss.variables), 2)
        self.assertEqual(ss.record.sample.shape, (1, 2))
        np.testing.assert_array_equal(ss.record.sample[0], [0, 0])

    @parameterized.expand(BQM_CLASSES)
    def test_variable_relabeling(self, BQM):
        """The sampler must accept non-integer/sequential variable names."""

        # use a small convex bqm
        bqm = BQM.from_ising({'x': 2, 'y': 2}, {'xy': -1})

        ss = SteepestDescentSampler().sample(bqm, **self.params)

        self.assertSetEqual(set(ss.variables), set(bqm.variables))
        self.assertDictEqual(ss.first.sample, {'x': -1, 'y': -1})

    @parameterized.expand(BQM_CLASSES)
    def test_unsortable_labels(self, BQM):
        """The sampler must accept unorderable variable names."""

        # use a small convex bqm
        bqm = BQM.from_ising({0: 2, 'a': 2}, {(0, 'a'): -1})

        ss = SteepestDescentSampler().sample(bqm, **self.params)

        self.assertSetEqual(set(ss.variables), set(bqm.variables))
        self.assertEqual(ss.first.energy, -5)

    @parameterized.expand(BQM_CLASSES)
    def test_reproducible_convergence(self, BQM):
        """On a convex problem, all samples must correspond to the global minimum."""

        # a convex section of hyperbolic paraboloid in the Ising space,
        # with global minimum at (-1,-1)
        bqm = BQM.from_ising({0: 2, 1: 2}, {(0, 1): -1})
        num_samples = 100

        # each sample is derived from a random initial state
        ss = SteepestDescentSampler().sample(
            bqm, num_reads=num_samples, **self.params)

        expected = np.tile([-1, -1], (num_samples, 1))

        np.testing.assert_array_equal(ss.record.sample, expected)

    @parameterized.expand(BQM_CLASSES)
    def test_initial_states_randomization(self, BQM):
        """Assuming uniform RNG, samples must follow a bimodal distribution."""

        # use a simple centrally symmetric hyperbolic paraboloid with
        # two minima in Ising space: (-1, 1) and (1, -1)
        bqm = BQM.from_ising({}, {'xy': 1})

        num = 1000
        tol = 0.10

        ss = SteepestDescentSampler().sample(
            bqm, num_reads=num, **self.params
        ).aggregate()

        # sanity check: two minima
        self.assertEqual(len(ss), 2)

        # check uniform distribution, with some tolerance `tol`
        inf = num / len(ss) - num * tol
        sup = num / len(ss) + num * tol
        for record in ss.record:
            self.assertTrue(inf < record.num_occurrences < sup)

    @parameterized.expand(BQM_CLASSES)
    def test_initial_states(self, BQM):
        """The sampler must deterministically converge from the initial state(s)."""

        # use a simple centrally symmetric hyperbolic paraboloid with
        # two minima in Ising space: (-1, 1) and (1, -1)
        bqm = BQM.from_ising({}, {(0, 1): 1})

        initial_states = dimod.SampleSet.from_samples(
            {0: -1, 1: -1}, vartype='SPIN', energy=0)

        # move along 0-dimension from (-1, -1) and settle in a local minimum (1, -1)
        ss = SteepestDescentSampler().sample(
            bqm, initial_states=initial_states, **self.params)

        self.assertEqual(len(ss), 1)
        self.assertDictEqual(ss.first.sample, {0: 1, 1: -1})

        # repeat this for 1000 samples, with initial state tiled
        num_reads = 1000
        ss = SteepestDescentSampler().sample(
            bqm, initial_states=initial_states,
            initial_states_generator='tile', num_reads=num_reads, **self.params
        ).aggregate()

        self.assertEqual(len(ss), 1)
        self.assertDictEqual(ss.first.sample, {0: 1, 1: -1})
        self.assertEqual(ss.first.num_occurrences, num_reads)

    @parameterized.expand(BQM_CLASSES)
    def test_initial_states_generation_and_validation(self, BQM):
        """Initial states are properly validated/expanded with the state generator."""

        sampler = SteepestDescentSampler()
        bqm = BQM.from_ising({'x': 1}, {})

        # num_reads inferred from initial_states
        init = dimod.SampleSet.from_samples([{'x': 1}, {'x': -1}],
                                            vartype='SPIN', energy=0)
        ss = sampler.sample(bqm, initial_states=init, **self.params)
        self.assertEqual(len(ss), len(init))

        # reads truncated with `num_reads`
        ss = sampler.sample(bqm, initial_states=init, num_reads=1, **self.params)
        self.assertEqual(len(ss), 1)

        # init states tiled according to `num_reads`
        init = dimod.SampleSet.from_samples({'x': 1}, vartype='SPIN', energy=0)
        ss = sampler.sample(
            bqm, num_reads=10, initial_states=init,
            initial_states_generator='tile', **self.params)
        self.assertEqual(len(ss), 10)
        self.assertEqual(list(ss.aggregate().samples()), [{'x': -1}])

        # tiling fails
        with self.assertRaises(ValueError):
            sampler.sample(
                bqm, initial_states=None, initial_states_generator='tile',
                **self.params)

        # tiling truncates
        init = dimod.SampleSet.from_samples([{'x': 1}, {'x': -1}],
                                            vartype='SPIN', energy=0)
        ss = sampler.sample(
            bqm, num_reads=1, initial_states=init,
            initial_states_generator='tile', **self.params)
        self.assertEqual(len(ss), 1)

        # none generator works
        init = dimod.SampleSet.from_samples([{'x': 1}, {'x': -1}],
                                            vartype='SPIN', energy=0)
        ss = sampler.sample(
            bqm, num_reads=1, initial_states=init,
            initial_states_generator='none', **self.params)
        self.assertEqual(len(ss), 1)

        # none generator fails
        with self.assertRaises(ValueError):
            sampler.sample(
                bqm, num_reads=1, initial_states=None,
                initial_states_generator='none', **self.params)

    @parameterized.expand([
        (([-1, -1], 'ab'), ),
        ((np.array([-1, -1]), 'ab'), ),
        ((np.array([-1, -1], dtype=np.int8), 'ab'), ),
    ])
    def test_initial_states_sample_like(self, initial_states):
        """Samples-like is accepted for initial_states."""

        # global minimum at (-1,-1)
        bqm = dimod.BQM.from_ising({'a': 2, 'b': 2}, {'ab': -1})

        ss = SteepestDescentSampler().sample(bqm, initial_states=initial_states)

        # result is identical to initial state, with zero downhill moves
        np.testing.assert_array_equal(ss.record.sample, [[-1, -1]])
        np.testing.assert_array_equal(ss.record.num_steps, [0])


class TestTimingInfo(unittest.TestCase):
    def setUp(self) -> None:
        empty = dimod.BQM(dimod.SPIN)
        one = dimod.BQM.from_ising({"a": 1}, {})
        two = dimod.BQM.from_ising({}, {("abc", (1, 2)): -1})

        sampler = SteepestDescentSampler()
        rng = np.random.default_rng(59921)

        self.sample_sets = []
        for bqm in [empty, one, two]:
            sample_set = sampler.sample(bqm, seed=rng.integers(2**30))
            self.sample_sets.append(sample_set)

        self.timing_keys = {"preprocessing_ns", "postprocessing_ns", "sampling_ns"}

    def test_keys_exist(self):
        for sample_set in self.sample_sets:
            with self.subTest(ss=sample_set):
                self.timing_keys.issubset(sample_set.info['timing'])

    def test_strictly_postive_timings(self):
        for sample_set in self.sample_sets:
            for category, duration in sample_set.info['timing'].items():
                self.assertGreater(duration, 0)


if __name__ == "__main__":
    unittest.main()