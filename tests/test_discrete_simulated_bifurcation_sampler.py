# Copyright 2025 D-Wave
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

import unittest
import inspect
import warnings
from itertools import product

import numpy as np
import dimod

from dwave.samplers.sbm import DiscreteSimulatedBifurcationSampler


class TestTimingInfo(unittest.TestCase):
    def setUp(self) -> None:
        empty = dimod.BQM(dimod.SPIN)
        one = dimod.BQM.from_ising({"a": 1}, {})
        two = dimod.BQM.from_ising({}, {("abc", (1, 2)): -1})

        sampler = DiscreteSimulatedBifurcationSampler()
        rng = np.random.default_rng(48448418563)

        self.sample_sets = []
        for bqm in [empty, one, two]:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                sample_set = sampler.sample(bqm, seed=int(rng.integers(2**30)))
            self.sample_sets.append(sample_set)

        self.timing_keys = {"preprocessing_ns", "postprocessing_ns", "sampling_ns"}

    def test_keys_exist(self):
        for sample_set in self.sample_sets:
            with self.subTest(ss=sample_set):
                self.assertTrue(self.timing_keys.issubset(sample_set.info["timing"]))

    def test_strictly_positive_timings(self):
        for sample_set in self.sample_sets:
            for category, duration in sample_set.info["timing"].items():
                self.assertGreater(duration, 0)


class TestSchedules(unittest.TestCase):
    def test_default_schedule(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        num_vars = 40
        h = {v: 0 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 10
        num_sweeps = 5

        resp = sampler.sample_ising(h, J, num_reads=num_reads, num_sweeps=num_sweeps)
        row, col = resp.record.sample.shape

        self.assertEqual(row, num_reads)
        self.assertEqual(col, num_vars)
        self.assertIs(resp.vartype, dimod.SPIN)

    def test_custom_a_schedule(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        num_vars = 8
        h = {v: 0 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        a_schedule = np.linspace(0, 1, 7)

        resp = sampler.sample_ising(
            h, J, num_reads=5, num_sweeps=len(a_schedule), a_schedule=a_schedule
        )
        self.assertEqual(resp.record.sample.shape, (5, num_vars))

    def test_single_sweep_schedule(self):
        # num_sweeps == 1 falls back to a constant schedule of length 1
        sampler = DiscreteSimulatedBifurcationSampler()
        J = {("a", "b"): -1, ("b", "c"): -1}
        resp = sampler.sample_ising({}, J, num_reads=3, num_sweeps=1)
        self.assertEqual(resp.record.sample.shape, (3, 3))


class TestDiscreteSimulatedBifurcationSampler(unittest.TestCase):

    def test_instantiation(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        dimod.testing.assert_sampler_api(sampler)

    def test_good_kwargs(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        kwargs = dict(inspect.signature(sampler.sample).parameters)
        kwargs.pop("bqm")
        kwargs.pop("kwargs")
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            kwargs_out = sampler.remove_unknown_kwargs(**kwargs)
        self.assertEqual(kwargs.keys(), kwargs_out.keys(), "Keyword arguments removed")

    def test_bad_kwargs(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        kwargs = {"foobar": None}
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            kwargs_out = sampler.remove_unknown_kwargs(**kwargs)
        self.assertFalse(kwargs_out, "Keyword arguments not removed")

    def test_c0(self):
        # c0 = 0.5*sqrt((N-1)/sum(J^2))
        J0 = np.random.random() - 1
        J = {("a", "b"): J0}
        bqm = dimod.BinaryQuadraticModel({}, J, 0, dimod.SPIN)
        response = DiscreteSimulatedBifurcationSampler().sample(bqm)
        c0 = response.info["c0"]
        self.assertAlmostEqual(c0, 0.5 / abs(J0))

        J1 = np.random.random() - 1
        J[("a", "c")] = J1
        bqm = dimod.BinaryQuadraticModel({}, J, 0, dimod.SPIN)
        response = DiscreteSimulatedBifurcationSampler().sample(bqm)
        c0 = response.info["c0"]
        self.assertAlmostEqual(c0, 0.5 * np.sqrt(2 / (J0**2 + J1**2)))

    def test_info_keys(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        response = sampler.sample_ising({}, {("a", "b"): -1}, num_reads=4)
        for key in ("c0", "x", "y", "timing"):
            self.assertIn(key, response.info)
        self.assertEqual(response.info["x"].shape, (4, 2))
        self.assertEqual(response.info["y"].shape, (4, 2))

    def test_sample_ising(self):
        h = {"a": 0, "b": 0}
        J = {("a", "b"): -1}

        resp = DiscreteSimulatedBifurcationSampler().sample_ising(h, J)

        row, col = resp.record.sample.shape

        self.assertEqual(col, 2)  # should get back two variables
        self.assertIs(resp.vartype, dimod.SPIN)  # should be ising

    def test_sample_qubo(self):
        # A QUBO maps to an Ising model with non-zero linear bias, which SBM
        # ignores (with a warning); we only check the returned vartype here.
        Q = {(0, 1): 1}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            resp = DiscreteSimulatedBifurcationSampler().sample_qubo(Q)

        row, col = resp.record.sample.shape

        self.assertEqual(col, 2)  # should get back two variables
        self.assertIs(resp.vartype, dimod.BINARY)  # should be qubo

    def test_basic_response(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        h = {"a": 0, "b": 0}
        J = {("a", "b"): -1}
        response = sampler.sample_ising(h, J)

        self.assertIsInstance(
            response, dimod.SampleSet, "Sampler returned an unexpected response type"
        )

    def test_energies_match_samples(self):
        # The returned energies should equal the bqm energies of the samples.
        sampler = DiscreteSimulatedBifurcationSampler()
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {("a", "b"): -1, ("b", "c"): 1, ("a", "c"): -1}
        )
        response = sampler.sample(bqm, num_reads=10)
        self.assertTrue(np.allclose(response.record.energy, bqm.energies(response)))

    def test_ferromagnet_ground_state(self):
        # Fully connected ferromagnet: ground state is all-aligned, energy -num_edges.
        sampler = DiscreteSimulatedBifurcationSampler()
        num_vars = 8
        h = {v: 0 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        response = sampler.sample_ising(h, J, num_reads=20, num_sweeps=100, seed=42)
        num_edges = len(J)
        self.assertEqual(response.first.energy, -num_edges)

    def test_num_reads(self):
        sampler = DiscreteSimulatedBifurcationSampler()

        h = {}
        J = {("a", "b"): 0.5, (0, "a"): -1, (1, "b"): 0.0}

        for num_reads in (1, 10, 100, 3223):
            response = sampler.sample_ising(h, J, num_reads=num_reads)
            row, col = response.record.sample.shape

            self.assertEqual(row, num_reads)
            self.assertEqual(col, 4)

    def test_default_num_reads(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        response = sampler.sample_ising({}, {("a", "b"): -1})
        self.assertEqual(len(response), 1)

    def test_empty_problem(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        empty = dimod.BQM(dimod.SPIN)
        response = sampler.sample(empty)
        self.assertEqual(response.record.sample.shape[1], 0)

    def test_h_ignored_warning(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            sampler.sample_ising({"a": 1, "b": -1}, {("a", "b"): -1})
        self.assertTrue(
            any("non-zero h" in str(w.message) for w in caught),
            "Expected a warning about ignored non-zero h values",
        )

    def test_no_warning_zero_h(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            sampler.sample_ising({"a": 0, "b": 0}, {("a", "b"): -1})
        self.assertFalse(
            any("non-zero h" in str(w.message) for w in caught),
            "Unexpected warning for zero h values",
        )

    def test_seed(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        num_vars = 4
        h = {v: 0 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 3

        previous_x = None
        for num_sweeps, seed in product((0, 3), (1, 2)):
            response0 = sampler.sample_ising(
                h, J, num_reads=num_reads, num_sweeps=num_sweeps, seed=seed
            )
            response1 = sampler.sample_ising(
                h, J, num_reads=num_reads, num_sweeps=num_sweeps, seed=seed
            )

            samples0 = response0.record.sample
            states_x = response0.info["x"]
            samples1 = response1.record.sample
            states_x1 = response1.info["x"]
            self.assertTrue(
                np.array_equal(samples0, samples1),
                "Same seed returned different results",
            )
            self.assertTrue(
                np.array_equal(states_x, states_x1),
                "Same seed returned different x results",
            )
            if seed == 2 and num_sweeps == 0:
                self.assertFalse(
                    np.array_equal(states_x, previous_x),
                    "Different seed returned same results",
                )
            else:
                previous_x = states_x

    def test_disconnected_problem(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        h = {}
        J = {
            # K_3
            (0, 1): -1,
            (1, 2): -1,
            (0, 2): -1,
            # disconnected K_3
            (3, 4): -1,
            (4, 5): -1,
            (3, 5): -1,
        }

        resp = sampler.sample_ising(h, J, num_sweeps=100, num_reads=100)

        row, col = resp.record.sample.shape

        self.assertEqual(row, 100)
        self.assertEqual(col, 6)
        self.assertIs(resp.vartype, dimod.SPIN)

    def test_interrupt_error(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        num_vars = 40
        h = {v: 0 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 100

        def f():
            raise NotImplementedError

        resp = sampler.sample_ising(
            h, J, num_reads=num_reads, num_sweeps=50, interrupt_function=f
        )

        self.assertEqual(len(resp), 1)

    def test_interrupt_stops_early(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        J = {("a", "b"): -1, ("b", "c"): -1}

        def stop():
            return True

        resp = sampler.sample_ising(
            {}, J, num_reads=10, num_sweeps=50, interrupt_function=stop
        )
        self.assertEqual(len(resp), 1)

    def test_interrupt_not_callable(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with self.assertRaises(TypeError):
            sampler.sample_ising({}, {("a", "b"): -1}, interrupt_function=5)


class TestInitialStates(unittest.TestCase):
    def test_initial_x_and_y(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        initial_x = 1 - 2 * np.random.random(size=(4, 2))
        initial_y = 1 - 2 * np.random.random(size=(4, 2))
        resp = sampler.sample_ising(
            {}, {("a", "b"): -1},
            initial_x=initial_x, initial_y=initial_y,
            num_reads=4, num_sweeps=10,
        )
        self.assertEqual(resp.record.sample.shape, (4, 2))

    def test_initial_x_without_y(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with self.assertRaises(ValueError):
            sampler.sample_ising(
                {}, {("a", "b"): -1}, initial_x=np.ones((2, 2)), num_reads=2
            )

    def test_initial_y_without_x(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with self.assertRaises(ValueError):
            sampler.sample_ising(
                {}, {("a", "b"): -1}, initial_y=np.ones((2, 2)), num_reads=2
            )

    def test_initial_shape_mismatch(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with self.assertRaises(ValueError):
            sampler.sample_ising(
                {}, {("a", "b"): -1},
                initial_x=np.ones((2, 2)), initial_y=np.ones((3, 2)),
                num_reads=2,
            )

    def test_initial_dtype_mismatch(self):
        sampler = DiscreteSimulatedBifurcationSampler()
        with self.assertRaises(ValueError):
            sampler.sample_ising(
                {}, {("a", "b"): -1},
                initial_x=np.ones((2, 2), dtype=np.float64),
                initial_y=np.ones((2, 2), dtype=np.float32),
                num_reads=2,
            )

    def test_initial_scale(self):
        # initial_scale bounds the magnitude of generated initial states.
        sampler = DiscreteSimulatedBifurcationSampler()
        resp = sampler.sample_ising(
            {}, {("a", "b"): -1}, num_reads=100, num_sweeps=0,
            initial_scale=0.5, seed=1,
        )
        self.assertLessEqual(np.max(np.abs(resp.info["x"])), 0.5)
        self.assertLessEqual(np.max(np.abs(resp.info["y"])), 0.5)

    def test_num_sweeps_zero_uses_sign(self):
        # With num_sweeps == 0 the samples are the signs of initial_x.
        sampler = DiscreteSimulatedBifurcationSampler()
        initial_x = np.array([[0.5, -0.5], [-0.5, 0.5]])
        initial_y = np.zeros((2, 2))
        resp = sampler.sample_ising(
            {}, {("a", "b"): -1},
            initial_x=initial_x, initial_y=initial_y,
            num_reads=2, num_sweeps=0,
        )
        variables = list(resp.variables)
        sample = resp.record.sample
        expected = np.sign(initial_x)
        # reorder expected columns to match returned variable order
        order = [["a", "b"].index(v) for v in variables]
        self.assertTrue(np.array_equal(sample, expected[:, order]))


if __name__ == "__main__":
    unittest.main()
