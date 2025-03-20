# Copyright 2024 D-Wave
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
import numpy as np
import copy
import itertools
import warnings
import json

import dimod

from dwave.samplers.sqa.sampler import RotorModelAnnealingSampler


class TestTimingInfo(unittest.TestCase):
    def setUp(self):
        empty = dimod.BQM(dimod.SPIN)
        one = dimod.BQM.from_ising({"a": 1}, {})
        two = dimod.BQM.from_ising({}, {("abc", (1, 2)): -1})

        sampler = RotorModelAnnealingSampler()
        rng = np.random.default_rng(48448418563)

        self.sample_sets = []
        for bqm in [empty, one, two]:
            sample_set = sampler.sample(
                bqm, seed=rng.integers(2**30), beta_range=[0.1, 1]
            )
            self.sample_sets.append(sample_set)

        self.timing_keys = {"preprocessing_ns", "postprocessing_ns", "sampling_ns"}

    def test_keys_exist(self):
        for sample_set in self.sample_sets:
            self.assertIn("timing", sample_set.info)
            with self.subTest(ss=sample_set):
                self.timing_keys.issubset(sample_set.info["timing"])

    def test_strictly_postive_timings(self):
        for sample_set in self.sample_sets:
            for category, duration in sample_set.info["timing"].items():
                self.assertGreater(duration, 0)


class TestSchedules(unittest.TestCase):

    def test_schedules(self):
        sampler = RotorModelAnnealingSampler()
        num_vars = 40
        h = {v: -1 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 10
        for schedule_type in ["geometric", "linear"]:
            resp = sampler.sample_ising(
                h, J, num_reads=num_reads, beta_schedule_type=schedule_type
            )

            row, col = resp.record.sample.shape

            self.assertEqual(row, num_reads)
            self.assertEqual(col, num_vars)  # should get back two variables
            self.assertIs(resp.vartype, dimod.SPIN)  # should be ising
            with self.assertRaises(ValueError):
                # Should not accept schedule:
                resp = sampler.sample_ising(
                    h,
                    J,
                    num_reads=num_reads,
                    beta_schedule_type=schedule_type,
                    Hp_field=[-1, 1],
                )
        with self.assertRaises(ValueError):
            sampler.sample_ising(h, J, num_reads=num_reads, beta_schedule_type="asd")

    def test_custom_schedule(self):
        sampler = RotorModelAnnealingSampler()
        num_vars = 40
        h = {v: -1 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 1
        with self.assertRaises(ValueError):
            # Missing schedule
            sampler.sample_ising(h, J, num_reads=num_reads, beta_schedule_type="custom")
        with self.assertRaises(ValueError):
            # Positivity
            sampler.sample_ising(
                h, J, num_reads=num_reads, beta_schedule_type="custom", Hp_field=[-1, 1]
            )
        with self.assertRaises(ValueError):
            # Numeric
            sampler.sample_ising(
                h,
                J,
                num_reads=num_reads,
                beta_schedule_type="custom",
                Hp_field=["asd", 1],
            )
        with self.assertRaises(ValueError):
            # Mismatch
            sampler.sample_ising(
                h,
                J,
                num_reads=num_reads,
                beta_schedule_type="custom",
                Hp_field=[0, 1],
                Hd_field=[0],
            )

        sampler.sample_ising(
            h, J, num_reads=num_reads, beta_schedule_type="custom", Hp_field=[0.1, 1]
        )
        sampler.sample_ising(
            h,
            J,
            num_reads=num_reads,
            beta_schedule_type="custom",
            Hp_field=[0.1, 1],
            Hd_field=[0.1, 0.1],
        )


class TestRotorModelAnnealingSampler(unittest.TestCase):

    def test_instantiation(self):
        sampler = RotorModelAnnealingSampler()
        dimod.testing.assert_sampler_api(sampler)

    def test_one_node_beta_range(self):
        h = {"a": -1}
        bqm = dimod.BinaryQuadraticModel(h, {}, 0, dimod.SPIN)
        response = RotorModelAnnealingSampler().sample(bqm)
        hot_beta, cold_beta = response.info["beta_range"]

        # beta is proportional to 1/temperature, therefore hot_beta < cold_beta
        self.assertLess(hot_beta, cold_beta)
        self.assertNotEqual(
            hot_beta, float("inf"), "Starting value of 'beta_range' is infinite"
        )
        self.assertNotEqual(
            cold_beta, float("inf"), "Final value of 'beta_range' is infinite"
        )

    def test_one_edge_beta_range(self):
        J = {("a", "b"): 1}
        bqm = dimod.BinaryQuadraticModel({}, J, 0, dimod.BINARY)
        response = RotorModelAnnealingSampler().sample(bqm)
        hot_beta, cold_beta = response.info["beta_range"]

        # beta is proportional to 1/temperature, therefore hot_beta < cold_beta
        self.assertLess(hot_beta, cold_beta)
        self.assertNotEqual(
            hot_beta, float("inf"), "Starting value of 'beta_range' is infinite"
        )
        self.assertNotEqual(
            cold_beta, float("inf"), "Final value of 'beta_range' is infinite"
        )

    def test_sample_ising(self):
        h = {"a": 0, "b": -1}
        J = {("a", "b"): -1}

        resp = RotorModelAnnealingSampler().sample_ising(h, J)

        row, col = resp.record.sample.shape

        self.assertEqual(col, 2)  # should get back two variables
        self.assertIs(resp.vartype, dimod.SPIN)  # should be ising

    def test_sample_qubo(self):
        Q = {(0, 1): 1}
        resp = RotorModelAnnealingSampler().sample_qubo(Q)

        row, col = resp.record.sample.shape

        self.assertEqual(col, 2)  # should get back two variables
        self.assertIs(resp.vartype, dimod.BINARY)  # should be qubo

    def test_basic_response(self):
        sampler = RotorModelAnnealingSampler()
        h = {"a": 0, "b": -1}
        J = {("a", "b"): -1}
        response = sampler.sample_ising(h, J)

        self.assertIsInstance(
            response, dimod.SampleSet, "Sampler returned an unexpected response type"
        )

    def test_num_reads(self):
        sampler = RotorModelAnnealingSampler()

        h = {}
        J = {("a", "b"): 0.5, (0, "a"): -1, (1, "b"): 0.0}

        for num_reads in (1, 10, 100, 3223, 10352):
            response = sampler.sample_ising(h, J, num_reads=num_reads)
            row, col = response.record.sample.shape

            self.assertEqual(row, num_reads)
            self.assertEqual(col, 4)

        for bad_num_reads in (0, -1, -100):
            with self.assertRaises(ValueError):
                sampler.sample_ising(h, J, num_reads=bad_num_reads)

        for bad_num_reads in (3.5, float("inf"), "string", [], {}):
            with self.assertRaises(TypeError):
                sampler.sample_ising(h, J, num_reads=bad_num_reads)

    def test_empty_problem(self):
        sampler = RotorModelAnnealingSampler()
        h = {"a": 0, "b": -1}
        J = {("a", "b"): -1}
        eh, eJ = {}, {}
        beta_range = [0.1, 1]  # To suppress warning.
        for h in (h, eh):
            for J in (J, eJ):
                _h = copy.deepcopy(h)
                _J = copy.deepcopy(J)
                sampler.sample_ising(_h, _J, beta_range=beta_range)
        # An empty problem does not allow for beta_range automation
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sampler.sample_ising(eh, eJ)

    def test_seed(self):
        sampler = RotorModelAnnealingSampler()
        num_vars = 40
        h = {v: -1 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 1000

        # test seed exceptions
        for bad_seed in (3.5, float("inf"), "string", [], {}):
            self.assertRaises(TypeError, sampler.sample_ising, {}, {}, seed=bad_seed)
        for bad_seed in (-1, -100, 2**65):
            self.assertRaises(ValueError, sampler.sample_ising, {}, {}, seed=bad_seed)

        # no need to do a bunch of sweeps, in fact the less we do the more
        # sure we can be that the same seed is returning the same result
        all_samples = []

        for seed in (1, 25, 2352, 736145, 5682453):
            response0 = sampler.sample_ising(
                h, J, num_reads=num_reads, num_sweeps=10, seed=seed
            )
            response1 = sampler.sample_ising(
                h, J, num_reads=num_reads, num_sweeps=10, seed=seed
            )

            samples0 = response0.record.sample
            samples1 = response1.record.sample

            self.assertTrue(
                np.array_equal(samples0, samples1),
                "Same seed returned different results",
            )

            for previous_sample in all_samples:
                self.assertFalse(
                    np.array_equal(samples0, previous_sample),
                    "Different seed returned same results",
                )

            all_samples.append(samples0)

    def test_disconnected_problem(self):
        sampler = RotorModelAnnealingSampler()
        h = {}
        J = {
            # K_3
            (0, 1): -1,
            (1, 2): -1,
            (0, 2): -1,
            # disonnected K_3
            (3, 4): -1,
            (4, 5): -1,
            (3, 5): -1,
        }

        resp = sampler.sample_ising(h, J, num_sweeps=1000, num_reads=100)

        row, col = resp.record.sample.shape

        self.assertEqual(row, 100)
        self.assertEqual(col, 6)  # should get back two variables
        self.assertIs(resp.vartype, dimod.SPIN)  # should be ising

    def test_interrupt_error(self):
        sampler = RotorModelAnnealingSampler()
        num_vars = 40
        h = {v: -1 for v in range(num_vars)}
        J = {(u, v): -1 for u in range(num_vars) for v in range(u, num_vars) if u != v}
        num_reads = 100

        def f():
            raise NotImplementedError

        resp = sampler.sample_ising(h, J, num_reads=num_reads, interrupt_function=f)

        self.assertEqual(len(resp), 1)

    def test_sampleset_initial_states(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {"ab": 1, "bc": 1, "ca": 1})
        initial_states = dimod.SampleSet.from_samples_bqm(
            {"a": 1, "b": -1, "c": 1}, bqm
        )

        response = RotorModelAnnealingSampler().sample(
            bqm,
            initial_states=initial_states,
            num_reads=1,
            proposal_acceptance_criteria="MetropolisNonErgodic",
        )

        self.assertEqual(len(response), 1)
        self.assertEqual(response.first.energy, -1)

    def test_initial_states_generator(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {"ab": -1, "bc": 1, "ac": 1})
        init = dimod.SampleSet.from_samples_bqm(
            [{"a": 1, "b": 1, "c": 1}, {"a": -1, "b": -1, "c": -1}], bqm
        )
        sampler = RotorModelAnnealingSampler()

        # 2 fixed initial state, 8 random
        resp = sampler.sample(bqm, initial_states=init, num_reads=10)
        self.assertEqual(len(resp), 10)

        # 2 fixed initial states, 8 random, explicit
        resp = sampler.sample(
            bqm, initial_states=init, initial_states_generator="random", num_reads=10
        )
        self.assertEqual(len(resp), 10)

        # all random
        resp = sampler.sample(bqm, initial_states_generator="random", num_reads=10)
        self.assertEqual(len(resp), 10)

        # all random
        resp = sampler.sample(bqm, num_reads=10)
        self.assertEqual(len(resp), 10)

        # zero-length init states in tuple format, extended by random samples
        zero_init_tuple = (np.empty((0, 3)), ["a", "b", "c"])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            resp = sampler.sample(bqm, initial_states=zero_init_tuple, num_reads=10)
        self.assertEqual(len(resp), 10)

        # explicit None for initial_states should use one random init state
        resp = sampler.sample(bqm, initial_states=None)
        self.assertEqual(len(resp), 1)

        # initial_states truncated to num_reads?
        resp = sampler.sample(
            bqm, initial_states=init, initial_states_generator="none", num_reads=1
        )
        self.assertEqual(len(resp), 1)

        resp = sampler.sample(
            bqm, initial_states=init, initial_states_generator="tile", num_reads=1
        )
        self.assertEqual(len(resp), 1)

        resp = sampler.sample(
            bqm, initial_states=init, initial_states_generator="random", num_reads=1
        )
        self.assertEqual(len(resp), 1)

        # 2 fixed initial states, repeated 5 times
        resp = sampler.sample(
            bqm, initial_states=init, initial_states_generator="tile", num_reads=10
        )
        self.assertEqual(len(resp), 10)

        # can't tile empty states
        with self.assertRaises(ValueError):
            sampler.sample(bqm, initial_states_generator="tile", num_reads=10)

        # not enough initial states
        with self.assertRaises(ValueError):
            sampler.sample(bqm, initial_states_generator="none", num_reads=3)

        # initial_states incompatible with the bqm
        init = dimod.SampleSet.from_samples({"a": 1, "b": 1}, vartype="SPIN", energy=0)
        with self.assertRaises(ValueError):
            sampler.sample(bqm, initial_states=init)

    def test_soft_num_reads(self):
        """Number of reads adapts to initial_states size, if provided."""

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {"ab": -1, "bc": 1, "ac": 1})
        init = dimod.SampleSet.from_samples_bqm(
            [{"a": 1, "b": 1, "c": 1}, {"a": -1, "b": -1, "c": -1}], bqm
        )
        sampler = RotorModelAnnealingSampler()

        # default num_reads == 1
        self.assertEqual(len(sampler.sample(bqm)), 1)
        self.assertEqual(len(sampler.sample(bqm, initial_states_generator="random")), 1)

        # with initial_states, num_reads == len(initial_states)
        self.assertEqual(len(sampler.sample(bqm, initial_states=init)), 2)

        # ... but explicit truncation works too
        self.assertEqual(len(sampler.sample(bqm, initial_states=init, num_reads=1)), 1)

        # if num_reads explicitly given together with initial_states, they are expanded
        self.assertEqual(len(sampler.sample(bqm, initial_states=init, num_reads=3)), 3)

        # if num_reads explicitly given together without initial_states, they are generated
        self.assertEqual(len(sampler.sample(bqm, num_reads=4)), 4)

    def test_0_num_sweeps(self):
        bqm = dimod.BinaryQuadraticModel({}, {"ab": 1}, 0, "SPIN")
        sampleset = dimod.SampleSet.from_samples_bqm(
            [{"a": 1, "b": -1}, {"a": -1, "b": 1}], bqm
        )

        result = RotorModelAnnealingSampler().sample(
            bqm, num_sweeps=0, initial_states=sampleset
        )

        self.assertTrue(np.array_equal(result.record.sample, sampleset.record.sample))
        self.assertEqual(len(result.record.sample), 2)

        result = RotorModelAnnealingSampler().sample(
            bqm,
            num_sweeps=0,
            num_reads=4,
            initial_states=sampleset,
            initial_states_generator="tile",
        )

        expected = np.tile(sampleset.record.sample, (2, 1))

        self.assertTrue(np.array_equal(result.record.sample, expected))
        self.assertEqual(len(result), 4)


class TestHeuristicResponse(unittest.TestCase):
    seed = 1981  # Avoid rare (confusing) probabilistic failures

    def test_job_shop_scheduling_with_linear(self):
        # This test is inherited from SimulatedAnnealingSampler() tests.
        # It is difficult to understand and might be refactored.
        #
        # Set up a job shop scheduling BQM
        #
        # Provide hardcode version of the bqm of "jobs"
        #   jobs = {'b': [(1, 1), (3, 1)],
        #           'o': [(2, 2), (4, 1)],
        #           'g': [(1, 2)]}
        #
        #   There are three jobs: 'b', 'o', 'g'
        #   Each tuple represents a task that runs on a particular machine for
        #   a given amount of time. I.e. (machine_id, duration_on_machine)
        #
        #   Variables below are labelled as '<job_name>_<task_index>,
        #   <task_start_time>'.
        linear = {
            "b_0,0": -2.0,
            "b_0,1": -2.0,
            "b_0,2": -2.0,
            "b_0,3": -2.0,
            "b_1,0": 0.125,
            "b_1,1": -1.5,
            "b_1,2": 0.0,
            "g_0,0": -1.875,
            "g_0,1": -1.5,
            "g_0,2": 0.0,
            "o_0,0": -2.0,
            "o_0,1": -2.0,
            "o_0,2": -2.0,
            "o_1,0": 0.03125,
            "o_1,1": -1.875,
            "o_1,2": -1.5,
            "o_1,3": 0.0,
        }

        quadratic = {
            ("b_0,0", "g_0,0"): 4,
            ("b_0,1", "b_0,0"): 4.0,
            ("b_0,1", "g_0,0"): 2,
            ("b_0,2", "b_0,0"): 4.0,
            ("b_0,2", "b_0,1"): 4.0,
            ("b_0,2", "b_1,2"): 2,
            ("b_0,2", "g_0,1"): 2,
            ("b_0,2", "g_0,2"): 4,
            ("b_0,3", "b_0,0"): 4.0,
            ("b_0,3", "b_0,1"): 4.0,
            ("b_0,3", "b_0,2"): 4.0,
            ("b_0,3", "b_1,2"): 2,
            ("b_0,3", "g_0,2"): 2,
            ("b_1,1", "b_0,1"): 2,
            ("b_1,1", "b_0,2"): 2,
            ("b_1,1", "b_0,3"): 2,
            ("b_1,1", "b_1,2"): 4.0,
            ("g_0,1", "b_0,1"): 4,
            ("g_0,1", "g_0,0"): 4.0,
            ("g_0,2", "g_0,0"): 4.0,
            ("g_0,2", "g_0,1"): 4.0,
            ("o_0,0", "o_1,1"): 2,
            ("o_0,1", "o_0,0"): 4.0,
            ("o_0,1", "o_1,1"): 2,
            ("o_0,1", "o_1,2"): 2,
            ("o_0,2", "o_0,0"): 4.0,
            ("o_0,2", "o_0,1"): 4.0,
            ("o_0,2", "o_1,1"): 2,
            ("o_1,2", "o_0,2"): 2,
            ("o_1,2", "o_1,1"): 4.0,
            ("o_1,3", "o_0,2"): 2,
            ("o_1,3", "o_1,1"): 4.0,
            ("o_1,3", "o_1,2"): 4.0,
        }

        jss_bqm = dimod.BinaryQuadraticModel(linear, quadratic, 9.0, dimod.BINARY)

        # Optimal energy
        optimal_solution = {
            "b_0,0": 1,
            "b_0,1": 0,
            "b_0,2": 0,
            "b_0,3": 0,
            "b_1,0": 0,
            "b_1,1": 1,
            "b_1,2": 0,
            "g_0,0": 0,
            "g_0,1": 1,
            "g_0,2": 0,
            "o_0,0": 1,
            "o_0,1": 0,
            "o_0,2": 0,
            "o_1,0": 0,
            "o_1,1": 0,
            "o_1,2": 1,
            "o_1,3": 0,
        }

        optimal_energy = jss_bqm.energy(optimal_solution)  # Evaluates to 0.5

        # Get heuristic solution
        sampler = RotorModelAnnealingSampler()
        response = sampler.sample(
            jss_bqm, beta_schedule_type="linear", num_reads=10, seed=self.seed
        )
        _, response_energy, _ = next(response.data())

        # Compare energies
        threshold = 0.1  # Arbitrary threshold - revisit when refactored?
        self.assertLess(response_energy, optimal_energy + threshold)

    def test_cubic_lattice_with_geometric(self):
        # This test is inherited from SimulatedAnnealingSampler() tests.
        # It is difficult to understand and might be refactored.
        #
        # Set up all lattice edges in a cube. Each edge is labelled by a 3-D
        # coordinate system
        def get_cubic_lattice_edges(N):
            for x, y, z in itertools.product(range(N), repeat=3):
                u = x, y, z
                yield u, ((x + 1) % N, y, z)
                yield u, (x, (y + 1) % N, z)
                yield u, (x, y, (z + 1) % N)

        # Add a J-bias to each edge
        np_rand = np.random.RandomState(128)
        J = {e: np_rand.choice((-1, 1)) for e in get_cubic_lattice_edges(12)}

        # Solve ising problem
        sampler = RotorModelAnnealingSampler()
        response = sampler.sample_ising(
            {},
            J,
            beta_schedule_type="geometric",
            num_reads=10,
            seed=self.seed,
            proposal_acceptance_criteria="MetropolisNonErgodic",
        )
        _, response_energy, _ = next(response.data())

        # lowest energy found = -3088 with a different benchmarking tool
        threshold = -3000  # Choice is not well-motivated (robust), revisit/refactor
        self.assertLess(
            response_energy,
            threshold,
            f"response_energy, {response_energy}, exceeds threshold, {threshold}",
        )


class TestAngularInformation(unittest.TestCase):
    def test_numpy_input_output(self):
        num_var = 2
        sampler = RotorModelAnnealingSampler()
        num_reads = 2
        h = {i: 1 for i in range(num_var)}
        J = {(i, i + 1): 1 for i in range(num_var - 1)}
        input_val = 33
        valid_states = {True: {0, 128}, False: {input_val, 128 + input_val}}
        costheta = {
            key: {np.cos(val / 128 * np.pi) for val in valid_states[key]}
            for key in valid_states
        }
        # Classical part of energy allowed values (3 decimal places):
        valid_energies1000 = {
            key: {
                round(1000 * (h[0] * ct0 + h[1] * ct1 + J[(0, 1)] * ct0 * ct1))
                for ct0 in costheta[key]
                for ct1 in costheta[key]
            }
            for key in costheta
        }

        for project_states in itertools.product([True, False], [True, False]):
            for num_sweeps in [0, 2]:
                initial_states = np.array(
                    [[input_val] * num_var] * num_reads, dtype=np.uint8
                )  # Variables, states
                ss = sampler.sample_ising(
                    h,
                    J,
                    project_states=project_states,
                    initial_states=initial_states,
                    num_sweeps=num_sweeps,
                    proposal_acceptance_criteria="MetropolisNonErgodic",
                )
                # rotor_state field created if final state not projected.
                self.assertTrue(("rotor_states" not in ss.info) == project_states[1])
                if not project_states[1]:
                    # If projected input, oscillates on [0, 128], else
                    # oscillates on [input_val, 128+input_val]
                    self.assertTrue(
                        set(np.unique(ss.info["rotor_states"])).issubset(
                            valid_states[project_states[0]]
                        )
                    )
                    if num_sweeps == 0 and not project_states[0]:
                        self.assertTrue(
                            np.all(ss.info["rotor_states"] == initial_states)
                        )
                # Energy is calculated prior to projection! Check 3 s.f.
                self.assertTrue(
                    set(np.round(1000 * ss.record.energy)).issubset(
                        valid_energies1000[project_states[0]]
                    )
                )


class TestCoreSpinUpdate(unittest.TestCase):
    sampler = RotorModelAnnealingSampler()
    # Avoid rare (but confusing) false alarms with lose significance threshold
    # and seed.
    seed = 2023
    prng = np.random.default_rng(seed)
    k = 3  # default significance threshold for distributional tests

    def normal_confidence_interval(self, mean_val, std_dev, num_samples, k=k):
        lower_bound = mean_val - k * std_dev / np.sqrt(num_samples)
        upper_bound = mean_val + k * std_dev / np.sqrt(num_samples)
        return lower_bound, upper_bound

    def binomial_confidence_interval(self, p, num_samples, k=k):
        # Spins flip with probability p
        # mean number of flips per sweep mu = num_var*p
        # variance, var = num_var*p*(1-p)
        # A k*sigma interval for number of flips is roughly mu +/-k*root(var)

        mu = num_samples * p
        sig = np.sqrt(num_samples * p * (1 - p))
        upper_bound = mu + k * sig
        lower_bound = mu - k * sig

        return lower_bound, upper_bound

    def test_Metropolis_ergodicity_breaking(self):
        # Default operation, Metropolis sequential order - deterministic in
        # Null BQM (flat energy landscape) given fixed initial condition.
        num_vars = 100
        # test result is independent of the realization, so no need for seed
        init_vector = 1 - 2 * np.random.randint(2, size=num_vars)
        bqm = dimod.BinaryQuadraticModel.from_ising({i: 0 for i in range(num_vars)}, {})
        initial_states = dimod.SampleSet.from_samples_bqm(
            {i: init_vector[i] for i in range(num_vars)}, bqm
        )
        Hp_field = [1]

        # Spins oscillate (ergodicity breaking):
        for num_sweeps in range(1, 3):
            response = RotorModelAnnealingSampler().sample(
                bqm,
                initial_states=initial_states,
                num_reads=1,
                num_sweeps_per_beta=num_sweeps,
                beta_schedule_type="custom",
                Hp_field=Hp_field,
                proposal_acceptance_criteria="MetropolisNonErgodic",
            )
            self.assertTrue(
                np.all(response.record.sample == (-1) ** num_sweeps * init_vector)
            )

    def test_central_limits_random_updates(self):
        # Gibbs update, and random ordering, produce probabilistic results
        # We can however be quite confident in central limits. To avoid rare
        # failures, the seed is set. If a failure is encountered (because the
        # pseudo random number generator is changed and we are simply unlucky
        # a new seed can be hard coded (and that should with high probability
        # resolve the problem in the absence of real bugs.
        num_vars = 10000
        bqm = dimod.BinaryQuadraticModel.from_ising({i: 0 for i in range(num_vars)}, {})
        initial_states = dimod.SampleSet.from_samples_bqm(
            {i: 1 for i in range(num_vars)}, bqm
        )
        k = 3  # Significance threshold
        beta_range = [0.1, 1]  # Bypass ill-conditioned case
        # Gibbs sequential order sweep (test of central limit):
        p = 0.5
        lower_bound, upper_bound = self.binomial_confidence_interval(p, num_vars)
        response = RotorModelAnnealingSampler().sample(
            bqm,
            initial_states=initial_states,
            num_reads=1,
            seed=self.seed,
            num_sweeps=1,
            proposal_acceptance_criteria="GibbsNonErgodic",
            beta_range=beta_range,
        )
        stat = np.sum(response.record.sample == 1)
        self.assertLess(stat, upper_bound)
        self.assertGreater(stat, lower_bound)
        # Metropolis random order sweep (test of central limit):
        # A spin will flip on a sweep if selected (with replacement) an odd
        # number of times. We anticipate roughly Poissonian statistics for
        # large num_var (large enough here). P(x = num selections) =
        # exp(-1)/x!, P(selected twice) = (1/2)^2 exp(-1)/2!, etc.
        # p = P(flipped) = exp(-1)*[1 + 1/2! + 1/4! ..] = exp(-1)*cosh(1)
        #   = 0.568
        # Roughly a Bernouilli random number, hence:
        p = np.cosh(1) * np.exp(-1)
        lower_bound, upper_bound = self.binomial_confidence_interval(p, num_vars)
        response = RotorModelAnnealingSampler().sample(
            bqm,
            initial_states=initial_states,
            num_reads=1,
            seed=self.seed,
            num_sweeps=1,
            randomize_order=True,
            beta_range=beta_range,
            proposal_acceptance_criteria="MetropolisNonErgodic",
        )
        stat = np.sum(response.record.sample == 1)
        self.assertLess(stat, upper_bound)
        self.assertGreater(stat, lower_bound)
        # Gibbs randomized order: 50: 50 on states sampled once.
        p = 1 / np.exp(1) + (1 - 1 / np.exp(1)) * 0.5
        lower_bound, upper_bound = self.binomial_confidence_interval(p, num_vars)
        response = RotorModelAnnealingSampler().sample(
            bqm,
            initial_states=initial_states,
            num_reads=1,
            seed=self.seed,
            num_sweeps=1,
            randomize_order=True,
            proposal_acceptance_criteria="GibbsNonErgodic",
            beta_range=beta_range,
        )
        stat = np.sum(response.record.sample == 1)
        self.assertLess(stat, upper_bound)
        self.assertGreater(stat, lower_bound)
        # Gibbs with energy signal, regardless of num_sweeps and initial
        # condition expect +1 state with probability
        # exp(beta_final)/[exp(beta_final) + exp(-beta_final)] on every
        # updated state (all states given sequential order)
        bqm = dimod.BinaryQuadraticModel.from_ising({i: 1 for i in range(num_vars)}, {})
        betas = [1, 1.5]
        p = np.exp(-betas[-1]) / (2 * np.cosh(betas[-1]))
        lower_bound, upper_bound = self.binomial_confidence_interval(p, num_vars, k)
        response = RotorModelAnnealingSampler().sample(
            bqm,
            initial_states=initial_states,
            num_reads=1,
            seed=self.seed,
            proposal_acceptance_criteria="GibbsNonErgodic",
            beta_schedule_type="custom",
            Hp_field=betas,
            num_sweeps_per_beta=1,
        )
        stat = np.sum(response.record.sample == 1)
        self.assertLess(stat, upper_bound)
        self.assertGreater(stat, lower_bound)

    def test_greedy_limit_independent_spins(self):
        num_vars = 10000
        # test result is independent of the realization, so no need for seed
        init_vector = np.random.normal(size=num_vars)
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {i: init_vector[i] for i in range(num_vars)}, {}
        )
        beta_schedule_type = "custom"
        beta_schedule = [float("inf")]
        # Check escape from/to trivial ground state, all methods:
        ground_state_vec = np.array(
            [-int(np.sign(bqm.linear[i])) for i in range(num_vars)]
        )
        ground_state = dimod.SampleSet.from_samples_bqm(
            {i: ground_state_vec[i] for i in range(num_vars)}, bqm
        )
        sky_state = dimod.SampleSet.from_samples_bqm(
            {i: -ground_state_vec[i] for i in range(num_vars)}, bqm
        )
        # All touched spins escape to the ground state.
        p = 1 - np.exp(-1)  # ~Probability index sampled atleast once:
        lower_bound, upper_bound = self.binomial_confidence_interval(p, num_vars)
        for randomize_order in [False, True]:
            for proposal_acceptance_criteria in [
                "MetropolisNonErgodic",
                "GibbsNonErgodic",
            ]:
                response = RotorModelAnnealingSampler().sample(
                    bqm,
                    initial_states=ground_state,
                    num_reads=1,
                    num_sweeps=1,
                    randomize_order=randomize_order,
                    proposal_acceptance_criteria=proposal_acceptance_criteria,
                    beta_schedule_type=beta_schedule_type,
                    Hp_field=beta_schedule,
                    seed=self.seed,
                )
                self.assertTrue(np.all(response.record.sample == ground_state_vec))
                response = RotorModelAnnealingSampler().sample(
                    bqm,
                    initial_states=sky_state,
                    num_reads=1,
                    num_sweeps=1,
                    randomize_order=randomize_order,
                    proposal_acceptance_criteria=proposal_acceptance_criteria,
                    beta_schedule_type=beta_schedule_type,
                    Hp_field=beta_schedule,
                    seed=self.seed,
                )
                if randomize_order:
                    # Partial recovery only (inline with sampled indices)
                    stat = np.sum(response.record.sample == ground_state_vec)
                    self.assertLess(stat, upper_bound)
                    self.assertGreater(stat, lower_bound)
                else:
                    # Recovers ground state
                    self.assertTrue(np.all(response.record.sample == ground_state_vec))

    def test_independent_angles(self):
        # Beginning from uniform angles
        # all methods should return uniform angles w.h.p
        num_reads = 4
        num_vars = 64
        disc = 256
        bqm = dimod.BinaryQuadraticModel.from_ising({i: 0 for i in range(num_vars)}, {})
        uid_states = self.prng.integers(
            256, size=(num_reads, num_vars), dtype=np.uint8
        )  # 256 disc. angles
        # Simple test for uniform with few samples.
        # p angle unrepresented = 1-exp(-1)
        lower_bound, upper_bound = self.binomial_confidence_interval(
            1 - np.exp(-(num_reads * num_vars) / disc), num_reads * num_vars
        )
        for randomize_order in [False, True]:
            for proposal_acceptance_criteria in [
                "MetropolisNonErgodic",
                "GibbsNonErgodic",
                "MetropolisUniform",
                "MetropolisTF",
            ]:
                response = RotorModelAnnealingSampler().sample(
                    bqm,
                    initial_states=uid_states.copy(),
                    randomize_order=randomize_order,
                    proposal_acceptance_criteria=proposal_acceptance_criteria,
                    beta_range=[0.1, 1],
                    project_states=(False, False),
                    seed=self.seed,
                )
                # angles with representation.
                stat = len(np.unique(response.info["rotor_states"]))
                self.assertLess(stat, upper_bound)
                self.assertGreater(stat, lower_bound)

        proposal_acceptance_criteria = "MetropolisUniform"
        randomize_order = False
        # Check that Delta function is immediately dispersed to u.i.d
        delta_states = 64 * np.ones(shape=(num_reads, num_vars), dtype=np.uint8)
        response = RotorModelAnnealingSampler().sample(
            bqm,
            initial_states=delta_states,
            randomize_order=randomize_order,
            proposal_acceptance_criteria=proposal_acceptance_criteria,
            beta_range=[0.1, 1],
            project_states=(False, False),
            seed=self.seed,
            num_sweeps=1,
        )
        # angles with representation.
        stat = len(np.unique(response.info["rotor_states"]))
        self.assertLess(stat, upper_bound)
        self.assertGreater(stat, lower_bound)

        proposal_acceptance_criteria = "MetropolisTF"
        randomize_order = False
        # Check that TF takes angles locally Hd/Hp=32 -> [32-delta,32+delta]
        # TO DO: Check 1 sided for large Hamiltonian (greedy)
        trans_fields = [0] * num_vars
        Hp_field = [1]
        num_reads = 1000
        root = 32
        delta_states = root * np.ones(shape=(num_reads, num_vars), dtype=np.uint8)
        for Gamma in [0, 1, 4, 16, 123, 255, 256]:
            Hd_field = [(Gamma + self.prng.random()) / 256]
            response = RotorModelAnnealingSampler().sample(
                bqm,
                initial_states=delta_states.copy(),
                randomize_order=randomize_order,
                proposal_acceptance_criteria=proposal_acceptance_criteria,
                beta_schedule_type="custom",
                Hp_field=Hp_field,
                Hd_field=Hd_field,
                project_states=(False, False),
                trans_fields=trans_fields,
                num_sweeps_per_beta=1,
                num_sweeps=1,
            )
            # Angles reachable
            us = np.unique(response.info["rotor_states"])
            if Hd_field[-1] == 0:
                osi = 0
            elif Hd_field[-1] < 1 / 256:
                osi = 1
            else:
                osi = int(Hd_field[-1] * 256 + 1) // 2
            # In range:
            self.assertTrue(
                set(us).issubset({(root + i) % 256 for i in range(-osi, osi + 1)})
            )
            # Expected absolute deviation should match (concentrate in the mean) at
            # (Hd/Hp)*(pi/2), pi/2 represented by 64
            mean_val = Hd_field[-1] * 64
            std_dev = Hd_field[-1] * 64  # Approx
            lower_bound, upper_bound = self.normal_confidence_interval(
                mean_val, std_dev, num_samples=num_reads * num_vars
            )

            stat = np.mean(
                np.abs(np.array(response.info["rotor_states"], dtype=float) - root)
            )
            self.assertLess(mean_val, upper_bound)
            self.assertGreater(mean_val, lower_bound)
            # print(np.histogram(response.info['rotor_states'],
            #     list(us) + [us[-1]+1]))
            # Could also check uniform over central values and symmetric.
            # For now, I've seen it with my eyes!

    def test_equilibrated_statistics_classical_limit(self):
        # Under sequential Gibbs sampling of independent spins 1 sweep
        # is sufficient to equilibrate at P(1) = (1+tanh(beta*Hd*h_i))/2.
        # We can test this at population level.

        # This is a multi-test, but can still be interpretted as k~3
        # sigma confidence interval with only 4 variables.
        num_vars = 4
        num_reads = 1000
        prng = np.random.default_rng(self.seed)
        init_vector = prng.random(size=num_vars)  # reproducible for clarity.
        Hp_field = [1, 2]
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {i: -init_vector[i] for i in range(num_vars)}, {}
        )
        lower_bound = np.zeros(num_vars)
        upper_bound = np.zeros(num_vars)
        for i in bqm.variables:
            p = np.exp(-Hp_field[-1] * bqm.linear[i]) / (
                2 * np.cosh(Hp_field[-1] * bqm.linear[i])
            )
            lower_bound[i], upper_bound[i] = self.binomial_confidence_interval(
                p, num_reads
            )

        response = RotorModelAnnealingSampler().sample(
            bqm,
            num_reads=num_reads,
            seed=self.seed,
            proposal_acceptance_criteria="GibbsNonErgodic",
            randomize_order=False,
            beta_schedule_type="custom",
            Hp_field=Hp_field,
            num_sweeps_per_beta=1,
        )
        stat = np.sum(response.record.sample == 1, axis=0)
        for i in range(num_vars):
            self.assertLess(stat[i], upper_bound[i])
            self.assertGreater(stat[i], lower_bound[i])

    def test_equilibrated_statistics_quantum_limit(self):
        # P(angle) ~ exp( Hd_array sum_i gamma_i sin(angle_i))

        num_vars = 4
        num_reads = 1000
        # prng = np.random.default_rng(self.seed)
        # trans_fields = prng.random(size=num_vars)  # reproducible
        # Add trans_fields later!
        Gamma = 1
        num_sweeps = 1000
        # Requires some time to build up correct number of breaks.
        # Can ignore possibility of non-equilibration after ~40 sweeps
        # (heuristic) at small Hd_field.
        Hp_field = [2] * num_sweeps
        Hd_field = [Gamma] * num_sweeps
        bqm = dimod.BinaryQuadraticModel.from_ising({i: 0 for i in range(num_vars)}, {})

        # Poisson-like statistics. mean ~ var
        states = np.arange(256) / 256 * (2 * np.pi)
        sin_states = np.abs(np.sin(states))
        Pstate = np.exp(Hd_field[-1] * (sin_states - 1))
        Pstate = Pstate / np.sum(Pstate)
        EX = np.sum(Pstate * sin_states)
        sigX = np.sqrt(np.sum(Pstate * sin_states**2) - EX**2)
        lower_bound, upper_bound = self.normal_confidence_interval(
            EX, sigX, num_reads * num_vars
        )

        response = RotorModelAnnealingSampler().sample(
            bqm,
            num_reads=num_reads,
            seed=self.seed,
            beta_schedule_type="custom",
            Hp_field=Hp_field,
            randomize_order=False,
            proposal_acceptance_criteria="MetropolisUniform",
            Hd_field=Hd_field,
            num_sweeps_per_beta=1,
            project_states=(True, False),
        )
        stat = np.mean(np.abs(np.sin(2 * np.pi * response.info["rotor_states"] / 256)))
        self.assertLess(stat, upper_bound)
        self.assertGreater(stat, lower_bound)

    def test_UniformVsTF_angles(self):
        # Check TF and Uniform are exactly matched in special case:
        num_reads = 10
        bqm = dimod.BinaryQuadraticModel.from_ising({0: 0.1}, {(0, 1): 0.1})
        samples = {}
        for pac in ["MetropolisUniform", "MetropolisTF"]:
            response = RotorModelAnnealingSampler().sample(
                bqm,
                num_reads=num_reads,
                seed=self.seed,
                beta_schedule_type="custom",
                Hp_field=[1, 1],
                proposal_acceptance_criteria="MetropolisUniform",
                Hd_field=[1.5, 1],
                num_sweeps_per_beta=1,
                project_states=(True, False),
            )
            samples[pac] = response.record.sample
        self.assertTrue(np.all(samples["MetropolisUniform"] == samples["MetropolisTF"]))

    def test_equilibrated_statistics_pair(self):
        # Non trivial Hamiltonian, check reproduces equilibrium:
        # 1*(0.8*sigZ_1 sigZ_2 - 0.1*sigZ_1 + 0.2*sigZ_2) - 1.1*(1.1 sigX_1 - 0.5*sigX_2)
        J = 0.8
        h = [-0.1, 0.2]
        trans_fields = [1.1, 0.5]

        theta = np.arange(256) / 128 * np.pi
        cos_theta = np.cos(theta)  # mag
        sin_theta = np.abs(np.sin(theta))  # alignment
        Hd_field = [1.1] * 100  # 100 sweeps, sufficient for equilibration
        Hp_field = [1] * 100
        EnQ = (
            -trans_fields[0] * sin_theta[:, np.newaxis]
            - trans_fields[1] * sin_theta[np.newaxis, :]
        )
        EnC = (
            h[0] * cos_theta[:, np.newaxis]
            + h[1] * cos_theta[np.newaxis, :]
            + J * cos_theta[:, np.newaxis] * cos_theta[np.newaxis, :]
        )
        measure = np.exp(
            -Hd_field[-1] * (EnQ - np.max(EnQ)) - Hp_field[-1] * (EnC - np.max(EnC))
        )
        Z = np.sum(measure)
        # Fairly sample???? Easier to just trust mixing of the chain for this example (no large barriers)
        # cdf = np.cumsum(measure)/Z
        # Under projection, quantities are linear:
        stats = {
            "mag1": cos_theta[:, np.newaxis],
            "mag2": cos_theta[np.newaxis, :],
            "corr": cos_theta[:, np.newaxis] * cos_theta[np.newaxis, :],
            "align1": sin_theta[:, np.newaxis],
            "align2": sin_theta[np.newaxis, :],
        }
        mu = {}
        sig = {}
        for stat in stats:
            mu[stat] = np.sum(measure * stats[stat]) / Z
            sig[stat] = np.sqrt(
                np.sum(measure * (stats[stat] * stats[stat])) / Z - mu[stat] * mu[stat]
            )

        num_reads = 1000
        for proposal_acceptance_criteria in ["MetropolisUniform", "MetropolisTF"]:
            response = RotorModelAnnealingSampler().sample_ising(
                {idx: val for idx, val in enumerate(h)},
                {(0, 1): J},
                num_reads=num_reads,
                seed=self.seed,
                beta_schedule_type="custom",
                Hp_field=Hp_field,
                trans_fields=trans_fields,
                Hd_field=Hd_field,
                randomize_order=False,
                project_states=(False, False),
                proposal_acceptance_criteria=proposal_acceptance_criteria,
                num_sweeps_per_beta=1,
            )
            # This is a multitest, but using 3 sigma threshold likely to pass
            equil_states = response.info["rotor_states"]

            for key, stat in stats.items():
                lower_bound, upper_bound = self.normal_confidence_interval(
                    mu[key], sig[key], num_samples=num_reads
                )
                shape = stat.shape
                if shape[0] == 1:
                    mean_val = np.mean(
                        [stat[0, rotor_state] for rotor_state in equil_states[:, 1]]
                    )
                elif shape[1] == 1:
                    mean_val = np.mean(
                        [stat[rotor_state, 0] for rotor_state in equil_states[:, 0]]
                    )
                else:
                    mean_val = np.mean(
                        [
                            stat[rotor_state[1], rotor_state[0]]
                            for rotor_state in equil_states
                        ]
                    )
                self.assertLess(mean_val, upper_bound)
                self.assertGreater(mean_val, lower_bound)


class TestSerializable(unittest.TestCase):

    def is_jsonable(self, x):
        try:
            json.dumps(x)
            return True
        except TypeError:
            return False

    def test_serializable(self):
        sampler = RotorModelAnnealingSampler()
        for is_ser in [False, True]:
            resp = sampler.sample_ising(
                {0: 1},
                {},
                beta_range=np.array([0.1, 2]),
                make_info_json_serializable=is_ser,
            )
            self.assertEqual(self.is_jsonable(resp.info), is_ser)


class TestStatistics(unittest.TestCase):
    def test_statistics(self):
        sampler = RotorModelAnnealingSampler()
        num_var = 10
        num_reads = 3
        num_sweeps = 4
        variables = np.arange(num_var)
        for ssi in [1, 2]:
            resp = sampler.sample_ising(
                {i: 1 for i in variables},
                {},
                num_sweeps=num_sweeps,
                num_reads=num_reads,
                schedule_sample_interval=ssi,
            )
            self.assertIn("statistics", resp.info)
            self.assertEqual(
                resp.info["statistics"].shape, (num_reads, num_sweeps // ssi, num_var)
            )
            self.assertTrue(
                set(np.unique(resp.info["statistics"])).issubset(
                    {i for i in range(256)}
                )
            )


if __name__ == "__main__":
    unittest.main()
