# Copyright 2018 D-Wave Systems Inc.
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

"""Test the TabuSampler python interface."""

import unittest

import dimod
import numpy as np

import dwave.samplers.tabu as tabu
from dwave.samplers.tabu.utils import tictoc


@dimod.testing.load_sampler_bqm_tests(tabu.TabuSampler)
class TestTabuSampler(unittest.TestCase):

    def test_instantiation(self):
        sampler = tabu.TabuSampler()
        dimod.testing.assert_sampler_api(sampler)

    def test_sample_basic(self):
        sampler = tabu.TabuSampler()

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {'ab': -1, 'bc': +1, 'ac': +1})

        resp = sampler.sample(bqm)

        dimod.testing.assert_response_energies(resp, bqm)

    def test_sample_num_reads(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {'ab': -1, 'bc': +1, 'ac': +1})

        resp = tabu.TabuSampler().sample(bqm, num_reads=57)
        dimod.testing.assert_response_energies(resp, bqm)
        self.assertEqual(sum(resp.record.num_occurrences), 57)

    def test_disconnected_problem(self):
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

        bqm = dimod.BinaryQuadraticModel.from_ising(h, J)
        resp = tabu.TabuSampler().sample(bqm)
        dimod.testing.assert_response_energies(resp, bqm)

    def test_empty(self):
        resp = tabu.TabuSampler().sample(
            dimod.BinaryQuadraticModel.empty(dimod.SPIN))
        dimod.testing.assert_response_energies(
            resp, dimod.BinaryQuadraticModel.empty(dimod.SPIN))

        resp = tabu.TabuSampler().sample(
            dimod.BinaryQuadraticModel.empty(dimod.BINARY))
        dimod.testing.assert_response_energies(
            resp, dimod.BinaryQuadraticModel.empty(dimod.BINARY))

        resp = tabu.TabuSampler().sample_qubo({})
        dimod.testing.assert_response_energies(
            resp, dimod.BinaryQuadraticModel.empty(dimod.BINARY))

        resp = tabu.TabuSampler().sample_ising({}, {})
        dimod.testing.assert_response_energies(
            resp, dimod.BinaryQuadraticModel.empty(dimod.SPIN))

    def test_single_variable_problem(self):
        bqm = dimod.BinaryQuadraticModel({'a': 1}, {}, 0.0, dimod.SPIN)
        resp = tabu.TabuSampler().sample(bqm)
        dimod.testing.assert_response_energies(resp, bqm)
        self.assertEqual(resp.first.energy, -1)

    def test_linear_problem(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({v: -1 for v in range(100)}, {})
        resp = tabu.TabuSampler().sample(bqm)
        dimod.testing.assert_response_energies(resp, bqm)

    def test_initial_states_smoketest(self):
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {'ab': -1, 'bc': 1, 'ac': 1})
        resp = tabu.TabuSampler().sample(bqm, initial_states=tabu.TabuSampler().sample(bqm))
        dimod.testing.assert_response_energies(resp, bqm)

    def test_initial_states(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {'ab': -1, 'bc': 1, 'ac': 1})
        init = dimod.SampleSet.from_samples(
            {'a': 0, 'b': 0, 'c': 0}, vartype=dimod.BINARY, energy=0)

        resp = tabu.TabuSampler().sample(bqm, initial_states=init)
        dimod.testing.assert_response_energies(resp, bqm)

    def test_initial_states_generator(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {'ab': -1, 'bc': 1, 'ac': 1})
        init = dimod.SampleSet.from_samples_bqm(
            [{'a': 1, 'b': 1, 'c': 1},
             {'a': -1, 'b': -1, 'c': -1}], bqm)

        # 2 fixed initial state, 8 random
        resp = tabu.TabuSampler().sample(
            bqm, initial_states=init, num_reads=10)
        self.assertEqual(len(resp), 10)

        # 2 fixed initial states, 8 random, explicit
        resp = tabu.TabuSampler().sample(
            bqm, initial_states=init, initial_states_generator='random', num_reads=10)
        self.assertEqual(len(resp), 10)

        # all random
        resp = tabu.TabuSampler().sample(
            bqm, initial_states_generator='random', num_reads=10)
        self.assertEqual(len(resp), 10)

        # all random
        resp = tabu.TabuSampler().sample(bqm, num_reads=10)
        self.assertEqual(len(resp), 10)

        # initial_states truncated to num_reads?
        resp = tabu.TabuSampler().sample(
            bqm, initial_states=init, initial_states_generator='none', num_reads=1)
        self.assertEqual(len(resp), 1)

        resp = tabu.TabuSampler().sample(
            bqm, initial_states=init, initial_states_generator='tile', num_reads=1)
        self.assertEqual(len(resp), 1)

        resp = tabu.TabuSampler().sample(
            bqm, initial_states=init, initial_states_generator='random', num_reads=1)
        self.assertEqual(len(resp), 1)

        # 2 fixed initial states, repeated 5 times
        resp = tabu.TabuSampler().sample(
            bqm, initial_states=init, initial_states_generator='tile', num_reads=10)
        self.assertEqual(len(resp), 10)

        # can't tile empty states
        with self.assertRaises(ValueError):
            resp = tabu.TabuSampler().sample(
                bqm, initial_states_generator='tile', num_reads=10)

        # not enough initial states
        with self.assertRaises(ValueError):
            resp = tabu.TabuSampler().sample(
                bqm, initial_states_generator='none', num_reads=3)

        # initial_states incompatible with the bqm
        init = dimod.SampleSet.from_samples(
            {'a': 1, 'b': 1}, vartype='SPIN', energy=0)
        with self.assertRaises(ValueError):
            resp = tabu.TabuSampler().sample(bqm, initial_states=init)

    def test_input_validation(self):
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {}, {'ab': -1, 'bc': 1, 'ac': 1})
        empty = dimod.BinaryQuadraticModel.empty(dimod.SPIN)

        # empty bqm
        self.assertEqual(len(tabu.TabuSampler().sample(empty)), 0)

        # invalid tenure type
        with self.assertRaises(TypeError):
            tabu.TabuSampler().sample(bqm, tenure=2.0)

        # invalid tenure range
        with self.assertRaises(ValueError):
            tabu.TabuSampler().sample(bqm, tenure=100)

        # invalid num_reads type
        with self.assertRaises(TypeError):
            tabu.TabuSampler().sample(bqm, num_reads=10.0)

        # invalid num_reads range
        with self.assertRaises(ValueError):
            tabu.TabuSampler().sample(bqm, num_reads=0)

        # invalid initial_states type
        with self.assertRaises(ValueError):
            tabu.TabuSampler().sample(bqm, initial_states=[])

        with self.assertRaises(ValueError):
            tabu.TabuSampler().sample(bqm, initial_states_generator='non-existing')

        # invalid initial_states length
        with self.assertRaises(ValueError):
            tabu.TabuSampler().sample(bqm, initial_states=[1, 1])

    def test_soft_num_reads(self):
        """Number of reads adapts to initial_states size, if provided."""

        bqm = dimod.BinaryQuadraticModel.from_ising({}, {'ab': -1, 'bc': 1, 'ac': 1})
        init = dimod.SampleSet.from_samples_bqm([{'a': 1, 'b': 1, 'c': 1},
                                                 {'a': -1, 'b': -1, 'c': -1}], bqm)

        # default num_reads == 1
        self.assertEqual(len(tabu.TabuSampler().sample(bqm)), 1)

        # with initial_states, num_reads == len(initial_states)
        self.assertEqual(len(tabu.TabuSampler().sample(bqm, initial_states=init)), 2)

        # if explicitly given, with initial_states, they are expanded
        self.assertEqual(len(tabu.TabuSampler().sample(bqm, initial_states=init, num_reads=3)), 3)

        # if explicitly given, without initial_states, they are generated
        self.assertEqual(len(tabu.TabuSampler().sample(bqm, num_reads=4)), 4)

    def test_seed(self):
        sampler = tabu.TabuSampler()

        # use a bqm with a huge plateau (extreme ground state degeneracy)
        # to give the tabu sampler many chances of flipping (exercising rng)
        bqm = dimod.generators.random.randint(10, dimod.SPIN, low=1, high=1)
        tenure = 5

        all_samples = []
        for seed in (1, 25, 2352):
            response0 = sampler.sample(
                bqm, num_reads=1, tenure=tenure, num_restarts=1, timeout=None, seed=seed)
            response1 = sampler.sample(
                bqm, num_reads=1, tenure=tenure, num_restarts=1, timeout=None, seed=seed)

            samples0 = response0.record.sample
            samples1 = response1.record.sample

            self.assertTrue(np.array_equal(samples0, samples1),
                            "Same seed returned different results")

            for previous_sample in all_samples:
                self.assertFalse(np.array_equal(samples0, previous_sample),
                                 "Different seed returned same results")

            all_samples.append(samples0)

    def test_timeout(self):
        sampler = tabu.TabuSampler()
        bqm = dimod.BinaryQuadraticModel.from_ising({}, {'ab': -1, 'bc': 1, 'ac': 1})

        with tictoc() as tt:
            response = sampler.sample(bqm, num_reads=1, timeout=500, seed=123)
        self.assertAlmostEqual(tt.dt, 0.5, places=1)

        with tictoc() as tt:
            response = sampler.sample(bqm, num_reads=3, timeout=200, seed=123)
        self.assertAlmostEqual(tt.dt, 0.6, places=1)

        # Run as simple-tabu-search with timeout:
        with tictoc() as tt:
            response = sampler.sample(bqm, num_reads=2, timeout=300, seed=123,
                                      num_restarts=0, lower_bound_z=2147483647)
        self.assertAlmostEqual(tt.dt, 0.6, places=1)

    def test_num_restarts(self):
        sampler = tabu.TabuSampler()
        bqm = dimod.generators.random.randint(10, 'SPIN', seed=123)
        target_restarts = 100

        response = sampler.sample(
            bqm, num_reads=1, timeout=None, num_restarts=target_restarts, seed=345)

        num_restarts = response.record['num_restarts']
        self.assertEqual(target_restarts, num_restarts)

    def test_energy_threshold(self):
        sampler = tabu.TabuSampler()
        bqm = dimod.generators.random.randint(100, 'SPIN', seed=123)
        energy_threshold = -400

        # Expect that energy_threshold is met before timeout
        with tictoc() as tt:
            response = sampler.sample(
                bqm, timeout=100000, energy_threshold=energy_threshold, seed=345)

        self.assertLessEqual(tt.dt, 1.0)

    def test_coeff_z(self):
        # bqm with excited state initialization should relax if and only if
        # updated:
        # num_var large enough that 'steepest ascent' subroutine does not
        # cover all variables:
        # single spin model (trivially solved by MIS), with initialization
        # in an excited state
        num_var = 11
        bqm = dimod.BinaryQuadraticModel.from_ising(
            {i: 1 for i in range(num_var)}, {})
        init = dimod.SampleSet.from_samples_bqm(
            [{i: bqm.linear[i] for i in range(num_var)}], bqm)
        sampler = tabu.TabuSampler()

        response = sampler.sample(bqm, num_reads=1, timeout=None,
                                  num_restarts=0,
                                  coefficient_z_first=0, lower_bound_z=0,
                                  initial_states=init)
        # No updates, insufficient to find global minima.
        self.assertEqual(response.record.energy[0], num_var)

        response = sampler.sample(bqm, num_reads=1, timeout=None,
                                  num_restarts=0,
                                  coefficient_z_first=0, lower_bound_z=1,
                                  initial_states=init)
        # Updated once, sufficient to find global minima.
        self.assertEqual(response.record.energy[0], -num_var)

        response = sampler.sample(bqm, num_reads=1, timeout=None,
                                  num_restarts=0,
                                  coefficient_z_first=1, lower_bound_z=0,
                                  initial_states=init)
        # Updated once, sufficient to find global minima.
        self.assertEqual(response.record.energy[0], -num_var)

        response = sampler.sample(bqm, num_reads=1, timeout=None,
                                  num_restarts=1,
                                  coefficient_z_first=0, lower_bound_z=0,
                                  coefficient_z_restart=1,
                                  initial_states=init)
        # Updated once, sufficient to find global minima.
        self.assertEqual(response.record.energy[0], -num_var)

        # subset steepest ascent (meaning descent, if we think of qubo
        # minimization) alone is insufficient to establish global minima, but
        # is sufficient to escape the global maxima:
        response = sampler.sample(bqm, num_reads=1, timeout=None,
                                  num_restarts=1,
                                  coefficient_z_first=0, lower_bound_z=0,
                                  coefficient_z_restart=0,
                                  initial_states=init)
        self.assertLess(response.record.energy[0], num_var)
        self.assertGreater(response.record.energy[0], -num_var)
