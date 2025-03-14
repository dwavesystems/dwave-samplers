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
import numpy as np

from itertools import product

import dimod

from dwave.samplers.sa.sampler import AnnealedImportanceSampling as AIS, log_sum_exp


class TestAnnealedImportanceSampling(unittest.TestCase):

    def test_logsumexp(self):
        input = [-1, 2.0, 3.14]
        self.assertAlmostEqual(log_sum_exp(input), np.log(np.exp(input).sum()))
        self.assertIsNone(log_sum_exp([]))

    def test_z_against_exact(self):
        # NOTE: keep system size `n` and energy scale small because we will make assertions on Z not logZ.
        # Z is unbiased even though logZ is more numerically stable.
        n = 7
        rng = np.random.default_rng(20250225)
        bqm: dimod.BQM = dimod.generators.ran_r(10, n, seed=int(rng.integers(9121)))
        for v in bqm.variables:
            bqm.set_linear(v, rng.standard_normal())
        bqm.normalize(0.01)

        states = list(product([-1, 1], repeat=n))
        sampler = AIS()
        sample_size = 100_000
        for target_beta, randomize, ac in product(
            [0.0, 0.5, 1.0, 3.14], [True, False], ["Metropolis", "Gibbs"]
        ):
            with self.subTest(
                "Single-sweep test for different sample size and target",
                sample_size=sample_size,
                target_beta=target_beta,
                acceptance_criterion=ac,
                randomize_order=randomize,
                num_sweeps=100,
            ):
                energies = bqm.energies(states)
                logz = log_sum_exp(-target_beta * energies)
                sample_set = sampler.sample(
                    bqm,
                    target_beta,
                    num_reads=sample_size,
                    num_sweeps=1,
                    randomize_order=randomize,
                    proposal_acceptance_criteria=ac,
                    seed=int(rng.integers(51823941)),
                )
                logdiff = sample_set.info["logz_estimate"] - logz
                ratio = np.exp(logdiff)
                self.assertAlmostEqual(ratio, 1, 3)

    def test_logz_estimate(self):
        rng = np.random.default_rng(20250224)
        n = 123
        bqm: dimod.BQM = dimod.generators.ran_r(10, n)
        for v, b in zip(bqm.variables, rng.choice([-1.2, 0.0, 3.14], n)):
            if b:
                bqm.set_linear(v, b)
        sampler = AIS()
        for sample_size, target_beta, randomize, ac in product(
            [1, 34], [0.0, 0.5, 1.0, 3.14], [True, False], ["Metropolis", "Gibbs"]
        ):
            with self.subTest(
                "Single-sweep test for different sample size and target",
                sample_size=sample_size,
                target_beta=target_beta,
                acceptance_criterion=ac,
                randomize_order=randomize,
            ):
                sample_set = sampler.sample(
                    bqm,
                    target_beta,
                    num_reads=sample_size,
                    num_sweeps=1,
                    proposal_acceptance_criteria=ac,
                    randomize_order=randomize,
                    seed=int(rng.integers(358858588)),
                )
                logw = sample_set.info["log_weights"]

                self.assertAlmostEqual(
                    log_sum_exp(logw) - np.log(sample_size),
                    sample_set.info["logz_estimate"],
                )

    def test_uniform_distribution(self):
        # NOTE: keep system size `n` and energy scale small because we will make assertions on Z not logZ.
        # Z is unbiased even though logZ is more numerically stable.
        n = 10
        rng = np.random.default_rng(20250225)
        bqm: dimod.BQM = 0 * dimod.generators.ran_r(10, n, seed=int(rng.integers(9121)))
        sampler = AIS()
        sample_size = 10000
        z = 2**n
        for target_beta, randomize, ac in product(
            [0.0, 0.5, 1.0, 3.14], [True, False], ["Metropolis", "Gibbs"]
        ):
            with self.subTest(
                "Single-sweep test for different sample size and target",
                sample_size=sample_size,
                target_beta=target_beta,
                acceptance_criterion=ac,
                randomize_order=randomize,
            ):
                sample_set = sampler.sample(
                    bqm,
                    target_beta,
                    num_reads=sample_size,
                    num_sweeps=1,
                    randomize_order=randomize,
                    proposal_acceptance_criteria=ac,
                    seed=int(rng.integers(51823941)),
                )
                self.assertAlmostEqual(np.exp(sample_set.info["logz_estimate"]), z)

    def test_one_sweep(self):
        # At one sweep, we should have (vanilla) importance sampling.
        rng = np.random.default_rng(20250224)
        n = 123
        bqm: dimod.BQM = dimod.generators.ran_r(10, n)
        for v, b in zip(bqm.variables, rng.choice([-1.2, 0.0, 3.14], n)):
            if b:
                bqm.set_linear(v, b)
        sampler = AIS()
        for sample_size, target_beta, randomize, ac in product(
            [1, 34], [0.0, 0.5, 1.0, 3.14], [True, False], ["Metropolis", "Gibbs"]
        ):
            with self.subTest(
                "Single-sweep test for different sample size and target",
                sample_size=sample_size,
                target_beta=target_beta,
                acceptance_criterion=ac,
                randomize_order=randomize,
            ):
                sample_set = sampler.sample(
                    bqm,
                    target_beta,
                    num_reads=sample_size,
                    num_sweeps=1,
                    proposal_acceptance_criteria=ac,
                    randomize_order=randomize,
                    seed=int(rng.integers(358858588)),
                )
                logw = sample_set.info["log_weights"]
                logw_true = -target_beta * bqm.energies(sample_set) + n * np.log(2)
                self.assertTrue(np.isclose(logw, logw_true).all())

    def test_no_sweep(self):
        self.assertRaises(ValueError, AIS().sample_ising, {0: 3}, {}, num_sweeps=0)
