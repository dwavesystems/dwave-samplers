# Copyright 2022 D-Wave Systems Inc.
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

import datetime
import sys
import typing

import dimod
import numpy as np

from dwave.samplers.random.cyrandom import sample


__all__ = ['RandomSampler']


class RandomSampler(dimod.Sampler):
    """A random sampler, useful as a performance baseline and for testing.

    Examples:

        >>> from dwave.samplers import RandomSampler
        >>> sampler = RandomSampler()

        Create a random binary quadratic model.

        >>> import dimod
        >>> bqm = dimod.generators.gnp_random_bqm(100, .5, 'BINARY')

        Get 20 random samples.

        >>> sampleset = sampler.sample(bqm, num_reads=20)

        Get the best 5 sample found in .1 seconds

        >>> sampleset = sampler.sample(bqm, time_limit=.1, max_num_samples=5)

    """

    parameters: typing.Mapping[str, typing.List] = dict(
        num_reads=[],
        time_limit=[],
        max_num_samples=[],
        seed=[],
        )
    """Keyword arguments accepted by the sampling methods.

    Examples:

        >>> from dwave.samplers import RandomSampler
        >>> sampler = RandomSampler()
        >>> sampler.parameters
        {'num_reads': [], 'time_limit': [], 'max_num_samples': [], 'seed': []}

    """

    properties: typing.Mapping[str, typing.Any] = dict(
        )
    """Information about the solver. Empty.

    Examples:

        >>> from dwave.samplers import RandomSampler
        >>> sampler = RandomSampler()
        >>> sampler.properties
        {}

    """

    def sample(self,
               bqm: dimod.BinaryQuadraticModel,
               *,
               num_reads: typing.Optional[int] = None,
               time_limit: typing.Optional[typing.Union[float, datetime.timedelta]] = None,
               max_num_samples: int = 1000,
               seed: typing.Union[None, int, np.random.Generator] = None,
               **kwargs,
               ) -> dimod.SampleSet:
        """Return random samples for a binary quadratic model.

        Args:
            bqm: Binary quadratic model to be sampled from.

            num_reads:
                The maximum number of random samples to be drawn.
                If neither ``num_reads`` nor ``time_limit`` are provided,
                ``num_reads`` is set to 1.
                If ``time_limit`` is provided, there is no limit imposed on
                the number of reads.

            time_limit:
                The maximum run time in seconds.
                Only the time to generate the samples, calculate the energies,
                and maintain the population is counted.

            max_num_samples:
                The maximum number of samples returned by the sampler.
                This limits the memory usage for small problems are large
                ``time_limits``.
                Ignored when ``num_reads`` is set.

            seed:
                Seed for the random number generator.
                Passed to :func:`numpy.random.default_rng()`.

        Returns:
            A sample set.
            Some additional information is provided in the
            :attr:`~dimod.SampleSet.info` dictionary:

                * **num_reads**: The total number of samples generated.

        """

        # we could count this towards preprocesing time but IMO it's fine to
        # skip for simplicity.
        self.remove_unknown_kwargs(**kwargs)

        # default case
        if num_reads is None and time_limit is None:
            num_reads = 1
            time_limit = float('inf')

        if num_reads is None:
            # we know that time_limit was set
            num_reads = sys.maxsize
        elif num_reads <= 0:
            raise ValueError("if given, num_reads must be a positive integer")
        else:
            # it was given, so max_num_samples is ignored
            max_num_samples = num_reads

        if time_limit is None:
            # num_reads is specified
            time_limit = float('inf')
        elif isinstance(time_limit, datetime.timedelta):
            time_limit = time_limit.total_seconds()
        elif time_limit <= 0:
            raise ValueError("if given, time_limit must be positive")

        if max_num_samples <= 0:
            raise ValueError("max_num_samples must be a positive integer")

        return sample(bqm,
                      num_reads=num_reads,
                      time_limit=time_limit,
                      max_num_samples=max_num_samples,
                      seed=seed,
                      )
