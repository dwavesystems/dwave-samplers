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

        Get the 20 best random samples found in .2 seconds of searching.

        >>> sampleset = sampler.sample(bqm, num_reads=20, time_limit=.2)

    """

    parameters: typing.Mapping[str, typing.List] = dict(
        num_reads=[],
        seed=[],
        time_limit=[],
        )
    """Keyword arguments accepted by the sampling methods.

    Examples:

        >>> from dwave.samplers import RandomSampler
        >>> sampler = RandomSampler()
        >>> sampler.parameters
        {'num_reads': [], 'seed': [], 'time_limit': []}

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
               num_reads: int = 10,
               seed: typing.Union[None, int, np.random.Generator] = None,
               time_limit: typing.Optional[float] = None,
               **kwargs,
               ) -> dimod.SampleSet:
        """Return random samples for a binary quadratic model.

        Args:
            bqm: Binary quadratic model to be sampled from.

            num_reads: The number of samples to be returned.

            seed:
                Seed for the random number generator.
                Passed to :func:`numpy.random.default_rng()`.

            time_limit:
                The maximum sampling time in seconds.
                If given and non-negative, samples are drawn until ``time_limit``.
                Only the best ``num_reads`` (or fewer) samples are kept.

        Returns:
            A sample set.
            Some additional information is provided in the
            :attr:`~dimod.SampleSet.info` dictionary:

                * **num_drawn**: The total number of samples generated.
                * **prepreocessing_time**: The time to parse the ``bqm`` and to
                  initialize the random number generator.
                * **sampling_time**: The time used to generate the samples
                  and calculate the energies. This is the number controlled by
                  ``time_limit``.
                * **postprocessing_time**: The time to construct the sample
                  set.

        """

        # we could count this towards preprocesing time but IMO it's fine to
        # skip for simplicity.
        self.remove_unknown_kwargs(**kwargs)

        return sample(bqm,
                      num_reads=num_reads,
                      seed=seed,
                      time_limit=-1 if time_limit is None else time_limit,
                      )
