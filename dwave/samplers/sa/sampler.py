# Copyright 2018 D-Wave Systems Inc.
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

import math

from numbers import Integral
from numpy.random import randint
from collections import defaultdict
from typing import List, Sequence, Tuple, Optional, Union
from time import perf_counter_ns
try:
    from typing import Literal
    BetaScheduleType = Literal['linear', 'geometric', 'custom']
except ImportError:
    BetaScheduleType = str


from dimod.core.initialized import InitialStateGenerator

import dimod
import numpy as np

from dwave.samplers.sa.simulated_annealing import simulated_annealing

import warnings

__all__ = ["Neal", "SimulatedAnnealingSampler"]


class SimulatedAnnealingSampler(dimod.Sampler, dimod.Initialized):
    """Simulated annealing sampler for binary quadratic models.

    `Simulated annealing <https://en.wikipedia.org/wiki/Simulated_annealing>`_
    can be used for heuristic optimization or approximate Boltzmann sampling. This
    implementation approaches the equilibrium distribution by performing updates
    at a sequence of decreasing temperatures, terminating at the target
    :math:`\\beta`.\ [#]_ By default each spin is updated once in a fixed order per point
    per :math:`\\beta` according to a Metropolis-Hastings update. When :math:`\\beta`
    is large the target distribution concentrates, at equilibrium, over ground
    states of the model. Samples are guaranteed to match the equilibrium for
    long, smooth :math:`\\beta` schedules.

    .. [#] :math:`\\beta` represents the inverse temperature, :math:`1/(k_B T)`,
       of a `Boltzmann distribution <https://en.wikipedia.org/wiki/Boltzmann_distribution>`_
       where :math:`T` is the thermodynamic temperature in kelvin and :math:`k_B`
       is Boltzmann's constant.

    For more information, see Kirkpatrick, S.; Gelatt Jr, C. D.; Vecchi, M. P.
    (1983). "Optimization by Simulated Annealing". Science. 220 (4598): 671–680

    Also aliased as :class:`.Neal`.

    .. note::
        Ocean SDK versions earlier than 8.0 supported also a
        ``SimulatedAnnealingSampler`` under the ``neal`` namespace in a
        ``dwave-neal`` package. If your code uses that obsoleted class,
        upgrade lines such as

        >>> from neal import SimulatedAnnealingSampler   # doctest: +SKIP

        to

        >>> from dwave.samplers import SimulatedAnnealingSampler

        or

        >>> from dwave.samplers.sa import SimulatedAnnealingSampler

    Examples:
        This example solves a simple Ising problem.

        >>> from dwave.samplers import SimulatedAnnealingSampler
        >>> sampler = SimulatedAnnealingSampler()
        >>> h = {'a': 0.0, 'b': 0.0, 'c': 0.0}
        >>> J = {('a', 'b'): 1.0, ('b', 'c'): 1.0, ('a', 'c'): 1.0}
        >>> sampleset = sampler.sample_ising(h, J, num_reads=10)
        >>> print(sampleset.first.energy)
        -1.0

    """

    parameters = None
    """Keyword arguments accepted by the sampling methods.

    See :meth:`.SimulatedAnnealingSampler.sample` for a description of the
    parameters.

    Examples:
        This example looks at a sampler's parameters and some of their values.

        >>> from dwave.samplers import SimulatedAnnealingSampler
        >>> sampler = SimulatedAnnealingSampler()
        >>> for kwarg in sorted(sampler.parameters):
        ...     print(kwarg)
        beta_range
        beta_schedule
        beta_schedule_type
        initial_states
        initial_states_generator
        interrupt_function
        num_reads
        num_sweeps
        num_sweeps_per_beta
        proposal_acceptance_criteria
        randomize_order
        seed
        >>> sampler.parameters['beta_range']
        []
        >>> sampler.parameters['beta_schedule_type']
        ['beta_schedule_options']

    """

    properties = None
    """A dict containing any additional information about the sampler.

    Examples:
        This example looks at the supported values of a sampler property.

        >>> from dwave.samplers import SimulatedAnnealingSampler
        >>> sampler = SimulatedAnnealingSampler()
        >>> sampler.properties['beta_schedule_options']
        ('linear', 'geometric', 'custom')

    """

    def __init__(self):
        # create a local copy in case folks for some reason want to modify them
        self.parameters = {'beta_range': [],
                           'num_reads': [],
                           'num_sweeps': [],
                           'num_sweeps_per_beta': [],
                           'beta_schedule_type': ['beta_schedule_options'],
                           'beta_schedule': [],
                           'seed': [],
                           'interrupt_function': [],
                           'initial_states': [],
                           'initial_states_generator': [],
                           'randomize_order': [],
                           'proposal_acceptance_criteria': [],
                           }
        self.properties = {'beta_schedule_options': ('linear', 'geometric',
                                                     'custom')}

    def sample(self, bqm: dimod.BinaryQuadraticModel,
               beta_range: Optional[Union[List[float], Tuple[float, float]]] = None,
               num_reads: Optional[int] = None,
               num_sweeps: Optional[int] = None,
               num_sweeps_per_beta: int = 1,
               beta_schedule_type: BetaScheduleType = "geometric",
               seed: Optional[int] = None,
               interrupt_function = None,
               beta_schedule: Optional[Union[Sequence[float], np.ndarray]] = None,
               initial_states: Optional[dimod.typing.SamplesLike] = None,
               initial_states_generator: InitialStateGenerator = "random",
               randomize_order: bool = False,
               proposal_acceptance_criteria: str = 'Metropolis',
               **kwargs) -> dimod.SampleSet:
        """Sample from a binary quadratic model.

        Args:
            bqm: Binary quadratic model to be sampled.

            beta_range:
                A 2-tuple or list defining the beginning and end of the
                :math:`\\beta`\ [#]_ schedule. The schedule is
                interpolated within this range according to the value specified
                by ``beta_schedule_type``. Default range is set based on the
                total bias associated with each node.

            num_reads:
                Number of reads. Each read is generated by one run of the
                simulated annealing algorithm. If ``num_reads`` is not explicitly
                given, it is selected to match the number of initial states
                given. If initial states are not provided, only one read is
                performed.

            num_sweeps:
                Number of sweeps used in annealing. If no value is provided
                and ``beta_schedule`` is None, the value defaults to 1000.

            num_sweeps_per_beta (int, optional, default=1)
                Number of sweeps to perform at each :math:`\\beta`. One sweep
                consists of a sequential Metropolis update of all spins.

            beta_schedule_type:
                :math:`\\beta` schedule type, or how the :math:`\\beta` values
                are interpolated  between the given ``beta_range``. Supported
                values are:

                * "linear"

                * "geometric"

                * "custom"

                "custom" is recommended for high-performance applications, which
                typically require optimizing :math:`\\beta` schedules beyond those
                of the "linear" and "geometric" options, with bounds beyond those
                provided by default. ``num_sweeps_per_beta`` and
                ``beta_schedule`` fully specify a custom schedule.

            beta_schedule:
                Sequence of :math:`\\beta` values swept. Format must be compatible
                with ``numpy.array(beta_schedule, dtype=float)``. Values should
                be non-negative.

            seed:
                Seed to use for the PRNG. Specifying a particular seed with a
                constant set of parameters produces identical results. If not
                provided, a random seed is chosen.

            initial_states:
                One or more samples, each defining an initial state for all the
                problem variables. Initial states are given one per read, but
                if fewer than ``num_reads`` initial states are defined,
                additional values are generated as specified by
                ``initial_states_generator``. See :func:`~dimod.as_samples` for a
                description of "samples-like".

            initial_states_generator:
                Defines the expansion of ``initial_states`` if fewer than
                ``num_reads`` are specified:

                * "none":
                    If the number of initial states specified is smaller than
                    ``num_reads``, raises ``ValueError``.

                * "tile":
                    Reuses the specified initial states if fewer than
                    ``num_reads`` or truncates if greater.

                * "random":
                    Expands the specified initial states with randomly generated
                    states if fewer than ``num_reads`` or truncates if greater.

            randomize_order:
                When `True`, each spin update selects a variable uniformly at random.
                This method is ergodic, obeys detailed balance and preserves symmetries
                of the model.

                When `False`, updates proceed sequentially through the labeled variables
                on each sweep so that all variables are updated once per sweep. This method:

                * can be non-ergodic in special cases when used with ``proposal_acceptance_critera=="Metropolis"``.

                * can introduce a dynamical bias as a function of variable order.

                * has faster per spin update than the True method.

            proposal_acceptance_criteria:
                When "Gibbs", each spin flip proposal is accepted according
                to the Gibbs criteria.
                When "Metropolis", each spin flip proposal is accepted according
                to the Metropolis-Hastings criteria.

            interrupt_function (function, optional):
                A function called with no parameters between each sample of
                simulated annealing. If the function returns True, simulated
                annealing terminates and returns with all of the samples and
                energies found so far.

        Returns:
            A `dimod.SampleSet` for the binary quadratic model.

            The `info` field of the sample set contains information about the sampling procedure:
            1. the beta range used,
            2. the beta schedule type used, and
            3. timing information (details below).

            Timing information is categorized into three: preprocessing, sampling, and
            postprocessing time. All timings are reported in units of nanoseconds.

            Preprocessing time is the total time spent converting the BQM variable type (if
            required), parsing input arguments, and determining an annealing schedule. Sampling time
            is the total time the algorithm spent in sampling states of the binary quadratic model.
            Postprocessing time is the total time spent reverting variable type and creating a
            `dimod.SampleSet`.

        Examples:
            This example runs simulated annealing on a binary quadratic model
            with various input parameters.

            >>> import dimod
            >>> from dwave.samplers import SimulatedAnnealingSampler
            ...
            >>> sampler = SimulatedAnnealingSampler()
            >>> bqm = dimod.BinaryQuadraticModel({'a': .5, 'b': -.5},
            ...                                  {('a', 'b'): -1}, 0.0,
            ...                                  dimod.SPIN)
            >>> # Run with default parameters
            >>> sampleset = sampler.sample(bqm)
            >>> # Run with specified parameters
            >>> sampleset = sampler.sample(bqm, seed=1234,
            ...                            beta_range=[0.1, 4.2],
            ...                            num_sweeps=20,
            ...                            beta_schedule_type='geometric')
            >>> # Reuse a seed
            >>> a1 = next((sampler.sample(bqm, seed=88)).samples())['a']
            >>> a2 = next((sampler.sample(bqm, seed=88)).samples())['a']
            >>> assert a1 == a2

            .. [#] :math:`\\beta` represents the inverse temperature, :math:`1/(k_B T)`, of a
               `Boltzmann distribution <https://en.wikipedia.org/wiki/Boltzmann_distribution>`_
               where :math:`T` is the thermodynamic temperature in kelvin and :math:`k_B` is
               Boltzmann's constant.

        """
        timestamp_preprocess = perf_counter_ns()
        # get the original vartype so we can return consistently
        original_vartype = bqm.vartype

        # convert to spin (if needed)
        if bqm.vartype is not dimod.SPIN:
            bqm = bqm.change_vartype(dimod.SPIN, inplace=False)

        # parse_initial_states could handle seed generation, but because we're
        # sharing it with the SA algo, we handle it out here
        if seed is None:
            seed = randint(2**31)
        elif not isinstance(seed, Integral):
            error_msg = ("'seed' should be None or an integer between 0 and 2^32 "
                         "- 1: value = {}".format(seed))
            raise TypeError(error_msg)
        elif not (0 <= seed < 2**31):
            error_msg = ("'seed' should be an integer between 0 and 2^32 - 1: "
                         "value = {}".format(seed))
            raise ValueError(error_msg)

        # parse the inputs
        parsed = self.parse_initial_states(
            bqm,
            num_reads=num_reads,
            initial_states=initial_states,
            initial_states_generator=initial_states_generator,
            seed=seed)

        num_reads = parsed.num_reads

        # read out the initial states and the variable order
        initial_states_array = np.ascontiguousarray(
            parsed.initial_states.record.sample)

        variable_order = parsed.initial_states.variables

        # read out the BQM
        ldata, (irow, icol, qdata), off = bqm.to_numpy_vectors(
            variable_order=variable_order)

        if interrupt_function and not callable(interrupt_function):
            raise TypeError("'interrupt_function' should be a callable")

        if not isinstance(num_sweeps_per_beta, Integral):
            error_msg = "'num_sweeps_per_beta' should be a positive integer: value = {}".format(num_sweeps_per_beta)
            raise TypeError(error_msg)
        if num_sweeps_per_beta < 1:
            error_msg = "'num_sweeps_per_beta' should be a positive integer: value = {}".format(num_sweeps_per_beta)
            raise ValueError(error_msg)

        # handle beta_schedule et al
        if beta_schedule_type == "custom":


            if beta_schedule is None:
                error_msg = "'beta_schedule' must be provided for beta_schedule_type = 'custom': value is None"
                raise ValueError(error_msg)
            else:
                try:
                    beta_schedule = np.array(beta_schedule,dtype=float)
                except:
                    raise ValueError('beta_schedule cannot be case as a numpy array of dtype=float')
                if num_sweeps is not None and num_sweeps!=len(beta_schedule)*num_sweeps_per_beta:
                    error_msg = "'num_sweeps' should be set to None, or a value consistent with 'beta_schedule' and 'num_sweeps_per_beta' for 'beta_schedule_type' = 'custom': value = ".format(num_sweeps)
                    raise ValueError(error_msg)
                if beta_range is not None and (beta_range[0]!=beta_schedule[0] or beta_range[-1]!=beta_schedule[-1]):
                    error_msg = "'beta_range' should be set to None, or a value consistent with 'beta_schedule', for 'beta_schedule_type'='custom'."
                    raise ValueError(error_msg)
                if np.min(beta_schedule)<0:
                    error_msg = "'beta_schedule' cannot include negative values."
                    raise ValueError(error_msg)
        else:
            if beta_schedule is not None:
                error_msg = "'beta_schedule' must be set to None for 'beta_schedule_type' not equal to 'custom'"
                raise ValueError(error_msg)
            if num_sweeps is None:
                num_sweeps = 1000

            num_betas, rem = divmod(num_sweeps,num_sweeps_per_beta)

            if rem > 0 or num_betas < 0:
                error_msg = "'num_sweeps' must be a non-negative value divisible by 'num_sweeps_per_beta'."
                raise ValueError(error_msg)

            if beta_range is None:
                beta_range = _default_ising_beta_range(bqm.linear, bqm.quadratic)
            elif len(beta_range) != 2 or min(beta_range) < 0:
                error_msg = "'beta_range' should be a 2-tuple, or 2 element list of positive numbers. The latter value is the target value."
                raise ValueError(error_msg)

            if num_betas == 1:
                #One sweep in the target model
                beta_schedule = np.array([beta_range[-1]],dtype=float)
            else:
                if beta_schedule_type == "linear":
                    # interpolate a linear beta schedule
                    beta_schedule = np.linspace(*beta_range, num=num_betas)
                elif beta_schedule_type == "geometric":
                    if min(beta_range) <= 0:
                        error_msg = "'beta_range' must contain non-zero values for 'beta_schedule_type' = 'geometric'"
                        raise ValueError(error_msg)
                    # interpolate a geometric beta schedule
                    beta_schedule = np.geomspace(*beta_range, num=num_betas)
                else:
                    raise ValueError("Beta schedule type {} not implemented".format(beta_schedule_type))

        timestamp_sample = perf_counter_ns()

        # run the simulated annealing algorithm
        samples, energies = simulated_annealing(
            num_reads, ldata, irow, icol, qdata,
            num_sweeps_per_beta, beta_schedule,
            seed, initial_states_array,
            randomize_order, proposal_acceptance_criteria,
            interrupt_function)
        timestamp_postprocess = perf_counter_ns()

        info = {
            "beta_range": beta_range,
            "beta_schedule_type": beta_schedule_type
        }
        response = dimod.SampleSet.from_samples(
            (samples, variable_order),
            energy=energies+bqm.offset,  # add back in the offset
            info=info,
            vartype=dimod.SPIN
        )

        response.change_vartype(original_vartype, inplace=True)

        # Developer note: the specific keys of the timing dict are chosen to be consistent with
        #                 other samplers' timing dict.
        response.info.update(dict(timing=dict(
            preprocessing_ns=timestamp_sample - timestamp_preprocess,
            sampling_ns=timestamp_postprocess - timestamp_sample,
            # Update timing info last to capture the full postprocessing time
            postprocessing_ns=perf_counter_ns() - timestamp_postprocess,
        )))

        return response


Neal = SimulatedAnnealingSampler


def _default_ising_beta_range(h, J,
                              max_single_qubit_excitation_rate = 0.01,
                              scale_T_with_N = True):
    """Determine the starting and ending beta from h J.

    The assumed usage of annealing is optimization. Defaults are chosen
    so that at the hot temperature is fast mixing, and so that at low
    temperature the rate of excitations is small, exploiting a simple
    approximation to the ground state energy space to determine a low
    temperature bound and a worst case approximation to energy
    distribution to determine the high temperature bound.

    Args:
        h (dict):
            External field of Ising model (also called linear bias)

        J (dict):
            Couplings of Ising model (also called quadratic biases)

        max_single_qubit_excitation_rate (float, optional, default = 0.01):
            Targeted single qubit excitation rate at final temperature.
            We can set this value small to lower the probability of
            excitations at the end of the anneal, combined with a simple
            approximation to the ground-state energy structure.
            For standard schedule types (such as geometric), a trade-off
            exists between this threshold and time to solution
            (exploiting multiple restarts).

        scale_T_with_N (bool, optional, default = False):
            The expected number of excitations at finite temperature scales
            with the number of variables. For non-vanishing probability
            of the ground state it is necessary for T to scale to zero
            with N. A simplified approximation to the ground-state energy
            structure is used - for high precision problems T may not scale.

    Assume each variable in J is also in h.

    Generally one should optimize min_beta, max_beta and the shape of the schedule
    to maximize some objective such as time to solution.
    The approximations here are simplified, and can be far from optimal in some cases. Use
    custom schedules to address such cases.
    Approximations used are of complexity O(number of biases), O(1) approximations are also
    possible, as are stronger approximations, but methods are not bottlenecks in practice.
    We use the maximum bias per spin to give an upper bound on the gap, allowing the initial
    sweeps to be fast mixing.
    We use an approximation to the minimum bias per spin to give a lower bound on the minimum
    energy gap, such that for the final sweeps states are highly unlikely to excite away from
    a global (or local) minimum.
    """
    if not 0 < max_single_qubit_excitation_rate < 1:
        raise ValueError('Targeted single qubit excitations rates must be in range (0,1)')

    # Approximate worst and best cases of the [non-zero] energy signal (effective field)
    # experienced per spin as function of neighbors: bias = h_i + sum_j Jij s_j:
    sum_abs_bias_dict = defaultdict(int, {k: abs(v) for k, v in h.items()})
    if sum_abs_bias_dict:
        min_abs_bias_dict = {k: v for k, v in sum_abs_bias_dict.items() if v != 0}
    else:
        min_abs_bias_dict = {}
    #This loop is slow, but is not a bottleneck for practical implementations of simulated annealing:
    for (k1, k2), v in J.items():
        for k in [k1,k2]:
            sum_abs_bias_dict[k] += abs(v)
            if v != 0: #This is slow, but this routine is not a bottleneck:
                if k in min_abs_bias_dict:
                    min_abs_bias_dict[k] = min(abs(v),min_abs_bias_dict[k])
                else:
                    min_abs_bias_dict[k] = abs(v)

    if not min_abs_bias_dict:
        #If dictionary is empty, a Null problem is being presented. This is
        #likely a user error, or edge case test.
        #We should also consider warning for min(min_abs_bias_dict.values())==0.
        #Metropolis-Hastings is not suitable for unbiased and uncoupled
        #variables, sampling of equilibrium is possible, but only if a uniform
        #random initial condition is used (this is default, but allows changes).
        warn_msg = ('All bqm biases are zero (all energies are zero), this is '
                    'likely a value error. Temperature range is set arbitrarily '
                    'to [0.1,1]. Metropolis-Hastings update is non-ergodic.')
        warnings.warn(warn_msg)
        return([0.1,1])


    # Selecting betas based on probability of flipping a qubit
    # Hot temp: We want to be able to quickly search the entire solution space (equilibrate).
    # Scale hot_beta so that all spins are able to flip with probability at least 50%, which
    # ensures mixing across all states is fast.
    # The most unlikely flip is related to the largest energy gap that must be overcome, with
    # Metropolis updates we require:
    #   0.50 = exp(-hot_beta * max_delta_energy) - (1)
    # This is solved as hot_beta = log(2)/max_delta_energy, max_delta energy is twice the
    # effective field, we take a worst case of the effective field to be conservative.
    # Max delta energy occurs when all biases are aligned, and contribute without frustration:
    max_effective_field = max(sum_abs_bias_dict.values(), default=0)

    if max_effective_field == 0:
        hot_beta = 1
    else:
        hot_beta = np.log(2) / (2*max_effective_field)

    # Cold temp: Towards the end of the annealing schedule, we want to eliminate excitations assuming
    # an optimization application; the samples relax to a local (or ideally global) minima.
    # If we arrive at the ground state (or local minima) on the penultimate iteration, we can require
    # that we have low probability to excite on the final iteration.
    # This means that the probability to excite any spin must be bounded by 0.01 on the final sweep
    #   0.01 = sum_i exp(-cold_beta * min_delta_energy_i) - (2)
    # We can approximate min_delta_energy_i as the cost of frustrating the smallest J or h. This approximation
    # can be tight for low precision problems, but can be improved for higher precision ones. In the worst
    # case (high connectity and high precision) determining min_delta_energy_i is NP-hard (related to a
    # knapsack problem) we settle for this approximation, which yields satisfactory results in many standard
    # model types.
    # (2) is a convex optimization problem, and easy to solve with a root finder. However, for brevity we can
    # further simplify, by assuming only spins experiencing minimal gaps are excitable
    #   0.01 ~ #minimal_gaps exp(- cold_beta min_i min_delta_energy_i) - (2)
    # Where #minimal_gaps counts the number of cases i for which min_delta_energy_i = min_i(min_delta_energy_i)
    # i.e. the number of spins which are most easily excited out of the ground state.
    # The solution is cold_beta = log(#minimal_gaps/0.01) / min_delta_energy.
    # min_delta_energy is twice the smallest (non-zero) effective field, we approximate this per spin as the
    # smallest associated bias term (h or J).
    if len(min_abs_bias_dict)==0:
        #Trivial problem
        cold_beta = hot_beta
    else:
        values_array = np.array(list(min_abs_bias_dict.values()),dtype=float)
        min_effective_field = np.min(values_array)
        if scale_T_with_N:
            #This is effective assuming finite precision, for pathological
            #or real valued problems this approximation may perform poorly
            number_min_gaps = np.sum(min_effective_field == values_array)
        else:
            number_min_gaps = 1
        cold_beta = np.log(number_min_gaps/max_single_qubit_excitation_rate) / (2*min_effective_field)

    return [hot_beta, cold_beta]

def default_beta_range(bqm):
    ising = bqm.spin
    return _default_ising_beta_range(ising.linear, ising.quadratic)
