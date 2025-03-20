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
"""
A dimod :term:`sampler` that uses the simulated annealing algorithm over path integral states

To understand the core C++ code, see the original dwave-pimc repository, which was published alongside https://doi.org/10.1038/s41467-021-20901-5
"""
from numbers import Integral
from typing import List, Sequence, Tuple, Optional, Union, Callable
from time import perf_counter_ns

import numpy as np
from numpy.random import randint

import dimod
from dimod.core.initialized import InitialStateGenerator
from dwave.samplers.sqa.pimc_annealing import simulated_annealing as path_annealing
from dwave.samplers.sqa.rotormc_annealing import simulated_annealing as rotor_annealing
from dwave.samplers.sa.sampler import _default_ising_beta_range


__all__ = ["PathIntegralAnnealingSampler", "RotorModelAnnealingSampler"]


class PathIntegralAnnealingSampler(dimod.Sampler, dimod.Initialized):
    """Simulated annealing sampler for path integral representations of BQMs

    A thermalized state of transverse-field Ising models can be understood
    by a distribution on paths (also called worldlines). These worldlines
    can be evolved via classical dynamics. The statistics of a quasi-static
    process can be established by this sampler with a fixed target model
    on a timescale dictated by the mixing time. Some types of quantum
    dynamics, such as incoherent single path tunneling scaling in time,
    can also be emulated.

    Examples:
        This example solves a simple Ising problem.

        >>> from dwave.samplers.sqa import PathIntegralAnnealingSampler
        >>> sampler = PathIntegralAnnealingSampler()
        >>> h = {'a': 0.0, 'b': 0.0, 'c': 0.0}
        >>> J = {('a', 'b'): 1.0, ('b', 'c'): 1.0, ('a', 'c'): 1.0}
        >>> sampleset = sampler.sample_ising(h, J, num_reads=10)
        >>> print(sampleset.first.energy)
        -1.0

    """

    parameters = None
    """dict: A dict where keys are the keyword parameters accepted by the 
    sampler methods (allowed kwargs) and values are lists of 
    :attr:`PathIntegralAnnealingSampler.properties` relevant to each parameter.

    See :meth:`.PathIntegralAnnealingSampler.sample` for a description of the 
    parameters.

    Examples:
        This example looks at a sampler's parameters and some of their values.

        >>> from dwave.samplers.sqa import PathIntegralAnnealingSampler
        >>> sampler = PathIntegralAnnealingSampler()
        >>> for kwarg in sorted(sampler.parameters):
        ...     print(kwarg)
        beta_range
        beta_schedule_type
        initial_states
        initial_states_generator
        interrupt_function
        num_reads
        num_sweeps
        num_sweeps_per_beta
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

        >>> from dwave.samplers.sqa import PathIntegralAnnealingSampler
        >>> sampler = PathIntegralAnnealingSampler()
        >>> sampler.properties['beta_schedule_options']
        ('linear', 'geometric', 'custom')

    """

    def __init__(self):
        # create a local copy in case folks for some reason want to modify them
        self.parameters = {
            "beta_range": [],
            "num_reads": [],
            "num_sweeps": [],
            "num_sweeps_per_beta": [],
            "beta_schedule_type": ["beta_schedule_options"],
            "seed": [],
            "interrupt_function": [],
            "initial_states": [],
            "initial_states_generator": [],
        }
        self.properties = {"beta_schedule_options": ("linear", "geometric", "custom")}

    def sample(
        self,
        bqm: dimod.BinaryQuadraticModel,
        *,
        beta_range: Optional[Union[List[float], Tuple[float, float]]] = None,
        num_reads: Optional[int] = None,
        num_sweeps: Optional[int] = None,
        num_sweeps_per_beta: int = 1,
        beta_schedule_type: str = "geometric",
        seed: Optional[int] = None,
        Hp_field: Optional[Union[Sequence[float], np.ndarray]] = None,
        Hd_field: Optional[Union[Sequence[float], np.ndarray]] = None,
        Gamma: float = 1,
        chain_coupler_strength: float = 1,
        qubits_per_chain: int = 1,
        qubits_per_update: int = 1,
        initial_states: Optional[dimod.typing.SamplesLike] = None,
        initial_states_generator: InitialStateGenerator = "random",
        make_info_json_serializable: bool = False,
        project_states: Tuple[bool, bool] = (True, True),
        num_breaks: Optional[np.ndarray] = None,
        breaks_in: Optional[np.ndarray] = None,
        breaks_buffer_out: Optional[np.ndarray] = None,
        schedule_sample_interval: Optional[int] = None,
        interrupt_function=None,
        **kwargs,
    ) -> dimod.SampleSet:
        """Sample from a binary quadratic model.

        Args:
            bqm:
                The binary quadratic model to be sampled.

            beta_range:
                A 2-tuple or list defining the beginning and end of the beta
                schedule, where beta is the inverse temperature. The schedule is
                interpolated within this range according to the value specified
                by ``beta_schedule_type``. Default range is set based on the
                total bias associated with each node.

            num_reads:
                Number of reads. Each read is generated by one run of the
                simulated annealing algorithm. If `num_reads` is not explicitly
                given, it is selected to match the number of initial states
                given. If initial states are not provided, only one read is
                performed.

            num_sweeps:
                Number of sweeps used in annealing. If no value is provided
                conforms to ``Hp_field`` and ``num_sweeps_per_beta``. If
                ``Hp_field`` is None the value is defaulted to 1000.

            num_sweeps_per_beta:
                Number of sweeps to perform at each beta. One sweep consists of
                a sequential Metropolis update of all spins.

            beta_schedule_type:
                Beta schedule type, or how the beta values are interpolated
                between the given ``beta_range``. Supported values are:

                * "linear"

                * "geometric"

                * "custom"

                "custom" is recommended for high-performance applications, which
                typically require optimizing beta schedules beyond those of the
                "linear" and "geometric" options, with bounds beyond those
                provided by default. ``num_sweeps_per_beta`` and
                ``Hp_field`` fully specify a custom schedule.

            Hp_field:
                Sequence of longitudinal field values swept. Format compatible with
                numpy.array(Hp_field,dtype=float) required. Values should
                be non-negative. Values should be scaled by the target inverse
                temperature.

            Hd_field:
                Sequence of transverse field values swept. Format compatible with
                numpy.array(Hd_field,dtype=float) required. Values should
                be non-negative. Values should be scaled by the target inverse
                temperature.

            Gamma:
                Transverse field applied to all qubits.

            chain_coupler_strength:
                A coupling strength applicable to all chains.
                TO DO: relax to variable (per chain and/or within chain)

            qubits_per_chain:
                The number of qubits in each chain. Qubits composing each chain must be
                sequentially labeled. All chains must be equal length and coupled as a
                path. The couplers composing the chain should not be contained in coupler
                list information.
                TO DO: relax to any disjoint list.

            qubits_per_update:
                When set 1, qubits are updated in a uniform random order.
                When set to any other value, chains are updated sequentially on each sweep.
                TO DO: better name and arg. checking. perhaps logical_updates =True/False

            seed:
                Seed to use for the PRNG. Specifying a particular seed with a
                constant set of parameters produces identical results. If not
                provided, a random seed is chosen.

            initial_states:
                One or more samples, each defining an initial state for all the
                problem variables. Initial states are given one per read, but
                if fewer than ``num_reads`` initial states are defined,
                additional values are generated as specified by
                ``initial_states_generator``. See func:``.as_samples`` for a
                description of "samples-like".

            initial_states_generator:
                Defines the expansion of ``initial_states`` if fewer than
                ``num_reads`` are specified:

                * "none":
                    If the number of initial states specified is smaller than
                    ``num_reads``, raises ValueError.

                * "tile":
                    Reuses the specified initial states if fewer than
                    ``num_reads`` or truncates if greater.

                * "random":
                    Expands the specified initial states with randomly generated
                    states if fewer than ``num_reads`` or truncates if greater.

            make_info_json_serializable:
                For json serialization of returned samplesets.

            project_states:
                First entry indicates if initial_states are classical.
                Second entry indicates if returned states should be projected
                (classical). Single-slice (projected) information is returned as
                samples, supplementary path-integral (break) information is
                returned in the info field.

            num_breaks:
                2d numpy array specifying number of breaks per spin. Variable
                ordering is determined by the sampleset.

            breaks_in:
                1d numpy array storing the ``np.sum(num_breaks)`` break
                positions sequentially over the initiali worldline samples.
                Note: This format may be revisited.

            breaks_buffer_out:
                1d numpy array storing the ``np.sum(num_breaks)`` break
                positions sequentially over the final worldline samples.
                Note: This format may be revisited.

            schedule_sample_interval:
                Number of schedule changes (steps in Hd, Hp) between
                samples. Samples are projected and stored sequentially
                in the sampleset as ``info['statistics']`` after ``num_sweeps_per_beta``
                updates at given schedule (Hd,Hp) values.

            interrupt_function:
                If provided, interrupt_function is called with no parameters
                between each sample of simulated annealing. If the function
                returns True, then simulated annealing will terminate and return
                with all of the samples and energies found so far.

        Returns:
            :obj:``dimod.Response``: A ``dimod`` :obj:``~dimod.Response`` object.

        Examples:
            This example runs path integral annealing on a binary quadratic
            model with various input parameters.

            >>> import dimod
            >>> from dwave.samplers.sqa import PathIntegralAnnealingSampler
            ...
            >>> sampler = PathIntegralAnnealingSampler()
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
            >>> a1 == a2
            True

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
            error_msg = (
                "'seed' should be None or an integer between 0 and 2^32 "
                "- 1: value = {}".format(seed)
            )
            raise TypeError(error_msg)
        elif not (0 <= seed < 2**31):
            error_msg = (
                "'seed' should be an integer between 0 and 2^32 - 1: "
                "value = {}".format(seed)
            )
            raise ValueError(error_msg)

        # parse the inputs
        parsed = self.parse_initial_states(
            bqm,
            num_reads=num_reads,
            initial_states=initial_states,
            initial_states_generator=initial_states_generator,
            seed=seed,
        )

        num_reads = parsed.num_reads

        # read out the initial states and the variable order
        initial_states_array = np.ascontiguousarray(parsed.initial_states.record.sample)

        variable_order = parsed.initial_states.variables

        # read out the BQM
        ldata, (irow, icol, qdata), off = bqm.to_numpy_vectors(
            variable_order=variable_order
        )

        if interrupt_function and not callable(interrupt_function):
            raise TypeError("'interrupt_function' should be a callable")

        if not isinstance(num_sweeps_per_beta, Integral):
            error_msg = (
                "'num_sweeps_per_beta' should be a positive integer: value = {}".format(
                    num_sweeps_per_beta
                )
            )
            raise TypeError(error_msg)
        if num_sweeps_per_beta < 1:
            error_msg = (
                "'num_sweeps_per_beta' should be a positive integer: value = {}".format(
                    num_sweeps_per_beta
                )
            )
            raise ValueError(error_msg)

        # handle beta_schedule et al
        if beta_schedule_type == "custom":

            if Hp_field is None:
                error_msg = "'Hp_field' must be provided for beta_schedule_type = 'custom': value is None"
                raise ValueError(error_msg)
            else:
                try:
                    Hp_field = np.array(Hp_field, dtype=float)
                except:
                    raise ValueError(
                        "Hp_field cannot be case as a numpy array of dtype=float"
                    )
                if (
                    num_sweeps is not None
                    and num_sweeps != len(Hp_field) * num_sweeps_per_beta
                ):
                    error_msg = "'num_sweeps' should be set to None, or a value consistent with 'Hp_field' and 'num_sweeps_per_beta' for 'beta_schedule_type' = 'custom': value = ".format(
                        num_sweeps
                    )
                    raise ValueError(error_msg)
                elif beta_range is not None and (
                    beta_range[0] != Hp_field[0] or beta_range[-1] != Hp_field[-1]
                ):
                    error_msg = "'beta_range' should be set to None, or a value consistent with 'Hp_field', for 'beta_schedule_type'='custom'."
                    raise ValueError(error_msg)
                elif np.min(Hp_field) < 0:
                    error_msg = "'Hp_field' cannot include negative values."
                    raise ValueError(error_msg)
        else:
            if Hp_field is not None:
                error_msg = "'Hp_field' must be set to None for 'beta_schedule_type' not equal to 'custom'"
                raise ValueError(error_msg)
            elif num_sweeps is None:
                num_sweeps = 100

            num_betas, rem = divmod(num_sweeps, num_sweeps_per_beta)

            if rem > 0 or num_betas < 0:
                error_msg = "'num_sweeps' must be a positive value divisible by 'num_sweeps_per_beta'."
                raise ValueError(error_msg)

            if beta_range is None:
                beta_range = _default_ising_beta_range(bqm.linear, bqm.quadratic)
            elif len(beta_range) != 2 or min(beta_range) < 0:
                error_msg = "'beta_range' should be a 2-tuple, or 2 element list of positive numbers. The latter value is the target value."
                raise ValueError(error_msg)

            if num_betas == 1:
                # One sweep in the target model
                Hp_field = np.array([beta_range[-1]], dtype=float)
            else:
                if beta_schedule_type == "linear":
                    # interpolate a linear beta schedule
                    Hp_field = np.linspace(*beta_range, num=num_betas)
                elif beta_schedule_type == "geometric":
                    if min(beta_range) <= 0:
                        error_msg = "'beta_range' must contain non-zero values for 'beta_schedule_type' = 'geometric'"
                        raise ValueError(error_msg)
                    # interpolate a geometric beta schedule
                    Hp_field = np.geomspace(*beta_range, num=num_betas)
                else:
                    raise ValueError(
                        "Beta schedule type {} not implemented".format(
                            beta_schedule_type
                        )
                    )

        num_betas = len(Hp_field)

        if Hd_field is None:
            Hd_field = np.zeros(num_betas)
        else:
            Hd_field = np.array(Hd_field, dtype=float)
            if len(Hd_field) < num_betas:
                raise ValueError("Hd_field incompatible with Hp_field")
        if num_breaks is None:
            num_breaks = np.zeros(shape=initial_states_array.shape, dtype=np.intc)
        elif num_breaks.shape != initial_states_array.shape:
            raise ValueError("Path integral initialization is inconsistent")

        if breaks_in is None:
            breaks_in = np.zeros(max(1, num_breaks.sum()))
        breaks_in = np.ascontiguousarray(breaks_in, dtype=np.intc)

        if breaks_buffer_out is None:
            if not project_states[1]:
                # Approximate sufficient memory requirement:
                # The expected number of breaks (and ~variance thereof) at
                # equilibrium is beta*HdField*tanh(beta*HdField) for a model
                # H = -HdField*sigX. These scale with num_samples and num_reads.
                # The presence of a non-trivial diagonal Hamiltonian suppresses
                # this value in typical case.
                # The buffer size requirement is determined by a complicated
                # extreme value distribution. This is complicated further by
                # variation of Hd when dynamics are at play.
                # The following is a heuristic that may fail. Failure is
                # indicated in the samplest info field by the key
                # ``memory_buffer_failure`` and a warning is raised.
                num_break_patterns = num_reads * bqm.num_variables
                # We want to bound the extreme value with high probability, this
                # lazy heuristic for now (perhaps forever!, logs and const. are cheap):
                if num_betas > 0:
                    out_requirement = int(
                        (1 + Hd_field[-1])
                        * (num_break_patterns + np.log(num_break_patterns))
                        + 1000
                    )
                else:
                    # Do nothing to input:
                    out_requirement = int(np.sum(num_breaks))
            else:
                out_requirement = 1  # bounds check bug
            breaks_buffer_out = np.zeros(out_requirement, dtype=np.intc)
        breaks_buffer_out = np.ascontiguousarray(breaks_buffer_out, dtype=np.intc)
        num_breaks = np.ascontiguousarray(num_breaks, dtype=np.intc)

        if schedule_sample_interval is None:
            schedule_sample_interval = 0  # Flag to indicate no sampling
        elif schedule_sample_interval > len(Hp_field) or schedule_sample_interval < 1:
            raise ValueError("schedule_sample_interval is incompatible with Hp_field")

        timestamp_sample = perf_counter_ns()
        # run the simulated annealing algorithm
        samples, energies, statistics, num_breaks, breaks_buffer_out = path_annealing(
            num_reads,
            ldata,
            irow,
            icol,
            qdata,
            num_sweeps_per_beta,
            Hp_field,
            Hd_field,
            Gamma,
            chain_coupler_strength,
            qubits_per_chain,
            qubits_per_update,
            seed,
            initial_states_array,
            project_states[0],
            project_states[1],
            num_breaks,
            breaks_in,
            breaks_buffer_out,
            schedule_sample_interval,
            interrupt_function,
        )
        timestamp_postprocess = perf_counter_ns()

        info = {
            "beta_range": beta_range,
            "beta_schedule_type": beta_schedule_type,
        }
        if schedule_sample_interval > 0:
            info["statistics"] = statistics

        if not project_states[1]:
            info["num_breaks"] = num_breaks
            num_breaks_total = np.sum(num_breaks)
            if breaks_buffer_out.size >= num_breaks_total:
                info["complete_break_specification"] = True
                info["break_positions"] = breaks_buffer_out[:num_breaks_total]
            else:
                info["complete_break_specification"] = False
                info["break_positions"] = breaks_buffer_out

        if make_info_json_serializable:
            info["beta_range"] = np.array(info["beta_range"]).tolist()
            for numpy_array_key in ["statistics", "num_breaks", "break_positions"]:
                if numpy_array_key in info:
                    info[numpy_array_key] = info[numpy_array_key].tolist()

        response = dimod.SampleSet.from_samples(
            (samples, variable_order),
            energy=energies + bqm.offset,  # add back in the offset
            info=info,
            vartype=dimod.SPIN,
        )

        response.change_vartype(original_vartype, inplace=True)

        response.info.update(
            timing=dict(
                preprocessing_ns=timestamp_sample - timestamp_preprocess,
                sampling_ns=timestamp_postprocess - timestamp_sample,
                # Update timing info last to capture the full postprocessing time
                postprocessing_ns=perf_counter_ns() - timestamp_postprocess,
            )
        )
        return response


class RotorModelAnnealingSampler(dimod.Sampler, dimod.Initialized):
    """Simulated annealing sampler for rotor model representations of BQMs

    A factorized state on coherent spins can be represented by a product
    state on 'rotors'. To reproduce the distribution projected in the
    computational basis one angle is a sufficient description per spin.
    The single state model can be heuristically relaxed to a thermal
    mixture of states in the space of angles. This model can capture
    some equilibrium and some dynamical properties of simple quantum
    annealing processes.
    See https://doi.org/10.48550/arXiv.1401.7087

    Examples:
        This example solves a simple Ising problem.

        >>> from dwave.samplers.sqa import RotorModelAnnealingSampler
        >>> sampler = RotorModelAnnealingSampler()
        >>> h = {'a': 0.0, 'b': 0.0, 'c': 0.0}
        >>> J = {('a', 'b'): 1.0, ('b', 'c'): 1.0, ('a', 'c'): 1.0}
        >>> sampleset = sampler.sample_ising(h, J, num_reads=10)
        >>> print(sampleset.first.energy)
        -1.0

    """

    parameters = None
    """dict: A dict where keys are the keyword parameters accepted by the
    sampler methods (allowed kwargs) and values are lists of
    :attr:`RotorModelAnnealingSampler.properties` relevant to each parameter.

    See :meth:`.RotorModelAnnealingSampler.sample` for a description of the
    parameters.

    Examples:
        This example looks at a sampler's parameters and some of their values.

        >>> from dwave.samplers.sqa import RotorModelAnnealingSampler
        >>> sampler = RotorModelAnnealingSampler()
        >>> for kwarg in sorted(sampler.parameters):
        ...     print(kwarg)
        beta_range
        beta_schedule_type
        initial_states
        initial_states_generator
        interrupt_function
        num_reads
        num_sweeps
        num_sweeps_per_beta
        seed
        >>> sampler.parameters['beta_range']
        []
        >>> sampler.parameters['beta_schedule_type']
        ['beta_schedule_options']

    """

    properties = None
    """dict: A dict containing any additional information about the sampler.

    Examples:
        This example looks at the values set for a sampler property.

        >>> from dwave.samplers.sqa import RotorModelAnnealingSampler
        >>> sampler = RotorModelAnnealingSampler()
        >>> sampler.properties['beta_schedule_options']
        ('linear', 'geometric', 'custom')

    """

    def __init__(self):
        # create a local copy in case folks for some reason want to modify them
        self.parameters = {
            "beta_range": [],
            "num_reads": [],
            "num_sweeps": [],
            "num_sweeps_per_beta": [],
            "beta_schedule_type": ["beta_schedule_options"],
            "seed": [],
            "interrupt_function": [],
            "initial_states": [],
            "initial_states_generator": [],
        }
        self.properties = {"beta_schedule_options": ("linear", "geometric", "custom")}

    def sample(
        self,
        bqm: dimod.BinaryQuadraticModel,
        *,
        beta_range: Optional[Union[List[float], Tuple[float, float]]] = None,
        num_reads: Optional[int] = None,
        num_sweeps: Optional[int] = None,
        num_sweeps_per_beta: int = 1,
        beta_schedule_type: str = "geometric",
        seed: Optional[int] = None,
        Hp_field: Optional[Union[Sequence[float], np.ndarray]] = None,
        Hd_field: Optional[Union[Sequence[float], np.ndarray]] = None,
        initial_states: Optional[Union[dimod.typing.SamplesLike, np.ndarray]] = None,
        initial_states_generator: InitialStateGenerator = "random",
        make_info_json_serializable: bool = False,
        project_states: Tuple[bool, bool] = (False, True),
        schedule_sample_interval: Optional[int] = None,
        randomize_order: bool = False,
        proposal_acceptance_criteria: str = "MetropolisUniform",
        trans_fields: Optional[Union[Sequence[float], np.ndarray]] = None,
        interrupt_function: Optional[Callable[[], bool]] = None,
        **kwargs,
    ) -> dimod.SampleSet:
        """Sample from a binary quadratic model using an implemented sample
        method.

        Args:
            bqm:
                The binary quadratic model to be sampled.

            trans_fields:
                Transverse field applied to all qubits, or a vector.

            beta_range:
                A 2-tuple or list defining the beginning and end of the beta
                schedule, where beta is the inverse temperature. The schedule is
                interpolated within this range according to the value specified
                by ``beta_schedule_type``. Default range is set based on the
                total bias associated with each node.

            num_reads:
                Number of reads. Each read is generated by one run of the
                simulated annealing algorithm. If `num_reads` is not explicitly
                given, it is selected to match the number of initial states
                given. If initial states are not provided, only one read is
                performed.

            num_sweeps:
                Number of sweeps used in annealing. If no value is provided
                and ``Hp_field`` is None the value is defaulted to 1000.

            num_sweeps_per_beta:
                Number of sweeps to perform at each beta. One sweep consists of
                a sequential Metropolis update of all spins.

            beta_schedule_type:
                Beta schedule type, or how the beta values are interpolated
                between the given ``beta_range``. Supported values are:

                * "linear"

                * "geometric"

                * "custom"

                "custom" is recommended for high-performance applications, which
                typically require optimizing beta schedules beyond those of the
                "linear" and "geometric" options, with bounds beyond those
                provided by default. ``num_sweeps_per_beta`` and
                ``Hp_field`` fully specify a custom schedule.

            Hp_field:
                Sequence of longitudinal field values swept. Format compatible with
                numpy.array(Hp_field, dtype=float) required. Values should
                be non-negative. Values should be scaled by the target inverse
                temperature.

            Hd_field:
                Sequence of transverse field values swept. Format compatible with
                numpy.array(Hd_field, dtype=float) required. Values should
                be non-negative. Values should be scaled by the target inverse
                temperature.

            chain_coupler_strength:
                A coupling strength applicable to all chains.
                TO DO: relax to variable (per chain and/or within chain)

            qubits_per_chain:
                The number of qubits in each chain. Qubits composing each chain must be
                sequentially labeled. All chains must be equal length and coupled as a
                path. The couplers composing the chain should not be contained in coupler
                list information.
                TO DO: relax to any disjoint list.

            qubits_per_update:
                When set 1, qubits are updated in a uniform random order.
                When set to any other value, chains are updated sequentially on each sweep.
                TO DO: better name and arg. checking. perhaps logical_updates =True/False

            seed:
                Seed to use for the PRNG. Specifying a particular seed with a
                constant set of parameters produces identical results. If not
                provided, a random seed is chosen.

            initial_states:
                One or more samples, each defining an initial state for all the
                problem variables. Initial states are given one per read, but
                if fewer than ``num_reads`` initial states are defined,
                additional values are generated as specified by
                ``initial_states_generator``. See func:``.as_samples`` for a
                description of "samples-like".

            initial_states_generator:
                Defines the expansion of ``initial_states`` if fewer than
                ``num_reads`` are specified:

                * "none":
                    If the number of initial states specified is smaller than
                    ``num_reads``, raises ValueError.

                * "tile":
                    Reuses the specified initial states if fewer than
                    ``num_reads`` or truncates if greater.

                * "random":
                    Expands the specified initial states with randomly generated
                    states if fewer than ``num_reads`` or truncates if greater.

            project_states:
                First entry indicates if initial_states are classical.
                This first entry is ignored unless the initial_states type is
                uint8_t, the state is assumed to be classical spins in these cases.
                Second entry indicates if full state information should be
                returned. When False the angle-states are returned as an
                an info field, the Sampleset samples are projected regardless.

            make_info_json_serializable:
                For json serialization of returned samplesets.

            randomize_order:
                When True, each spin update selects a variable uniformly at random.
                When False, updates proceed sequentially through the labeled variables
                on each sweep so that all variables are updated once per sweep.
                The True method is ergodic and obeys detailed balance at all temperatures.
                Symmetries of the Boltzmann distribution(s) are not broken by the
                update order.
                The False method
                    - when combined with ``metropolis_update=True`` can be
                    non-ergodic in the limits of zero or infinite temperature,
                    and converge slowly near these limits.
                    - can introduce a dynamical bias as a function of variable
                    labeling convention.
                    - has faster per spin update than the True method.

            proposal_acceptance_criteria:
                Supported dynamics:

                * `MetropolisNonErgodic'
                Angle reflection proposal is accepted by the Metropolis-Hastings criteria  (MHC).

                * `MetropolisUniform`
                A random uniformly distributed angle is proposed and accepted by the MHC [arxiv:1401.70871]

                * 'MetropolisTF'
                A random angle phi, uniformly distributed on [-X, X], is selected, X=min(1, Hd_field/Hp_field))pi.
                The proposal theta_i + X is accepted by the MHC.
                Proposals x>pi are interpretted as 2pi-x, proposals x<0 are interpreted as -x (obeys
                detailed balance and homogenizes rates). [Phys. Rev. Applied 15, 014029]

                * `GibbsNonErgodic`
                Proposals per 'MetropolisNonErgodic', accepted according to Gibbs criteria.

            schedule_sample_interval:
                Number of schedule changes (steps in Hd, Hp) between
                samples. Samples are projected and stored sequentially
                in the sampleset as ``info['statistics']`` after ``num_sweeps_per_beta``
                updates at given schedule (Hd, Hp) values.

            interrupt_function:
                If provided, interrupt_function is called with no parameters
                between each sample of simulated annealing. If the function
                returns True, then simulated annealing will terminate and return
                with all of the samples and energies found so far.

        Returns:
            :obj:``dimod.Response``: A ``dimod`` :obj:``~dimod.Response`` object.

        Examples:
            This example runs simulated annealing on a binary quadratic model
            with some different input parameters.

            >>> import dimod
            >>> from dwave.samplers.sqa import RotorModelAnnealingSampler
            ...
            >>> sampler = RotorModelAnnealingSampler()
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
            >>> a1 == a2
            True

        """
        timestamp_preprocess = perf_counter_ns()

        # get the original vartype so we can return consistently
        original_vartype = bqm.vartype

        # convert to spin (if needed)
        if bqm.vartype is not dimod.SPIN:
            bqm = bqm.change_vartype(dimod.SPIN, inplace=False)

        if trans_fields is None:
            trans_fields = np.ones(bqm.num_variables)
        elif isinstance(trans_fields, (list, tuple, np.ndarray)):
            if len(trans_fields) != bqm.num_variables:
                raise ValueError(
                    "bqm.num_variables and trans_fields lengths incompatible"
                )
            trans_fields = np.array(trans_fields)
        else:
            # Uniform field:
            trans_fields = trans_fields * np.ones(bqm.num_variables)

        # parse_initial_states could handle seed generation, but because we're
        # sharing it with the SA algo, we handle it out here
        if seed is None:
            seed = randint(2**31)
        elif not isinstance(seed, Integral):
            error_msg = (
                "'seed' should be None or an integer between 0 and 2^32 "
                "- 1: value = {}".format(seed)
            )
            raise TypeError(error_msg)
        elif not (0 <= seed < 2**31):
            error_msg = (
                "'seed' should be an integer between 0 and 2^32 - 1: "
                "value = {}".format(seed)
            )
            raise ValueError(error_msg)
        # Required for reproducible projection of samples:
        prng = np.random.default_rng(seed)

        # parse the inputs
        if type(initial_states) == np.ndarray:
            # I don't trust this function yet, ugly compromize
            parsed = self.parse_initial_states(
                bqm,
                num_reads=num_reads,
                initial_states=None,
                initial_states_generator=initial_states_generator,
                seed=seed,
            )
            if num_reads is None:
                num_reads = initial_states.shape[0]
            if (num_reads, bqm.num_variables) != initial_states.shape:
                raise ValueError(
                    "Shape of initial_states specified as"
                    "ndarray should match (num_reads, num_vars)"
                )

        else:
            parsed = self.parse_initial_states(
                bqm,
                num_reads=num_reads,
                initial_states=initial_states,
                initial_states_generator=initial_states_generator,
                seed=seed,
            )
            num_reads = parsed.num_reads
            initial_states = parsed.initial_states.record.sample
        # read out the initial states and the variable order

        variable_order = parsed.initial_states.variables
        if initial_states.dtype == np.uint8:
            if project_states[0]:
                initial_states_array = np.ascontiguousarray(
                    128
                    * (
                        np.cos(np.pi * (initial_states / 128))
                        < 1 - 2 * prng.random(size=initial_states.shape)
                    ),
                    dtype=np.uint8,
                )
            else:
                initial_states_array = np.ascontiguousarray(initial_states)
        else:
            # Assume spins, and cast to uint8_t
            if not set(np.unique(initial_states)).issubset({-1, 1}):
                raise ValueError(
                    "initial_states dtype is neither np.uint8 nor consistent with spins {-1, 1}"
                )
            initial_states_array = np.ascontiguousarray(
                64 - 64 * initial_states, dtype=np.uint8
            )
        # read out the BQM
        ldata, (irow, icol, qdata), off = bqm.to_numpy_vectors(
            variable_order=variable_order
        )

        if interrupt_function and not callable(interrupt_function):
            raise TypeError("'interrupt_function' should be a callable")

        if not isinstance(num_sweeps_per_beta, Integral):
            error_msg = (
                "'num_sweeps_per_beta' should be a positive integer: value = {}".format(
                    num_sweeps_per_beta
                )
            )
            raise TypeError(error_msg)
        elif num_sweeps_per_beta < 1:
            error_msg = (
                "'num_sweeps_per_beta' should be a positive integer: value = {}".format(
                    num_sweeps_per_beta
                )
            )
            raise ValueError(error_msg)
        # handle beta_schedule et al
        if beta_schedule_type == "custom":
            if Hp_field is None:
                error_msg = "'Hp_field' must be provided for beta_schedule_type = 'custom': value is None"
                raise ValueError(error_msg)
            else:
                try:
                    Hp_field = np.array(Hp_field, dtype=float)
                except:
                    raise ValueError(
                        "Hp_field cannot be case as a numpy array of dtype=float"
                    )
                if (
                    num_sweeps is not None
                    and num_sweeps != len(Hp_field) * num_sweeps_per_beta
                ):
                    raise ValueError(
                        "'num_sweeps' should be set to None, or a value consistent with 'Hp_field' "
                        "and 'num_sweeps_per_beta' for 'beta_schedule_type' = 'custom': value = "
                        f"{num_sweeps} != {Hp_field.size}x{num_sweeps_per_beta}"
                    )
                elif beta_range is not None and (
                    beta_range[0] != Hp_field[0] or beta_range[-1] != Hp_field[-1]
                ):
                    raise ValueError(
                        "'beta_range' should be set to None, or a value consistent "
                        "with 'Hp_field', for 'beta_schedule_type'='custom'."
                    )
                elif np.min(Hp_field) < 0:
                    raise ValueError("'Hp_field' cannot include negative values.")
        elif num_sweeps == 0:
            Hp_field = []
        else:
            if Hp_field is not None:
                error_msg = "'Hp_field' must be set to None for 'beta_schedule_type' not equal to 'custom'"
                raise ValueError(error_msg)
            elif num_sweeps is None:
                num_sweeps = 100

            num_betas, rem = divmod(num_sweeps, num_sweeps_per_beta)

            if rem > 0 or num_betas < 0:
                error_msg = "'num_sweeps' must be a positive value divisible by 'num_sweeps_per_beta'."
                raise ValueError(error_msg)

            if beta_range is None:
                beta_range = _default_ising_beta_range(bqm.linear, bqm.quadratic)
            elif len(beta_range) != 2 or min(beta_range) < 0:
                error_msg = "'beta_range' should be a 2-tuple, or 2 element list of positive numbers. The latter value is the target value."
                raise ValueError(error_msg)
            if num_betas == 1:
                # One sweep in the target model
                Hp_field = np.array([beta_range[-1]], dtype=float)
            else:
                if beta_schedule_type == "linear":
                    # interpolate a linear beta schedule
                    Hp_field = np.linspace(*beta_range, num=num_betas)
                elif beta_schedule_type == "geometric":
                    if min(beta_range) <= 0:
                        error_msg = "'beta_range' must contain non-zero values for 'beta_schedule_type' = 'geometric'"
                        raise ValueError(error_msg)
                    # interpolate a geometric beta schedule
                    Hp_field = np.geomspace(*beta_range, num=num_betas)
                else:
                    raise ValueError(
                        "Beta schedule type {} not implemented".format(
                            beta_schedule_type
                        )
                    )

        num_betas = len(Hp_field)

        if Hd_field is None:
            Hd_field = np.zeros(num_betas)
        else:
            Hd_field = np.array(Hd_field, dtype=float)
            if len(Hd_field) < num_betas:
                raise ValueError("Hd_field incompatible with Hp_field")

        if schedule_sample_interval is None:
            schedule_sample_interval = 0  # Flag to indicate no sampling
        elif schedule_sample_interval > len(Hp_field) or schedule_sample_interval < 1:
            raise ValueError("schedule_sample_interval is incompatible with Hp_field")

        timestamp_sample = perf_counter_ns()
        # run the simulated annealing algorithm
        samples, energies, statistics = rotor_annealing(
            num_reads,
            ldata,
            irow,
            icol,
            qdata,
            trans_fields,
            num_sweeps_per_beta,
            Hp_field,
            Hd_field,
            seed,
            initial_states_array,
            randomize_order,
            proposal_acceptance_criteria,
            schedule_sample_interval,
            interrupt_function,
        )
        timestamp_postprocess = perf_counter_ns()

        info = {
            "beta_range": beta_range,
            "beta_schedule_type": beta_schedule_type,
        }
        if schedule_sample_interval > 0:
            info["statistics"] = statistics

        if not project_states[1]:
            info["rotor_states"] = samples.copy()
        if make_info_json_serializable:
            info["beta_range"] = np.array(info["beta_range"]).tolist()
            for numpy_array_key in ["statistics", "rotor_states"]:
                if numpy_array_key in info:
                    info[numpy_array_key] = info[numpy_array_key].tolist()
        # states are always projected to SPINS to play nice(st)
        # with SampleSet and ocean.
        samples = np.array(
            1
            - 2
            * (np.cos(np.pi * samples / 128) < 1 - 2 * prng.random(size=samples.shape)),
            dtype=np.int8,
        )

        response = dimod.SampleSet.from_samples(
            (samples, variable_order),
            energy=energies + bqm.offset,  # add back in the offset
            info=info,
            vartype="SPIN",
        )

        response.change_vartype(original_vartype, inplace=True)

        # Developer note: the specific keys of the timing dict are chosen to be consistent with
        #                 other samplers' timing dict.
        response.info.update(
            dict(
                timing=dict(
                    preprocessing_ns=timestamp_sample - timestamp_preprocess,
                    sampling_ns=timestamp_postprocess - timestamp_sample,
                    # Update timing info last to capture the full postprocessing time
                    postprocessing_ns=perf_counter_ns() - timestamp_postprocess,
                )
            )
        )

        return response
