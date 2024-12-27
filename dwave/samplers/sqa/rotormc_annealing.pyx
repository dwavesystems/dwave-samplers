# distutils: language = c++
# distutils: include_dirs = dwave/samplers/sqa/src/
# distutils: sources = dwave/samplers/sqa/src/cpu_rotormc.cpp
# cython: language_level = 3

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

from libcpp cimport bool
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np


cdef extern from "cpu_rotormc.h":
    ctypedef bool (*callback)(void *function)
    ctypedef enum Proposal:
        GibbsNonErgodic, MetropolisNonErgodic,
        MetropolisUniform, MetropolisTF
    
    int general_simulated_annealing(
            np.uint8_t* samples,
            double* energies,
            const int num_samples,
            const vector[double] & h,
            const vector[int] & coupler_starts,
            const vector[int] & coupler_ends,
            const vector[double] & coupler_weights,
            const vector[double] & trans_fields,
            const int sweeps_per_beta,
            const vector[double] & Hp_field,
            const vector[double] & Hd_field,
            const unsigned long long seed,
            const bool randomize_order,
            const Proposal proposal_acceptance_criteria,
            np.uint8_t *statistics,
            const int schedule_sample_interval,
            callback interrupt_callback,
            void *interrupt_function) nogil


def simulated_annealing(num_samples, h, coupler_starts, coupler_ends,
                        coupler_weights, trans_fields, sweeps_per_beta,
                        Hp_field, Hd_field,
                        seed,
                        np.ndarray[np.uint8_t, ndim=2, mode="c"] states_numpy,
                        randomize_order,
                        proposal_acceptance_criteria,
                        schedule_sample_interval,
                        interrupt_function=None):
    """Wraps `general_simulated_annealing` from `cpu_rotormc.cpp`. Accepts
    an Ising problem defined on a general graph and returns samples
    using simulated annealing of rotor states.

    Parameters
    ----------
    num_samples : int
        Number of samples to get from the sampler.

    h : list(float)
        The h or field values for the problem.

    coupler_starts : list(int)
        A list of the start variable of each coupler. For a problem
        with the couplers (0, 1), (1, 2), and (3, 1), `coupler_starts`
        should be [0, 1, 3].

    coupler_ends : list(int)
        A list of the end variable of each coupler. For a problem
        with the couplers (0, 1), (1, 2), and (3, 1), `coupler_ends`
        should be [1, 2, 1].

    coupler_weights : list(float)
        A list of the J values or weight on each coupler, in the same
        order as `coupler_starts` and `coupler_ends`.

    trans_fields : list(float)
        The transverse field values for the problem.

    sweeps_per_beta : int
        The number of sweeps to perform at each beta value provided in
        `Hp_field`. The total number of sweeps per sample is
        sweeps_per_beta * len(Hp_field).

    Hp_field : list(float)
        A list of the different longitudinal (problem) field values to run sweeps at.

    Hd_field : list(float)
        A list of the different transverse (driver) field values to run sweeps at.

    chain_coupler_strength : float
        A coupling strength applicable to all chains. 
        TO DO: relax to variable (per chain and/or within chain)

    qubits_per_chain : int
        The number of qubits in each chain. Qubits composing each chain must be 
        sequentially labeled. All chains must be equal length and coupled as a
        path. The couplers composing the chain should not be contained in coupler 
        list information. 
        TO DO: relax to any disjoint list.
        
    qubits_per_update : int
        When set 1, qubits are updated in a uniform random order.
        When set to any other value, chains are updated sequentially on each sweep.
        TO DO: better name and arg. checking. perhaps logical_updates =True/False

    seed : 32 bit unsigned int > 0
        The seed to use for the PRNG. Must be a positive integer
        greater than 0. If the same seed is used and the rest of the
        parameters are the same, the returned samples will be
        identical.

    states_numpy : np.ndarray[np.uint8_t, ndim=2, mode="c"], values in (-1, 1)
        The initial seeded states of the simulated annealing runs. Should be of
        a contiguous numpy.ndarray of shape (num_samples, num_variables).

    schedule_sample_interval: int
        Number of schedule changes (steps in Hd, Hp) between 
        samples. If the schedule index mod schedule_sample_interval is 0 then
        sampling is performed after ``num_sweeps_per_beta`` updates.
        Full angle states are stored in uint8 format. Projections can
        be performed subsequently in python when required (see sampler.py).

    interrupt_function: function
        Should accept no arguments and return a bool. The function is
        called between samples and if it returns True, simulated annealing
        will return early with the samples it already has.

    Returns
    -------
    samples : numpy.ndarray
        A 2D numpy array where each row is a sample.

    energies: np.ndarray
        The energies.

    statistics: np.ndarray
        Matrix of statistics:
        [independent processes x MCMC steps tracked]
        by [vector of statistics].
        For now, statistics are mid-anneal samples only, and so the
        type is np.uint8_t for compression purposes.

    """
    num_vars = len(h)

    # in the case that we either need no samples or there are no variables,
    # we can safely return an empty array (and set energies to 0)
    if num_samples*num_vars == 0:
        annealed_states = np.empty((num_samples, num_vars), dtype=np.uint8)
        stats = np.empty(0, dtype=np.uint8) 
        return annealed_states, np.zeros(num_samples, dtype=np.double), stats

    # allocate ndarray for energies
    energies_numpy = np.empty(num_samples, dtype=np.float64)
    cdef double[:] energies = energies_numpy
    
    # Possible feature enhancement: allow alternative statistics,
    # log-intervals, or accumulation; but difficult to accomodate all possible
    # experiments
    num_collection_points = 0
    if schedule_sample_interval:
        num_statistics = num_vars
        num_collection_points = 1 + (Hd_field.size - 1)//schedule_sample_interval
        stat_size = (num_samples, num_collection_points, num_statistics)
    else:
        stat_size = (num_samples, 1, 1) # Non-empty owing to addressing issue..
    cdef np.ndarray[np.uint8_t, ndim=3, mode='c'] statistics_numpy = np.empty(stat_size, dtype=np.uint8)
    # explicitly convert all Python types to C while we have the GIL
    
    cdef np.uint8_t* _states = &states_numpy[0, 0]
    cdef double* _energies = &energies[0]
    cdef np.uint8_t* _statistics = &statistics_numpy[0, 0, 0]
    cdef int _num_samples = num_samples
    cdef vector[double] _h = h
    cdef vector[int] _coupler_starts = coupler_starts
    cdef vector[int] _coupler_ends = coupler_ends
    cdef vector[double] _coupler_weights = coupler_weights
    cdef vector[double] _trans_fields = trans_fields
    cdef int _sweeps_per_beta = sweeps_per_beta
    cdef vector[double] _Hp_field = Hp_field
    cdef vector[double] _Hd_field = Hd_field
    cdef unsigned long long _seed = seed
    cdef bool _randomize_order = randomize_order
    cdef Proposal _proposal
    if proposal_acceptance_criteria == 'GibbsNonErgodic':
        _proposal = GibbsNonErgodic
    elif proposal_acceptance_criteria == 'MetropolisNonErgodic':
        _proposal = MetropolisNonErgodic
    elif proposal_acceptance_criteria == 'MetropolisUniform':
        _proposal = MetropolisUniform
    elif proposal_acceptance_criteria == 'MetropolisTF':
        _proposal = MetropolisTF
    else:
        raise ValueError('Unknown proposal_acceptance_criteria')
    cdef int _schedule_sample_interval = schedule_sample_interval
    cdef void* _interrupt_function
    if interrupt_function is None:
        _interrupt_function = NULL
    else:
        _interrupt_function = <void *>interrupt_function

    with nogil:
        num = general_simulated_annealing(_states,
                                          _energies,
                                          _num_samples,
                                          _h,
                                          _coupler_starts,
                                          _coupler_ends,
                                          _coupler_weights,
                                          _trans_fields,
                                          _sweeps_per_beta,
                                          _Hp_field,
                                          _Hd_field,
                                          _seed,
                                          _randomize_order,
                                          _proposal,
                                          _statistics,
                                          _schedule_sample_interval,
                                          interrupt_callback,
                                          _interrupt_function)

    # discard the noise if we were interrupted
    return states_numpy[:num], energies_numpy[:num], statistics_numpy[:num]

cdef bool interrupt_callback(void * const interrupt_function) noexcept with gil:
    try:
        return (<object>interrupt_function)()
    except Exception:
        # if an exception occurs, treat as an interrupt
        return True
