# cython: language_level = 3

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

from libcpp cimport bool
from libcpp.vector cimport vector

import numpy as np


cdef extern from "cpu_dsbm.h":
    ctypedef bool (*callback)(void *function)
    
    int general_discrete_simulated_bifurcation_machine(
            double* samples_x,
            double* samples_y,
            const int num_samples,
	    const int num_vars,
            const vector[int] & coupler_starts,
            const vector[int] & coupler_ends,
            const vector[double] & coupler_weights,
	    const vector[double] & a_schedule,
            const double a0,
	    const double c0,
	    const double Delta_t,
            callback interrupt_callback,
            void *interrupt_function) nogil


def dsbm(num_samples, num_vars, coupler_starts, coupler_ends,
    coupler_weights, a0, c0, Delta_t,
    a_schedule,
    double[:, ::1] initial_x,
    double[:, ::1] initial_y,
    interrupt_function=None):
    """Wraps `general_discrete_simulated_bifurcation_machine` from `cpu_dsbm.cpp`. Accepts
    an Ising problem defined on a general graph and returns samples
    using discrete_simulated_bifurcation_machine

    Parameters
    ----------
    num_samples : int
        Number of samples to get from the sampler.

    num_vars : int
        Number of variables in the problem

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

    initial_x : np.ndarray[np.uint8_t, ndim=2, mode="c"], values in (-1, 1)
        The initial seeded states of the simulated annealing runs. Should be of
        a contiguous numpy.ndarray of shape (num_samples, num_variables).

    initial_y : np.ndarray[np.uint8_t, ndim=2, mode="c"], values in (-1, 1)
        The initial seeded states of the simulated annealing runs. Should be of
        a contiguous numpy.ndarray of shape (num_samples, num_variables).

    interrupt_function: function
        Should accept no arguments and return a bool. The function is
        called between samples and if it returns True, simulated annealing
        will return early with the samples it already has.

    Returns
    -------
    sample : int
        Number of samples processed

    """
    # in the case that we either need no samples or there are no variables,
    # we can safely return an empty array (and set energies to 0)
    if num_samples*num_vars == 0:
        return num_samples
    
    cdef double* _states_x = &initial_x[0, 0]
    cdef double* _states_y = &initial_y[0, 0]
    cdef int _num_samples = num_samples
    cdef int _num_vars = num_vars
    cdef vector[int] _coupler_starts = coupler_starts
    cdef vector[int] _coupler_ends = coupler_ends
    cdef vector[double] _coupler_weights = coupler_weights
    cdef vector[double] _a_schedule = a_schedule
    cdef float _a0 = a0
    cdef float _c0 = c0
    cdef float _Delta_t = Delta_t
    cdef void* _interrupt_function
    if interrupt_function is None:
        _interrupt_function = NULL
    else:
        _interrupt_function = <void *>interrupt_function

    with nogil:
        num = general_discrete_simulated_bifurcation_machine(_states_x,
						_states_y,	
                                          	_num_samples,	
                                          	_num_vars,	
                                          	_coupler_starts,
                                          	_coupler_ends,
                                          	_coupler_weights,
						_a_schedule,
						_a0,
						_c0,
						_Delta_t,
                                          	interrupt_callback,
                                          	_interrupt_function)

    # Return by reference
    return num

cdef bool interrupt_callback(void * const interrupt_function) noexcept with gil:
    try:
        return (<object>interrupt_function)()
    except Exception:
        # if an exception occurs, treat as an interrupt
        return True
