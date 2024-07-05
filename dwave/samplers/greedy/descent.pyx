# distutils: language = c++
# distutils: include_dirs = greedy/src/
# distutils: include_dirs = dwave/samplers/greedy/src/
# distutils: sources = dwave/samplers/greedy/src/descent.cpp

# Copyright 2019 D-Wave Systems Inc.
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

from libcpp cimport bool
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

cimport dwave.samplers.greedy.decl as decl

def steepest_gradient_descent(num_samples,
                              linear_biases,
                              coupler_starts, coupler_ends, coupler_weights,
                              np.ndarray[np.int8_t, ndim=2, mode="c"] states_numpy,
                              large_sparse_opt=False):

    """Wraps `steepest_gradient_descent` from `descent.cpp`. Accepts
    an Ising problem defined on a general graph and returns samples
    using steepest gradient descent seeded with initial states from
    `states_numpy`.

    Parameters
    ----------
    num_samples : int
        Number of samples to get from the sampler.

    linear_biases : list(float)
        The linear biases or field values for the problem.

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

    states_numpy : np.ndarray[np.int8_t, ndim=2, mode="c"], values in (-1, 1)
        The initial seeded states of the gradient descent runs. Should be of
        a contiguous numpy.ndarray of shape (num_samples, num_variables).

    large_sparse_opt : bool
        When set to True, large-and-sparse problem graph optimizations are used.

    Returns
    -------
    samples : numpy.ndarray
        A 2D numpy array where each row is a sample.

    energies: numpy.ndarray
        Sample energies.

    num_steps: numpy.ndarray
        Number of downhill steps per sample.
    """
    num_vars = len(linear_biases)

    # short-circuit null edge cases
    if num_samples == 0 or num_vars == 0:
        states = np.empty((num_samples, num_vars), dtype=np.int8)
        return (states,
                np.zeros(num_samples, dtype=np.double),
                np.zeros(num_samples, dtype=np.uint32))

    # allocate ndarray for energies
    energies_numpy = np.empty(num_samples, dtype=np.float64)
    cdef double[:] energies = energies_numpy

    # allocate ndarray for steps
    num_steps_numpy = np.empty(num_samples, dtype=np.uint32)
    cdef unsigned[:] num_steps = num_steps_numpy

    # explicitly convert all Python types to C while we have the GIL
    cdef np.int8_t* _states = &states_numpy[0, 0]
    cdef double* _energies = &energies[0]
    cdef unsigned* _num_steps = &num_steps[0]
    cdef int _num_samples = num_samples
    cdef bool _large_sparse_opt = large_sparse_opt
    cdef vector[double] _linear_biases = linear_biases

    # TODO: in dimod 0.10+, coupler indices default to int64, but we downcast
    # them still to int. We can update this when we upgrade to a more efficient
    # bqm data interface. Until then, beware of the 2B num_variables limit.
    cdef vector[int] _coupler_starts = coupler_starts
    cdef vector[int] _coupler_ends = coupler_ends
    cdef vector[double] _coupler_weights = coupler_weights

    with nogil:
        decl.steepest_gradient_descent(
            _states, _energies, _num_steps, _num_samples,
            _linear_biases, _coupler_starts, _coupler_ends, _coupler_weights,
            _large_sparse_opt)

    return states_numpy, energies_numpy, num_steps_numpy
