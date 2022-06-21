# distutils: language = c++
# cython: language_level = 3

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

cimport cython

from cpython.pycapsule cimport PyCapsule_IsValid, PyCapsule_GetPointer
from libc.time cimport time, time_t
from libcpp.algorithm cimport sort
from libcpp.vector cimport vector
from posix.time cimport clock_gettime, timespec, CLOCK_REALTIME

import dimod
cimport dimod
import numpy as np
cimport numpy as np
cimport numpy.random

cdef extern from *:
    """
    #if defined(_WIN32) || defined(_WIN64)

    #include <Windows.h>

    double realtime_clock() {
        LARGE_INTEGER frequency;
        LARGE_INTEGER now;

        QueryPerformanceFrequency(&frequency);
        QueryPerformanceCounter(&now);

        return now.QuadPart / frequency.QuadPart;
    }

    #else

    double realtime_clock() {
        struct timespec ts;
        clock_gettime(CLOCK_MONOTONIC, &ts);
        return ts.tv_sec + ts.tv_nsec / 1e9;
    }

    #endif
    """
    double realtime_clock()


cdef struct state_t:
    np.float64_t energy
    vector[np.int8_t] sample


cdef bint comp_state(const state_t& a, const state_t& b) nogil:
    return a.energy < b.energy


cdef state_t get_sample(dimod.cyBQM_float64 cybqm,
                        numpy.random.bitgen_t* bitgen,
                        bint is_spin = False,
                        ):
    # developer note: there is a bunch of potential optimization here
    cdef state_t state

    # generate the sample
    state.sample.reserve(cybqm.num_variables())
    cdef Py_ssize_t i
    for i in range(cybqm.num_variables()):
        state.sample.push_back(bitgen.next_uint32(bitgen.state) % 2)

    if is_spin:
        # go back through and convert to spin
        for i in range(state.sample.size()):
            state.sample[i] = 2 * state.sample[i] - 1

    state.energy = cybqm.data().energy(state.sample.begin())

    return state


@cython.boundscheck(False)
@cython.wraparound(False)
def sample(bqm, Py_ssize_t num_reads, object seed, np.float64_t time_limit):

    cdef double preprocessing_start_time = realtime_clock()

    cdef Py_ssize_t i, j  # counters for use later

    # Get Cython access to the BQM. We could template to avoid the copy,
    # but honestly everyone just uses float64 anyway so...
    cdef dimod.cyBQM_float64 cybqm = dimod.as_bqm(bqm, dtype=float).data
    cdef bint is_spin = bqm.vartype is dimod.SPIN

    # Get Cython access to the rng
    rng = np.random.default_rng(seed)
    cdef numpy.random.bitgen_t *bitgen
    cdef const char *capsule_name = "BitGenerator"
    capsule = rng.bit_generator.capsule
    if not PyCapsule_IsValid(capsule, capsule_name):
        raise ValueError("Invalid pointer to anon_func_state")
    bitgen = <numpy.random.bitgen_t *> PyCapsule_GetPointer(capsule, capsule_name)

    cdef double sampling_start_time = realtime_clock()

    cdef double sampling_stop_time
    if time_limit < 0:
        sampling_stop_time = float('inf')
    else:
        sampling_stop_time = sampling_start_time + time_limit

    # try sampling
    cdef Py_ssize_t num_drawn = 0
    cdef vector[state_t] samples
    for i in range(num_reads):
        samples.push_back(get_sample(cybqm, bitgen, is_spin))
        num_drawn += 1

        if realtime_clock() > sampling_stop_time:
            break

    if time_limit >= 0:
        while realtime_clock() < sampling_stop_time:

            samples.push_back(get_sample(cybqm, bitgen, is_spin))
            sort(samples.begin(), samples.end(), comp_state)
            samples.pop_back()

            num_drawn += 1

    cdef double postprocessing_start_time = realtime_clock()

    if time_limit < 0:
         # for consistency we sort in this case as well, though we count
         # it towards postprocessing since it's not necessary
        sort(samples.begin(), samples.end(), comp_state)

    record = np.rec.array(np.empty(num_reads,
                      dtype=[('sample', np.int8, (cybqm.num_variables(),)),
                             ('energy', float),
                             ('num_occurrences', int)]))

    record['num_occurrences'][:] = 1

    cdef np.float64_t[:] energies_view = record['energy']
    for i in range(num_reads):
        energies_view[i] = samples[i].energy

    cdef np.int8_t[:, :] sample_view = record['sample']
    for i in range(num_reads):
        for j in range(cybqm.num_variables()):
            sample_view[i, j] = samples[i].sample[j]

    sampleset = dimod.SampleSet(record, bqm.variables, info=dict(), vartype=bqm.vartype)

    sampleset.info.update(
        num_drawn=num_drawn,
        prepreocessing_time=sampling_start_time-preprocessing_start_time,
        sampling_time=postprocessing_start_time-sampling_start_time,
        postprocessing_time=realtime_clock()-preprocessing_start_time,
        )

    return sampleset
