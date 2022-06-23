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
from libcpp.algorithm cimport make_heap, push_heap, pop_heap
from libcpp.utility cimport pair
from libcpp.vector cimport vector

import dimod
cimport dimod
import numpy as np
cimport numpy as np
cimport numpy.random

# chrono is not included in Cython's libcpp. So we do it more manually
cdef extern from *:
    """
    #include <chrono>

    double realtime_clock() {
        auto t = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(t.time_since_epoch()).count();
    }
    """
    double realtime_clock()


# it would be nicer to use a struct, but OSX throws segfaults when sorting
# Cython-created structs. We should retest that when we switch to Cython 3.
ctypedef pair[np.float64_t, vector[np.int8_t]] state_t


cdef state_t get_sample(dimod.cyBQM_float64 cybqm,
                        numpy.random.bitgen_t* bitgen,
                        bint is_spin = False,
                        ):
    # developer note: there is a bunch of potential optimization here
    cdef state_t state

    # generate the sample
    state.second.reserve(cybqm.num_variables())
    cdef Py_ssize_t i
    for i in range(cybqm.num_variables()):
        state.second.push_back(bitgen.next_uint32(bitgen.state) % 2)

    if is_spin:
        # go back through and convert to spin
        for i in range(state.second.size()):
            state.second[i] = 2 * state.second[i] - 1

    state.first = cybqm.data().energy(state.second.begin())

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

    cdef vector[state_t] samples
    for i in range(num_reads):
        samples.push_back(get_sample(cybqm, bitgen, is_spin))

    cdef Py_ssize_t num_drawn = num_reads

    if time_limit >= 0 and realtime_clock() < sampling_stop_time:
        make_heap(samples.begin(), samples.end())

        while realtime_clock() < sampling_stop_time:
            samples.push_back(get_sample(cybqm, bitgen, is_spin))
            push_heap(samples.begin(), samples.end())
            pop_heap(samples.begin(), samples.end())
            samples.pop_back()

            num_drawn += 1

    cdef double postprocessing_start_time = realtime_clock()

    record = np.rec.array(np.empty(num_reads,
                      dtype=[('sample', np.int8, (cybqm.num_variables(),)),
                             ('energy', float),
                             ('num_occurrences', int)]))

    record['num_occurrences'][:] = 1

    cdef np.float64_t[:] energies_view = record['energy']
    for i in range(num_reads):
        energies_view[i] = samples[i].first

    cdef np.int8_t[:, :] sample_view = record['sample']
    for i in range(num_reads):
        for j in range(cybqm.num_variables()):
            sample_view[i, j] = samples[i].second[j]

    sampleset = dimod.SampleSet(record, bqm.variables, info=dict(), vartype=bqm.vartype)

    sampleset.info.update(
        num_drawn=num_drawn,
        prepreocessing_time=sampling_start_time-preprocessing_start_time,
        sampling_time=postprocessing_start_time-sampling_start_time,
        postprocessing_time=realtime_clock()-preprocessing_start_time,
        )

    return sampleset
