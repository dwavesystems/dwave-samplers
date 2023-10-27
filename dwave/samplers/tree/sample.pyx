# distutils: language = c++
# cython: language_level = 3
# distutils: include_dirs = dwave/samplers/tree/src/include/
#
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

from typing import Tuple
from libcpp cimport bool
from libc.stdlib cimport free


import numpy as np
import dimod
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel

cimport cython
cimport numpy as np
from cython.operator cimport dereference as deref
from dimod cimport cyBQM_float64
from dimod.libcpp cimport BinaryQuadraticModel as cppBinaryQuadraticModel

from dwave.samplers.tree.orang cimport samples_type, PyArray_ENABLEFLAGS

cdef extern from "src/include/sample.hpp":
    void sampleBQM[V, B](cppBinaryQuadraticModel[B, V]& refBQM,
                         int* var_order,
                         double beta,
                         int low,
                         double max_complexity,
                         int num_samples,
                         bool marginals,
                         int seed,
                         double* log_pf,
                         int** samples_data, int* samples_rows, int* samples_cols,
                         double** single_mrg_data, int* single_mrg_len,
                         double** pair_mrg_data, int* pair_mrg_rows, int* pair_mrg_cols,
                         int** pair_data, int* pair_rows, int* pair_cols) except +

@cython.annotation_typing(False)  # maintain Cython2 behavior
def sample_bqm_wrapper(bqm: BinaryQuadraticModel,
                       beta: float,
                       max_complexity: float,
                       order: list,
                       marginals: bool = False,
                       num_reads: int = 1,
                       seed: float = None) -> Tuple[np.ndarray, dict]:
    """Cython wrapper for :func:`sampleBQM`.

    Args:
        bqm:
            Binary quadratic model to sample from. Variables must be linearly indexed.

        beta:
            `Boltzmann distribution`_ inverse temperature parameter.

        max_complexity:
            Upper bound on algorithm's complexity.

        order:
            List of variables representing the variable elimination order.

        marginals:
            Determines whether or not to compute the marginals. If True, they
            will be included in the returned dict.

        num_reads:
            Number of samples to return.

        seed:
            Random number generator seed. Negative values will cause a time-based
            seed to be used.

    Returns:
        The samples and marginals.
    """
    if not bqm.num_variables:
        raise ValueError("bqm must have at least one variable.")

    if len(order) != bqm.num_variables or set(order) != bqm.variables:
        raise ValueError("order must contain the variables in bqm.")

    cdef cyBQM_float64 cybqm = bqm.data

    cdef double _beta = beta
    cdef double _max_complexity = max_complexity
    cdef int low = -1 if bqm.vartype is dimod.SPIN else 0

    cdef int _seed
    if seed is None:
        _seed = np.random.randint(np.iinfo(np.intc).max, dtype=np.intc)
    else:
        _seed = seed

    cdef int[:] elimination_order = np.asarray(order, dtype=np.intc)
    cdef int* elimination_order_ptr = &elimination_order[0]

    # pass in pointers so that sample_bqm can fill things in
    cdef double logpf

    # samples
    cdef int srows, scols
    cdef samples_type* samples_pointer

    # marginals
    cdef double* single_marginals_pointer
    cdef int smlen

    cdef double* pair_marginals_pointer
    cdef int pmrows, pmcols

    cdef int* pair_pointer
    cdef int prows, pcols

    sampleBQM(deref(cybqm.cppbqm),
              elimination_order_ptr,
              _beta,
              low,
              _max_complexity,
              num_reads,
              marginals,
              _seed,
              &logpf,
              &samples_pointer, &srows, &scols,
              &single_marginals_pointer, &smlen,
              &pair_marginals_pointer, &pmrows, &pmcols,
              &pair_pointer, &prows, &pcols)

    # create a numpy array without making a copy then tell numpy it needs to
    # free the memory
    samples = np.asarray(<samples_type[:srows, :scols]> samples_pointer)
    PyArray_ENABLEFLAGS(samples, np.NPY_OWNDATA)

    # convert the samples to spin if necessary
    cdef size_t i, j
    if low == -1:
        for i in range(srows):
            for j in range(scols):
                if samples[i, j] == 0:
                    samples[i, j] = -1

    if marginals:
        variable_marginals = np.asarray(<double[:smlen]> single_marginals_pointer)
        PyArray_ENABLEFLAGS(variable_marginals, np.NPY_OWNDATA)

        if pmrows * pmcols:
            interaction_marginals = np.asarray(<double[:pmrows, :pmcols]> pair_marginals_pointer)
            PyArray_ENABLEFLAGS(interaction_marginals, np.NPY_OWNDATA)
        else:
            interaction_marginals = np.empty(shape=(pmrows, pmcols), dtype=np.double)

        if prows * pcols:
            interactions = np.asarray(<int[:prows, :pcols]> pair_pointer)
            PyArray_ENABLEFLAGS(interactions, np.NPY_OWNDATA)
        else:
            interactions = np.empty(shape=(prows, pcols), dtype=np.double)

        marginal_data = dict(variable_marginals=variable_marginals,
                             interactions=interactions,
                             interaction_marginals=interaction_marginals,
                             log_partition_function=logpf,
                             )
    else:
        marginal_data = dict(log_partition_function=logpf)

    marginal_data['log_partition_function'] += -beta*bqm.offset

    return samples, marginal_data
