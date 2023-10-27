# distutils: language = c++
# cython: language_level = 3
# distutils: include_dirs = dwave/samplers/tree/src/include/

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

import dimod
import numpy as np

cimport cython
cimport numpy as np
from cython.operator cimport dereference as deref
from dimod cimport cyBQM_float64
from dimod.libcpp cimport BinaryQuadraticModel as cppBinaryQuadraticModel
from dimod.binary.binary_quadratic_model import BinaryQuadraticModel

from dwave.samplers.tree.orang cimport energies_type, samples_type, PyArray_ENABLEFLAGS

cdef extern from "src/include/solve.hpp":
    void solveBQM[V, B](cppBinaryQuadraticModel[B, V]& refBQM,
                        int* var_order,
                        double beta,
                        int low,
                        double max_complexity,
                        int max_solutions,
                        double** energies_data, int* energies_len,
                        int** sols_data, int* sols_rows, int* sols_cols) except +

samples_dtype = np.intc  # needs to be consistent with samples_type
energies_dtype = np.double  # needs to be consistent with energies_type 

@cython.annotation_typing(False)  # maintain Cython2 behavior
def solve_bqm_wrapper(bqm: BinaryQuadraticModel,
                      order: list,
                      max_complexity: int,
                      max_solutions: int = 1):
    """Cython wrapper for :func:`solveBQM`.

    Args:
        bqm:
            Binary quadratic model to sample from. Variables must be linearly indexed.

        order:
            List of variables representing the variable elimination order.

        max_complexity:
            Upper bound on algorithm's complexity.

        max_solutions:
            Maximum number of solutions to find.

    Returns:
        The samples and marginals.
    """

    if not bqm.num_variables:
        raise ValueError("bqm must have at least one variable.")

    if len(order) != bqm.num_variables or set(order) != bqm.variables:
        raise ValueError("order must contain the variables in bqm.")

    cdef cyBQM_float64 cybqm = bqm.data
    cdef double beta = -1  # solving
    cdef double _max_complexity = max_complexity
    cdef int low = -1 if bqm.vartype is dimod.SPIN else 0

    cdef int[:] elimination_order = np.asarray(order, dtype=np.intc)
    cdef int* elimination_order_ptr = &elimination_order[0]

    # Pass in a pointer so that solveBQM can fill it in. Note that this is
    # a design choice inherited from orang's c++ implementation. In the future
    # we may want to change it.
    cdef int num_energies, srows, scols
    cdef energies_type* energies_pointer
    cdef samples_type* samples_pointer

    solveBQM(deref(cybqm.cppbqm),
             elimination_order_ptr,
             beta,
             low,
             _max_complexity,
             max_solutions,
             &energies_pointer, &num_energies,
             &samples_pointer, &srows, &scols
            )

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

    energies = np.asarray(<energies_type[:num_energies]> energies_pointer)

    PyArray_ENABLEFLAGS(energies, np.NPY_OWNDATA)

    energies += bqm.offset

    return samples, energies
