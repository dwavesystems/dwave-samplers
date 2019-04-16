# distutils: language = c++
# cython: language_level = 3
#
# =============================================================================
# Copyright 2019 D-Wave Systems Inc.
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
#
# =============================================================================
import dimod
import numpy as np
cimport numpy as np

from orang.conversions cimport coo_tables
from orang.orang cimport tables_type, energies_type, samples_type, solveTables

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

samples_dtype = np.intc  # needs to be consistent with samples_type
energies_dtype = np.double  # needs to be consistent with energies_type 

def solve_coo(linear, quadratic, offset, vartype, max_complexity, order=None,
              max_solutions=1):
    """Find ground states of a binary quadratic model in coo format.

    Args:
        linear (array_like):
            A 1D array-like iterable of linear biases.

        quadratic (tuple[array_like, array_like, array_like]):
            A 3-tuple of 1D array_like vectors of the form (row, col, bias).

        offset (numeric, optional):
            Constant offset for the binary quadratic model.

        vartype (:class:`.Vartype`):
            Variable type for the binary quadratic model. Accepted input values:
            * :class:`.Vartype.SPIN`
            * :class:`.Vartype.BINARY`

        max_complexity (int): Upper bound on the exponent of space requirements
            at each step of the algorithm. Should be larger than the treewidth.

        order (array-like, optional):
            The elimination order of the variables. If not provided the
            variables are eliminated in range-order.

        max_solutions (int, optional, default=1):
            The maximum number of solutions to return.

    Returns:
        tuple: A 2-tuple containing the samples and energies as numpy arrays.


    """

    cdef int num_variables = len(linear)
    cdef double low = -1 if vartype is dimod.SPIN else 0
    cdef double beta = -1  # solving

    if not num_variables:
        # nothing to sample from
        raise ValueError("coo must contain at least one variable")

    cdef tables_type tables = coo_tables(linear, quadratic, low, beta)

    cdef int[:] elimination_order
    if order is None:
        elimination_order = np.arange(num_variables, dtype=np.intc)
    else:
        elimination_order = np.asarray(order, dtype=np.intc)

    samples, energies = solve_tables(tables, num_variables, low,
                                     elimination_order, max_complexity,
                                     max_solutions)

    energies += offset

    return samples, energies


cdef solve_tables(tables_type tables, int num_variables, double low,
               int[:] elimination_order, double max_complexity,
               int max_solutions):
    """Find ground states from the given tables."""

    cdef int elimination_order_length = elimination_order.shape[0]
    cdef int* elimination_order_pointer = &elimination_order[0]

    # Pass in a pointer so that solveTables can fill it in. Note that this is
    # a design choice inherited from orang's c++ implementation. In the future
    # we may want to change it.
    cdef int num_energies, srows, scols
    cdef energies_type* energies_pointer
    cdef samples_type* samples_pointer

    solveTables(tables, num_variables,
                 elimination_order_pointer, elimination_order_length,
                 max_complexity, max_solutions,
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

    return samples, energies
