# distutils: language = c++
# cython: language_level = 3
import math

from libc.stdint cimport uint32_t, uint16_t
from libcpp cimport bool
from libc.stdlib cimport free
from libcpp.vector cimport vector

from cython.operator cimport dereference as deref

import dimod

import numpy as np
cimport numpy as np


cdef extern from "boost/shared_ptr.hpp" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr() except +
        shared_ptr(T) except +
        shared_ptr(T*) except +
        T operator*()  # boost::detail::sp_dereference< T >::type operator*()

cdef extern from "base.h" namespace "orang":
    ctypedef uint32_t Var
    ctypedef vector[Var] VarVector

    ctypedef uint16_t DomIndex
    ctypedef vector[DomIndex] DomIndexVector

cdef extern from "table.h" namespace "orang":
    cdef cppclass Table[Y]:
        Table() except +
        Table(const VarVector& scope, const DomIndexVector& domSizes) except +
        Y& operator[](size_t)
        size_t size()
        ctypedef shared_ptr[Table[Y]] smartptr

cdef extern from "python-api.h":
    void solve_tables(
        vector[Table[double].smartptr] tables, int num_vars,
        int* voData, int voLen,
        double maxComplexity, int maxSolutions,
        double** energiesData, int* energiesLen,
        int** solsData, int* solsRows, int* solsCols)


cdef Table[double].smartptr make_table(VarVector interaction, double bias, vartype):
    """Return the shared pointer to a table representing an interaction."""

    cdef int num_interactions = len(interaction)

    # specify the size of the domain for each variable (index-linked to interaction)
    cdef DomIndexVector domain = DomIndexVector(num_interactions, 2)  # [2] * len(interaction)

    cdef Table[double].smartptr t = Table[double].smartptr(new Table[double](interaction, domain))

    cdef int table_size = 2 ** num_interactions
    cdef size_t idx

    if vartype is dimod.SPIN:
        for idx in range(table_size):
            deref(t)[idx] = math.pow(-1, num_interactions - bin(idx).count("1")) * bias
    else:
        assert vartype is dimod.BINARY
        for idx in range(table_size - 1):
            deref(t)[idx] = 0.0
        deref(t)[table_size - 1] = bias

    return t


def solve(poly, int num_variables, double complexity,
          vartype,
          int num_reads):
    """
    Args:
        poly: keys should be index-labelled, unique tuples


    """
    if not poly:
        return np.empty((0, 0), dtype=np.intc), np.empty(0, dtype=np.double)

    # set up the tables
    cdef vector[Table[double].smartptr] tables
    for interaction, bias in poly.items():
        tables.push_back(make_table(interaction, bias, vartype))

    cdef np.ndarray[int] variable_order = np.arange(num_variables, dtype=np.intc)

    if vartype is dimod.SPIN:
        low = -1
    else:
        low = 0

    cdef int variable_order_length = variable_order.shape[0]
    cdef int* variable_order_pointer = &variable_order[0]

    # return values
    cdef int num_energies, srows, scols
    cdef double* energies
    cdef int* samples

    solve_tables(tables, num_variables,
                 variable_order_pointer, variable_order_length,
                 complexity, num_reads,
                 &energies, &num_energies,
                 &samples, &srows, &scols
                 )

    np_energies = np.empty(num_energies, dtype=np.double)
    cdef double[::1] energies_view = np_energies
    cdef int i
    for i in range(num_energies):
        energies_view[i] = energies[i]

    free(energies)

    np_samples = np.full((srows, scols), low, dtype=np.intc)
    cdef int [:, ::1] samples_view = np_samples
    cdef int idx = 0
    cdef int row, col
    for row in range(srows):
        for col in range(scols):
            if samples[idx] > 0:
                samples_view[row, col] = 1
            idx = idx + 1

    free(samples)

    return np_samples, np_energies
