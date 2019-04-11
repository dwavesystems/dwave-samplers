# distutils: language = c++
# cython: language_level = 3

import dimod

import numpy as np
cimport numpy as np

from libc.stdint cimport uint32_t, uint16_t
from libc.stdlib cimport free
from libcpp.vector cimport vector

from cython cimport view
from cython.operator cimport dereference as deref

cdef extern int __builtin_popcount(unsigned int) nogil
# cdef extern int __builtin_popcountll(unsigned long long) nogil

cdef extern from "boost/shared_ptr.hpp" namespace "boost":
    cdef cppclass shared_ptr[T]:
        shared_ptr() except +
        shared_ptr(T) except +
        shared_ptr(T*) except +
        T operator*()  except + # boost::detail::sp_dereference< T >::type operator*()

cdef extern from "base.h" namespace "orang":
    ctypedef uint32_t Var
    ctypedef vector[Var] VarVector

    ctypedef uint16_t DomIndex
    ctypedef vector[DomIndex] DomIndexVector

cdef extern from "numpy/arrayobject.h":
    # we need this to give numpy ownership of arrays generated in the C code
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef extern from "table.h" namespace "orang":
    cdef cppclass Table[Y]:
        Table() except +
        Table(const VarVector& scope, const DomIndexVector& domSizes) except +
        Y& operator[](size_t)
        size_t size()
        ctypedef shared_ptr[Table[Y]] smartptr
        void set_value(size_t, Y)


cdef extern from "python-api.h":
    void solve_tables(
        vector[Table[double].smartptr] tables, int num_vars,
        int* voData, int voLen,
        double maxComplexity, int maxSolutions,
        double** energiesData, int* energiesLen,
        int** solsData, int* solsRows, int* solsCols)


cdef Table[double].smartptr make_table(VarVector interaction, double bias, int low):
    """Return the shared pointer to a table representing an interaction.

    Note: interaction must be sorted
    """

    # dev note: should this be size_t?
    cdef int num_interactions = len(interaction)

    # specify the size of the domain for each variable (index-linked to
    # the interactions). Because we're over binary variables
    # domain = [2] * len(interactions)
    cdef DomIndexVector domain = DomIndexVector(num_interactions, 2)

    cdef Table[double].smartptr t = Table[double].smartptr(new Table[double](interaction, domain))

    cdef size_t table_size = 2 ** num_interactions
    cdef size_t idx
    cdef double energy

    table = deref(t)

    # specify the energy for all possible combinations of the interaction
    # variables. We take the configuration of the variables and convert it
    # to an integer to get the index. We treat -1 as 0.
    # So (-1, 1, -1) => idx == 2
    for idx in range(table_size):
        energy = bias * low ** (num_interactions - __builtin_popcount(idx))
        print('in', idx, energy)
        deref(t)[idx] = energy
        print('right after', deref(t)[idx])
        table[idx] = energy
        print('after after', deref(t)[idx], table[idx])

    print(num_interactions)
    print(table_size)
    for idx in range(table_size):
        print('out', idx, deref(t)[idx])

    return t


cdef class cyTables:
    cdef vector[Table[double].smartptr] tables
    cdef int low  # 0 for binary, -1 for spin
    cdef int num_variables
    cdef double offset

    def __init__(self, int num_variables, vartype):
        if vartype is dimod.SPIN:
            self.low = -1
        else:
            self.low = 0

        self.num_variables = num_variables

        self.offset = 0

    def add_interaction(self, VarVector interaction, double bias):
        # interaction must be sorted
        self.tables.push_back(make_table(interaction, bias, self.low))

    def add_offset(self, double offset):
        self.offset += offset

    def __len__(self):
        # number of interactions in the table
        return self.tables.size()


class Tables(cyTables, object):

    @dimod.decorators.vartype_argument('vartype')
    def __init__(self, num_variables, vartype):
        super(Tables, self).__init__(num_variables, vartype)

        self.vartype = vartype

    @classmethod
    def from_coo(cls, linear, quadratic, offset, vartype):

        # We could do this stuff as an argument type but we expect this object to
        # be used from python so doing it this way gets us much more flexibility
        cdef double[:] ldata = np.asarray(linear, np.double)
        cdef double[:] qdata = np.asarray(quadratic[2], np.double)
        cdef unsigned int[:] irow = np.asarray(quadratic[0], np.uint32)
        cdef unsigned int[:] icol = np.asarray(quadratic[1], np.uint32)

        cdef int num_variables = ldata.shape[0]
        cdef int num_interactions = qdata.shape[0]

        tables = cls(num_variables, vartype)

        cdef unsigned int idx, u, v
        for idx in range(num_variables):
            tables.add_interaction([idx], ldata[idx])

        for idx in range(num_interactions):
            u = irow[idx]
            v = icol[idx]
            if u < v:
                tables.add_interaction([u, v], qdata[idx])
            else:
                tables.add_interaction([v, u], qdata[idx])

        tables.add_offset(offset)

        return tables

def solve(cyTables tables, elimination_order, double complexity, int num_reads):


    cdef int[:] order = np.asarray(elimination_order, np.intc) 
    cdef int elimination_order_length = order.shape[0]
    cdef int* elimination_order_pointer = &order[0]
 
    # Pass in a pointer so that solve_tables can fill it in. Note that this is
    # a design choice inherited from orang's c++ implementation. In the future
    # we may want to change it.
    cdef int num_energies, srows, scols
    cdef double* energies_pointer
    cdef int* samples_pointer

    solve_tables(tables.tables, tables.num_variables,
                 elimination_order_pointer, elimination_order_length,
                 complexity, num_reads,
                 &energies_pointer, &num_energies,
                 &samples_pointer, &srows, &scols
                 )

    # create a numpy array without making a copy
    energies = np.asarray(<double[:num_energies]> energies_pointer)
    samples = np.asarray(<int[:srows, :scols]> samples_pointer)

    # tell numpy that it needs to free the memory when the array is garbage
    # collected
    PyArray_ENABLEFLAGS(energies, np.NPY_OWNDATA)
    PyArray_ENABLEFLAGS(samples, np.NPY_OWNDATA)

    return samples, energies
