# distutils: language = c++
# cython: language_level = 3

import dimod

import numpy as np
cimport numpy as np

from libc.stdint cimport uint32_t, uint16_t
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


cdef Table[double].smartptr make_table(VarVector interaction, double bias, int low):
    """Return the shared pointer to a table representing an interaction."""

    # dev note: should this be size_t?
    cdef int num_interactions = len(interaction)

    # specify the size of the domain for each variable (index-linked to
    # the interactions). Because we're over binary variables
    # domain = [2] * len(interactions)
    cdef DomIndexVector domain = DomIndexVector(num_interactions, 2)

    cdef Table[double].smartptr t = Table[double].smartptr(new Table[double](interaction, domain))

    cdef size_t table_size = 2 ** num_interactions
    cdef unsigned int idx

    # specify the energy for all possible combinations of the interaction
    # variables. We take the configuration of the variables and convert it
    # to an integer to get the index. We treat -1 as 0.
    # So (-1, 1, -1) => idx == 2
    for idx in range(table_size):
        deref(t)[idx] = bias * low ** (num_interactions - __builtin_popcount(idx))

    return t


cdef class cyTables:
    cdef vector[Table[double].smartptr] tables
    cdef int low  # 0 for binary, -1 for spin
    cdef int num_variables

    def __init__(self, int num_variables, vartype):
        if vartype is dimod.SPIN:
            self.low = -1
        else:
            self.low = 0

        self.num_variables = num_variables

    def _add_variable(self, int variable, double bias):
        self.tables.push_back(make_table([variable], bias, self.low))

    def _add_interaction(self, interaction, double bias):
        self.tables.push_back(make_table(interaction, bias, self.low))


class Tables(cyTables, object):

    @dimod.decorators.vartype_argument('vartype')
    def __init__(self, num_variables, vartype):
        super(Tables, self).__init__(num_variables, vartype)

        self.vartype = vartype

    # @classmethod
    # def from_ising(cls, h, J):
    #     bqm = dimod.BQM.from_ising(h, J)

    #     variables = dimod.variables.Variables(bqm.variables)

    #     ldata, (irow, icol, qdata), off = bqm.to_numpy_vectors(variable_order=variables)

    #     return cls.from_coo(ldata, irow, icol, qdata, off, dimod.SPIN, variables)

    @classmethod
    def from_coo(cls, double[:] ldata, long[:] irow, long[:] icol, double[:] qdata,
                 double offset, vartype):

        cdef int idx

        cdef int num_variables = ldata.shape[0]
        cdef int num_interactions = qdata.shape[0]

        tables = cls(num_variables, vartype)

        for idx in range(num_variables):
            tables._add_interaction([idx], ldata[idx])

        for idx in range(num_interactions):
            interaction = [irow[idx], icol[idx]]
            tables._add_interaction(interaction, qdata[idx])

        return tables

def solve(cyTables tables, int[:] elimination_order, double complexity, int num_reads):

    cdef int elimination_order_length = elimination_order.shape[0]
    cdef int* elimination_order_pointer = &elimination_order[0]
 
    # Need to create values that solve_tables can fill in. Note that this is
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

    # create a numpy array without making a copy. Numpy manages memory and
    # cleanup
    energies = np.asarray(<double[:num_energies]> energies_pointer)
    samples = np.asarray(<int[:srows, :scols]> samples_pointer)

    return samples, energies
