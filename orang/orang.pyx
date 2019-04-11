# distutils: language = c++
# cython: language_level = 3

import dimod
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector


cdef extern from "boost/shared_ptr.hpp" namespace "boost":
    cdef cppclass shared_ptr[T]:
        pass
        # shared_ptr() except +
        # shared_ptr(T) except +
        # shared_ptr(T*) except +
        # T operator*()  except + # boost::detail::sp_dereference< T >::type operator*()

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

cdef extern from "table.h" namespace "orang":
    cdef cppclass Table[Y]:
        pass
        # Table() except +
        # Table(const VarVector& scope, const DomIndexVector& domSizes) except +
        # Y& operator[](size_t)
        # size_t size()
        ctypedef shared_ptr[Table[Y]] smartptr
        # void set_value(size_t, Y)

cdef extern from "conversions.h":
    vector[Table[double].smartptr] cooTables(
        size_t numLinear,
        const double* lVals,
        size_t numQuadratic,
        const unsigned int* iRow, const unsigned int* iCol, const double* qVals,
        double low,
        double beta) except +

cdef extern from "python-api.h":
    void solve_tables(
        vector[Table[double].smartptr] tables, int num_vars,
        int* voData, int voLen,
        double maxComplexity, int maxSolutions,
        double** energiesData, int* energiesLen,
        int** solsData, int* solsRows, int* solsCols) except +


def solve_coo(linear, quadratic, offset, vartype, complexity, order=None, max_solutions=1):

    ldata = np.asarray(linear, dtype=np.double)
    irow = np.asarray(quadratic[0], dtype=np.uintc)
    icol = np.asarray(quadratic[1], dtype=np.uintc)
    qdata = np.asarray(quadratic[2], dtype=np.double)

    num_variables = ldata.shape[0]
    if not num_variables:
        raise ValueError('coo bqm must contain at least one variable')

    num_interactions = qdata.shape[0]
    if icol.shape[0] != num_interactions or irow.shape[0] != num_interactions:
        raise ValueError('inconsistent quadratic vectors')

    if order is None:
        order = np.arange(len(linear), dtype=np.intc)
    else:
        order = np.asarray(order, dtype=np.intc)

    samples, energies =  _solve_coo(
        ldata,
        irow, icol, qdata,
        -1 if vartype is dimod.SPIN else 0,
        order,
        complexity,
        max_solutions)

    energies += offset

    return samples, energies


cdef _solve_coo(double[:] ldata,
               unsigned int[:] irow,
               unsigned int[:] icol,
               double[:] qdata,
               double low,
               int[:] elimination_order, double complexity,
               int max_solutions):

    cdef size_t num_variables = ldata.shape[0]
    cdef size_t num_interactions = qdata.shape[0]

    cdef double beta = -1  # solving

    cdef double* ldata_pointer = &ldata[0]

    cdef double* qdata_pointer = &qdata[0] if num_interactions else NULL
    cdef unsigned int* irow_pointer = &irow[0] if num_interactions else NULL
    cdef unsigned int* icol_pointer = &icol[0] if num_interactions else NULL

    # convert the problem to a form that orang can understand
    cdef vector[Table[double].smartptr] tables = cooTables(
                           num_variables,
                           ldata_pointer,
                           num_interactions,
                           irow_pointer, icol_pointer, qdata_pointer,
                           low, beta)

    return _solve(tables, num_variables, low,
                  elimination_order, complexity,
                  max_solutions)


cdef _solve(vector[Table[double].smartptr] tables, int num_variables, double low,
               int[:] elimination_order, double complexity,
               int max_solutions):

    cdef int elimination_order_length = elimination_order.shape[0]
    cdef int* elimination_order_pointer = &elimination_order[0]

    # Pass in a pointer so that solve_tables can fill it in. Note that this is
    # a design choice inherited from orang's c++ implementation. In the future
    # we may want to change it.
    cdef int num_energies, srows, scols
    cdef double* energies_pointer
    cdef int* samples_pointer

    solve_tables(tables, num_variables,
                 elimination_order_pointer, elimination_order_length,
                 complexity, max_solutions,
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

    # convert the samples to spin if necessary
    cdef size_t i, j
    if low == -1:
        for i in range(srows):
            for j in range(scols):
                if samples[i, j] == 0:
                    samples[i, j] = -1

    return samples, energies
