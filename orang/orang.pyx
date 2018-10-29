# distutils: language = c++
from libc.stdlib cimport free

import numpy as np
cimport numpy as np

import dimod
import dwave_networkx as dnx


cdef extern from "python-api.h":
    void solve_qubo(
        double* qData, int qRows, int qCols,
        int* voData, int voLen,
        double maxComplexity, int maxSolutions,
        double** energiesData, int* energiesLen,
        int** solsData, int* solsRows, int* solsCols)


def _solve(double[:] qdata, int num_variables,
            int num_reads,
            int[:] variable_order,
            double complexity):

    cdef double* qdata_pointer = &qdata[0]

    cdef int variable_order_length = variable_order.shape[0]
    cdef int* variable_order_pointer = &variable_order[0]

    # return values
    cdef int num_energies, srows, scols
    cdef double* energies
    cdef int* samples

    solve_qubo(qdata_pointer, num_variables, num_variables,
               variable_order_pointer, variable_order_length,
               complexity, num_reads,
               &energies, &num_energies,
               &samples, &srows, &scols
               )

    cdef int i
    np_energies = np.fromiter((energies[i] for i in range(num_energies)), dtype=np.double)

    np_samples = np.fromiter((samples[i] for i in range(srows*scols)), dtype=np.intc)

    free(energies)
    free(samples)

    return np_samples.reshape((srows, scols)), np_energies



def sample(bqm, num_reads, order=None):

    num_variables = len(bqm)

    if order is None:
        tw, order = dnx.min_fill_heuristic(bqm.adj) 
        

    variables = list(bqm.variables)  # need an ordered list


    Q = bqm.binary.to_numpy_matrix(variable_order=variables)
    qdata = np.asarray(Q.flatten(order='C'), dtype=np.double)  # make sure C ordered an correct dtype

    variable_order = np.arange(num_variables, dtype=np.intc)  # to match the given order

    samples, energies = _solve(qdata, num_variables, num_reads, variable_order, tw+1.)

    if bqm.vartype is dimod.SPIN:
        samples = 2 * samples - 1

    return dimod.SampleSet.from_samples((samples, variables), bqm.vartype, energy=energies+bqm.binary.offset)
