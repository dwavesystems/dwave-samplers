# distutils: language = c++

cimport cython

from libc.stdlib cimport free
from libcpp cimport bool

import numpy as np
cimport numpy as np


cdef extern from "python-api.h":
    void solve_qubo(
        double* qData, int qRows, int qCols,
        int* voData, int voLen,
        double maxComplexity, int maxSolutions,
        double** energiesData, int* energiesLen,
        int** solsData, int* solsRows, int* solsCols)

    void sample_qubo(
        double* qData, int qRows, int qCols,
        int* voData, int voLen,
        double maxComplexity,
        int numSamples,
        bool marginals,
        double beta,
        int rngSeed,
        double* logPf,
        int** samplesData, int* samplesRows, int* samplesCols,
        double** singleMrgData, int* singleMrgLen,
        double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
        int** pairData, int* pairRows, int* pairCols);


def sample(double[:] qdata, int num_variables, int low,
           int num_reads,
           int[:] variable_order, double complexity,
           bool marginals,
           double beta,
           int seed):
    """

    Args:

        offset

        qdata (double[:]): The qubo matrix flattened.

        num_variables (int): The size of the dimension of the unflattend qubo matrix.

        variable_order (int[:]): The order in which the variables should be eliminated.

        complexity (double): The maximum complexity, should be greater than or equal to
            the treewidth + 1.

        marginals: whether or not to compute marginals (boolean).
            Marginals are always returned but will be empty if this
            parameter is False.

        beta: Boltzmann distribution inverse temperature parameter

        seed: random number generator seed.  Negative values will cause
            a time-based seed to be used.

    Returns:
        (log_pf, samples, single_mrg, pair_mrg)
        log_pf: natural log of the partition function
        samples: numpy array of samples (each row is a sample)
        single_mrg: 1D numpy array of single-variable marginals.  Each
            value is the probability that the corresponding variable is 1.
        pair_mrg: numpy array of pairwise marginals.  Corresponding
            indices are in pairs return value.  Each row contains four
            values: [m_i0_j0, m_i1_j0, m_i0_j1, m_i1_j1] where m_ix_jy
            is the probability that variables i and j have respective
            values x and y, where [i, j] is the corresponding row of the
            pairs return value.
        pairs: numpy array of pairwise marginal indices.  These are
            precisely the pairs with nonzero Q coefficient."

    """
    if num_variables == 0:
        raise NotImplementedError
        return np.empty((0, 0), dtype=np.intc), np.empty(0, np.double)

    cdef double* qdata_pointer = &qdata[0]

    cdef int variable_order_length = variable_order.shape[0]
    cdef int* variable_order_pointer = &variable_order[0]

    # return values

    cdef double logPf

    cdef int* samples
    cdef int srows, scols

    cdef double* single_marginals_data
    cdef int smlen

    cdef double* pair_marginals_data
    cdef int pmrows, pmcols

    cdef int* pair_data
    cdef int prows, pcols

    sample_qubo(qdata_pointer, num_variables, num_variables,
                variable_order_pointer, variable_order_length,
                complexity,
                num_reads,
                marginals,
                beta,
                seed,
                &logPf,
                &samples, &srows, &scols,
                &single_marginals_data, &smlen,
                &pair_marginals_data, &pmrows, &pmcols,
                &pair_data, &prows, &pcols
                )


    np_samples = np.full((srows, scols), low, dtype=np.intc)
    cdef int idx = 0
    cdef int row, col
    for row in range(srows):
        for col in range(scols):
            if samples[idx] > 0:
                np_samples[row, col] = 1
            idx = idx + 1

    free(samples)

    cdef int i
    variable_marginals = np.empty(smlen, dtype=np.double)
    for i in range(smlen):
        variable_marginals[i] = single_marginals_data[i]

    free(single_marginals_data)

    interaction_marginals = np.empty((pmrows, pmcols), dtype=np.double)
    idx = 0
    for row in range(pmrows):
        for col in range(pmcols):
            interaction_marginals[row, col] = pair_marginals_data[idx]
            idx = idx + 1

    free(pair_marginals_data)

    pairs = np.empty((prows, pcols), dtype=np.intc)
    idx = 0
    for row in range(prows):
        for col in range(pcols):
            pairs[row, col] = pair_data[idx]
            idx = idx + 1

    free(pair_data)

    return np_samples, float(logPf), variable_marginals, interaction_marginals, pairs

@cython.wraparound(False)
def solve(double offset, double[:] qdata, int num_variables, int low,
          int num_reads,
          int[:] variable_order, double complexity):
    """

    Args:
        offset (double): The constant energy offset

        qdata (double[:]): The qubo matrix flattened.

        num_variables (int): The size of the dimension of the unflattend qubo matrix.

        low (int): -1 for SPIN, 0 for BINARY

        variable_order (int[:]): The order in which the variables should be eliminated.

        complexity (double): The maximum complexity, should be greater than or equal to
            the treewidth + 1.

    Returns:
        samples, energies as numpy arrays.

    Notes:
        not a lot of checking is done

    """

    if num_variables == 0:
        return np.empty((0, 0), dtype=np.intc), np.empty(0, np.double)

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

    np_energies = np.empty(num_energies, dtype=np.double)
    cdef int i
    for i in range(num_energies):
        np_energies[i] = energies[i] + offset

    free(energies)

    np_samples = np.full((srows, scols), low, dtype=np.intc)
    cdef int idx = 0
    cdef int row, col
    for row in range(srows):
        for col in range(scols):
            if samples[idx] > 0:
                np_samples[row, col] = 1
            idx = idx + 1

    free(samples)

    return np_samples, np_energies
