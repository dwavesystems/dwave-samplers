# distutils: language = c++
# cython: language_level = 3

import dimod
import numpy as np
cimport numpy as np

from libcpp cimport bool
from libc.stdlib cimport free

from orang.conversions cimport coo_tables
from orang.orang cimport tables_type, samples_type, sampleTables

cdef extern from "numpy/arrayobject.h":
    void PyArray_ENABLEFLAGS(np.ndarray arr, int flags)

samples_dtype = np.intc  # needs to be consistent with samples_type


def sample_coo(linear, quadratic, offset, vartype, beta, max_complexity, order=None,
               num_reads=1, marginals=False, seed=None):

    cdef int num_variables = len(linear)
    cdef double dlow = -1 if vartype is dimod.SPIN else 0
    cdef double b = beta

    cdef int seed_
    if seed is None:
        seed_ = np.random.randint(np.iinfo(np.intc).max, dtype=np.intc)
    else:
        seed_ = seed

    if not num_variables:
        # nothing to sample from
        raise ValueError("coo must contain at least one variable")

    cdef tables_type tables = coo_tables(linear, quadratic, dlow, b)

    cdef int[:] elimination_order
    if order is None:
        elimination_order = np.arange(num_variables, dtype=np.intc)
    else:
        elimination_order = np.asarray(order, dtype=np.intc)

    cdef int ilow = -1 if vartype is dimod.SPIN else 0

    samples, marginals = sample_tables(tables, num_variables, ilow,
                                       elimination_order, max_complexity,
                                       num_reads, marginals, seed_)

    # need to add the offset back in to the log partition function
    marginals['log_partition_function'] += -beta*offset

    return samples, marginals


cdef sample_tables(tables_type tables, int num_variables, int low,
                   int[:] elimination_order, double max_complexity,
                   int num_reads, bool marginals, int seed):

    cdef int elimination_order_length = elimination_order.shape[0]
    cdef int* elimination_order_pointer = &elimination_order[0]

    # pass in pointers so that sampleTables can fill things in
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

    sampleTables(tables, num_variables, low,
                 elimination_order_pointer, elimination_order_length, max_complexity,
                 num_reads,
                 marginals,
                 seed,
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

    return samples, marginal_data
