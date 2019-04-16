# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np

from orang.orang cimport tables_type, bias_type, index_type, cooTables

index_dtype = np.uintc  # needs to be consistent with index_type
bias_dtype = np.double  # needs to be consistent with bias_type

cdef tables_type coo_tables(linear, quadratic, double low, double beta):

    cdef bias_type[:] ldata = np.asarray(linear, dtype=bias_dtype)
    cdef bias_type* ldata_pointer = &ldata[0]
    cdef size_t num_variables = ldata.shape[0]
    if not num_variables:
        raise ValueError('coo bqm must contain at least one variable')

    cdef index_type[:] irow = np.asarray(quadratic[0], dtype=index_dtype)
    cdef index_type[:] icol = np.asarray(quadratic[1], dtype=index_dtype)
    cdef bias_type[:] qdata = np.asarray(quadratic[2], dtype=bias_dtype)

    cdef size_t num_interactions = qdata.shape[0]
    if icol.shape[0] != num_interactions or irow.shape[0] != num_interactions:
        raise ValueError('inconsistent quadratic vectors')

    cdef bias_type* qdata_pointer = &qdata[0] if num_interactions else NULL
    cdef index_type* irow_pointer = &irow[0] if num_interactions else NULL
    cdef index_type* icol_pointer = &icol[0] if num_interactions else NULL

    return cooTables(num_variables, ldata_pointer,
                     num_interactions, irow_pointer, icol_pointer, qdata_pointer,
                     low,
                     beta)
