# distutils: language = c++
# cython: language_level = 3

import numpy as np
cimport numpy as np

from orang.orang cimport tables_type, bias_type, index_type, cooTables

index_dtype = np.uintc  # needs to be consistent with index_type
bias_dtype = np.double  # needs to be consistent with bias_type

cdef tables_type coo_tables(linear, quadratic, double low, double beta):
    """Create the tables for a coo-formatted binary quadratic model.

    Args:
        linear (array_like):
            A 1D array-like iterable of linear biases.

        quadratic (tuple[array_like, array_like, array_like]):
            A 3-tuple of 1D array_like vectors of the form (row, col, bias).

        low (double):
            The low value for the binary values. -1 for SPIN, 0 for BINARY.

        beta (double):
            The inverse temperature.

    Returns:
        tables_type

    Notes:
        To avoid copies, biases should be passed in as NumPy arrays with dtype
        :class:`numpy.double` and the indices should be passed in as NumPy
        arrays with dtype :class:`numpy.intc`.

    """

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
