# distutils: language = c++
# cython: language_level = 3

from orang.orang cimport tables_type

cdef tables_type coo_tables(linear, quadratic, double low, double beta)
