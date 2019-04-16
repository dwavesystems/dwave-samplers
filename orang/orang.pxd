# distutils: language = c++
# cython: language_level = 3

cimport numpy as np

from libcpp cimport bool
from libcpp.vector cimport vector

ctypedef double bias_type
ctypedef double energies_type
ctypedef unsigned int index_type
ctypedef int samples_type


cdef extern from "boost/shared_ptr.hpp" namespace "boost":
    cdef cppclass shared_ptr[T]:
        pass

cdef extern from "table.h" namespace "orang":
    cdef cppclass Table[Y]:
        ctypedef shared_ptr[Table[Y]] smartptr

ctypedef vector[Table[bias_type].smartptr] tables_type

cdef extern from "conversions.h":
    tables_type cooTables(
        size_t numLinear,
        const bias_type* lVals,
        size_t numQuadratic,
        const index_type* iRow, const index_type* iCol, const bias_type* qVals,
        double low,
        double beta) except +

cdef extern from "python-api.h":
    void solveTables(
        tables_type tables, int num_vars,
        int* voData, int voLen, double maxComplexity,
        int maxSolutions,
        energies_type** energiesData, int* energiesLen,
        samples_type** solsData, int* solsRows, int* solsCols) except +

    void sampleTables(
        tables_type tables, int num_vars, int low,
        int* voData, int voLen, double maxComplexity,
        int numSamples,
        bool marginals,
        int rngSeed,
        double* logPf,
        samples_type** samplesData, int* samplesRows, int* samplesCols,
        double** singleMrgData, int* singleMrgLen,
        double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
        int** pairData, int* pairRows, int* pairCols) except +
