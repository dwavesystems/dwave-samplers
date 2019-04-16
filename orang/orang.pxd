# distutils: language = c++
# cython: language_level = 3
#
# =============================================================================
# Copyright 2019 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
# =============================================================================
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
