/**
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
*/
#ifndef ORANG_PYTHON_API_H_INCLUDED
#define ORANG_PYTHON_API_H_INCLUDED

#include <vector>

#include <table.h>

typedef std::vector<orang::Table<double>::smartptr> Tables;

void sample_ising(
  double* hData, int hLen,
  double* jData, int jRows, int jCols,
  int* voData, int voLen,
  double maxComplexity, int numSamples, bool marginals, double beta, int rngSeed,
  double* logPf,
  int** samplesData, int* samplesRows, int* samplesCols,
  double** singleMrgData, int* singleMrgLen,
  double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
  int** pairData, int* pairRows, int* pairCols);

void sample_qubo(
  double* qData, int qRows, int qCols,
  int* voData, int voLen,
  double maxComplexity, int numSamples, bool marginals, double beta, int rngSeed,
  double* logPf,
  int** samplesData, int* samplesRows, int* samplesCols,
  double** singleMrgData, int* singleMrgLen,
  double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
  int** pairData, int* pairRows, int* pairCols);

void sample_tables(
    Tables tables, int numVars,
    int low,
    int* voData, int voLen, double maxComplexity,
    int numSamples,
    bool marginals,
    int rngSeed,
    double* logPf,
    int** samplesData, int* samplesRows, int* samplesCols,
    double** singleMrgData, int* singleMrgLen,
    double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
    int** pairData, int* pairRows, int* pairCols
    );

void solve_ising(
  double* hData, int hLen,
  double* jData, int jRows, int jCols,
  int* voData, int voLen,
  double maxComplexity, int maxSolutions,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols);

void solve_qubo(
  double* qData, int qRows, int qCols,
  int* voData, int voLen,
  double maxComplexity, int maxSolutions,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols);

void solveTables(
  Tables tables, int num_vars,
  int* voData, int voLen, double maxComplexity,
  int maxSolutions,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols);

void sampleTables(
  Tables tables, int numVars,
  int low,
  int* voData, int voLen, double maxComplexity,
  int numSamples,
  bool marginals,
  int rngSeed,
  double* logPf,
  int** samplesData, int* samplesRows, int* samplesCols,
  double** singleMrgData, int* singleMrgLen,
  double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
  int** pairData, int* pairRows, int* pairCols);

#endif
