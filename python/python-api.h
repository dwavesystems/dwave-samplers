#ifndef ORANG_PYTHON_API_H_INCLUDED
#define ORANG_PYTHON_API_H_INCLUDED

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

#endif
