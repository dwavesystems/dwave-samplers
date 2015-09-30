#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstddef>
#include <cstdlib>

#include <boost/iterator/indirect_iterator.hpp>

#include <orang/base.h>
#include <orang/combine.h>
#include <orang/table.h>
#include <orang/task.h>
#include <orang/treedecomp.h>
#include <orang/buckettree.h>
#include <orang/operations/min.h>

#include "python-api.h"
#include "conversions.h"

using std::size_t;
using std::max;
using std::numeric_limits;
using std::vector;

using boost::make_indirect_iterator;

using orang::DomIndexVector;
using orang::VarVector;
using orang::Table;
using orang::TreeDecomp;
using orang::BucketTree;
using orang::MinSolution;
using orang::MinSolutionSet;

typedef orang::Task<orang::MinOperations<double, orang::Plus<double> > > SolveTask;

namespace {


void solve(SolveTask& task,
  int* voData, int voLen,
  double maxComplexity, int maxSolutions, int z,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols) {

  VarVector varOrder = varOrderVec(voLen, voData, task.numVars());
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  if (!(decomp.complexity() <= maxComplexity)) throw std::runtime_error("complexity exceeded");

  bool solvable = maxSolutions > 0;
  BucketTree<SolveTask> bucketTree(task, decomp, DomIndexVector(task.numVars()), solvable, false);
  double baseValue = bucketTree.problemValue();

  if (solvable) {
    task.maxSolutions(maxSolutions);
    MinSolutionSet<double> solutionSet = bucketTree.solve();
    int numSolutions = static_cast<int>(solutionSet.solutions().size());
    int numVars = static_cast<int>(task.numVars());

    *solsRows = numSolutions;
    *solsCols = numVars;
    if (numeric_limits<size_t>::max() / sizeof(**solsData) / *solsRows / *solsCols == 0) {
      throw std::invalid_argument("solution size too large");
    }
    *solsData = static_cast<int*>(malloc(*solsRows * *solsCols * sizeof(**solsData)));
    if (!solsData) throw std::bad_alloc();
    MallocPtr sdmp(*solsData);

    *energiesLen = numSolutions;
    *energiesData = static_cast<double*>(malloc(*energiesLen * sizeof(**energiesData)));
    if (!energiesData) throw std::bad_alloc();
    sdmp.release();

    int s[2] = {z, 1};

    MinSolutionSet<double>::solution_set::const_iterator solsIter = solutionSet.solutions().begin();
    for (int i = 0; i < numSolutions; ++i) {
      (*energiesData)[i] = baseValue + solsIter->value;
      for (int j = 0; j < numVars; ++j) {
        (*solsData)[i * numVars + j] = s[solsIter->solution[j]];
      }
      ++solsIter;
    }

  } else {
    *solsRows = 0;
    *solsCols = 0;
    *solsData = static_cast<int*>(malloc(1));
    MallocPtr sdmp(solsData);

    *energiesLen = 1;
    *energiesData = static_cast<double*>(malloc(sizeof(double)));
    if (!energiesData) throw std::bad_alloc();
    sdmp.release();
    **energiesData = baseValue;

  }
}

} // namespace {anonymous}

void solve_ising(
  double* hData, int hLen,
  double* jData, int jRows, int jCols,
  int* voData, int voLen,
  double maxComplexity, int maxSolutions,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols) {

  int minVars = max(hLen, max(jRows, jCols));
  vector<Table<double>::smartptr> tables = isingTables(hLen, hData, jRows, jCols, jData, -1.0);
  SolveTask task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), 1, minVars);

  solve(task, voData, voLen, maxComplexity, maxSolutions, -1,
    energiesData, energiesLen, solsData, solsRows, solsCols);
}

void solve_qubo(
  double* qData, int qRows, int qCols,
  int* voData, int voLen,
  double maxComplexity, int maxSolutions,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols) {

  int minVars = max(qRows, qCols);
  vector<Table<double>::smartptr> tables = quboTables(qRows, qCols, qData, -1.0);
  SolveTask task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), 1, minVars);

  solve(task, voData, voLen, maxComplexity, maxSolutions, 0,
    energiesData, energiesLen, solsData, solsRows, solsCols);
}
