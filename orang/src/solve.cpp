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
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstddef>
#include <cstdlib>

#include <boost/iterator/indirect_iterator.hpp>

#include <base.h>
#include <combine.h>
#include <table.h>
#include <task.h>
#include <treedecomp.h>
#include <buckettree.h>
#include <operations/min.h>

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
typedef std::vector<orang::Table<double>::smartptr> Tables;

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
  Tables tables = isingTables(hLen, hData, jRows, jCols, jData, -1.0);
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
  Tables tables = quboTables(qRows, qCols, qData, -1.0);
  SolveTask task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), 1, minVars);

  solve(task, voData, voLen, maxComplexity, maxSolutions, 0,
    energiesData, energiesLen, solsData, solsRows, solsCols);
}

void solveTables(
  Tables tables, int num_vars,
  int* voData, int voLen,  // elimination order
  double maxComplexity, int maxSolutions,
  double** energiesData, int* energiesLen,
  int** solsData, int* solsRows, int* solsCols) {

  SolveTask task(make_indirect_iterator(tables.begin()),
                 make_indirect_iterator(tables.end()),
                 1,
                 num_vars);

  solve(task, voData, voLen, maxComplexity, maxSolutions, 0,
    energiesData, energiesLen, solsData, solsRows, solsCols);
}
