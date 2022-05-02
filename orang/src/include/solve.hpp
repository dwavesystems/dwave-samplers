// Copyright 2021 D-Wave Systems Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <new>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

#include <cstddef>
#include <cstdlib>

#include <base.h>
#include <table.h>
#include <task.h>
#include <combine.h>
#include <treedecomp.h>
#include <buckettree.h>
#include <operations/min.h>

#include "utils.hpp"

using std::size_t;
using std::numeric_limits;
using std::vector;

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
           int* var_order,
           int num_vars,
           double max_complexity,
           int max_solutions,
           int z,
           double** energies_data, int* energies_len,
           int** sols_data, int* sols_rows, int* sols_cols
) {
  VarVector var_order_vect = varOrderVec(num_vars, var_order, task.numVars());
  TreeDecomp decomp(task.graph(), var_order_vect, task.domSizes());
  if (!(decomp.complexity() <= max_complexity)) throw std::runtime_error("complexity exceeded");

  bool solvable = max_solutions > 0;
  BucketTree<SolveTask> bucket_tree(task, decomp, DomIndexVector(task.numVars()), solvable, false);
  double base_value = bucket_tree.problemValue();

  if (solvable) {
    task.maxSolutions(max_solutions);
    MinSolutionSet<double> solution_set = bucket_tree.solve();
    int num_solutions = static_cast<int>(solution_set.solutions().size());

    // todo: isn't num_vars and task.numVars the same?
    num_vars = static_cast<int>(task.numVars());

    *sols_rows = num_solutions;
    *sols_cols = num_vars;
    if (numeric_limits<size_t>::max() / sizeof(**sols_data) / *sols_rows / *sols_cols == 0) {
      throw std::invalid_argument("solution size too large");
    }
    *sols_data = static_cast<int*>(malloc(*sols_rows * *sols_cols * sizeof(**sols_data)));
    if (!sols_data) throw std::bad_alloc();
    MallocPtr sdmp(*sols_data);

    *energies_len = num_solutions;
    *energies_data = static_cast<double*>(malloc(*energies_len * sizeof(**energies_data)));
    if (!energies_data) throw std::bad_alloc();
    sdmp.release();

    int s[2] = {z, 1};

    MinSolutionSet<double>::solution_set::const_iterator sols_iter = solution_set.solutions().begin();
    for (int i = 0; i < num_solutions; ++i) {
      (*energies_data)[i] = base_value + sols_iter->value;
      for (int j = 0; j < num_vars; ++j) {
        (*sols_data)[i * num_vars + j] = s[sols_iter->solution[j]];
      }
      ++sols_iter;
    }
  } else {
    *sols_rows = 0;
    *sols_cols = 0;
    *sols_data = static_cast<int*>(malloc(1));
    MallocPtr sdmp(sols_data);

    *energies_len = 1;
    *energies_data = static_cast<double*>(malloc(sizeof(double)));
    if (!energies_data) throw std::bad_alloc();
    sdmp.release();
    **energies_data = base_value;
  }
}

} // namespace {anonymous}

template <class V, class B>
void solveBQM(dimod::BinaryQuadraticModel<B, V> &bqm,
              int* var_order,
              double beta,
              int low,
              double max_complexity,
              int max_solutions,
              double** energies_data, int* energies_len,
              int** sols_data, int* sols_rows, int* sols_cols
) {
  vector<Table<double>::smartptr> tables = getTables(bqm, beta, low);

  int num_vars = bqm.num_variables();
  SolveTask task(tables.begin(), tables.end(), 1, num_vars);

  solve(task,
        var_order,
        num_vars,
        max_complexity,
        max_solutions,
        0,
        energies_data,
        energies_len,
        sols_data,
        sols_rows,
        sols_cols);
}
