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

#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>
#include <random>

#include <cstddef>
#include <cstdlib>

#include <base.h>
#include <table.h>
#include <task.h>
#include <treedecomp.h>
#include <buckettree.h>
#include <operations/logsumprod.h>

#include "dimod/quadratic_model.h"
#include "utils.hpp"

using std::size_t;
using std::copy;
using std::map;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::vector;

using orang::Var;
using orang::DomIndexVector;
using orang::VarVector;
using orang::Table;
using orang::TreeDecomp;
using orang::BucketTree;
using orang::TableMerger;

namespace {

class Rng {
  std::mt19937 engine_;
  std::uniform_real_distribution<> distribution_;

public:
  Rng(std::mt19937& engine) :
    engine_(engine) {}

  double operator()() {
    return distribution_(engine_);
  }
};

typedef orang::Task<orang::LogSumProductOperations<Rng> > SampleTask;
typedef std::vector<orang::Table<double>::smartptr> Tables;

typedef pair<Var, Var> VarPair;

struct PairMrgVals {
  double values[4];
};

typedef map<VarPair, PairMrgVals> PairMrgMap;

class Normalizer {
private:
  double log_pf_;

public:
  Normalizer(double log_pf) : log_pf_(log_pf) {}
  double operator()(double x) const { return exp(x - log_pf_); }
};

unsigned int randomSeed(int seed) {
  if (seed >= 0) {
    return seed;
  } else {
    return static_cast<unsigned>(std::rand());
  }
}

vector<double> singleMarginals(const BucketTree<SampleTask>& bucket_tree) {
  vector<double> mrg(bucket_tree.task().numVars());

  VarVector vars1(1);
  TableMerger<SampleTask> merge_tables(bucket_tree.task());
  SampleTask::marginalizer_smartptr marginalizer = bucket_tree.task().marginalizer();
  for (const auto &nt: bucket_tree.nodeTables()) {
    vars1[0] = nt.nodeVar;
    SampleTask::table_smartptr m_table = merge_tables(vars1, nt.tables.begin(),
        nt.tables.end(), *marginalizer);

    Normalizer normalize((*marginalizer)(0, *m_table));
    mrg[nt.nodeVar] = normalize((*m_table)[1]);
  }

  return mrg;
}

PairMrgMap pairMarginals(const BucketTree<SampleTask>& bucket_tree) {
  PairMrgMap mrg;

  for (const auto &t: bucket_tree.task().tables()) {
    if (t->vars().size() == 2) {
      VarPair p(t->vars()[0].index, t->vars()[1].index);
      mrg[p];
    }
  }

  VarVector vars2(2);
  TableMerger<SampleTask> merge_tables(bucket_tree.task());
  SampleTask::marginalizer_smartptr marginalizer = bucket_tree.task().marginalizer();
  for (const auto &nt: bucket_tree.nodeTables()) {
    for (const auto &v : nt.sepVars) {
      VarPair p(min(nt.nodeVar, v), max(nt.nodeVar, v));
      if (mrg.find(p) == mrg.end()) continue;
      vars2[0] = p.first;
      vars2[1] = p.second;
      SampleTask::table_smartptr m_table = merge_tables(vars2, nt.tables.begin(),
          nt.tables.end(), *marginalizer);

      Normalizer normalize((*marginalizer)(0, *m_table));
      PairMrgVals& mv = mrg[p];
      mv.values[0] = normalize((*m_table)[0]);
      mv.values[1] = normalize((*m_table)[1]);
      mv.values[2] = normalize((*m_table)[2]);
      mv.values[3] = normalize((*m_table)[3]);
    }
  }
  return mrg;
}

void sample(SampleTask& task,
            int* var_order,
            int z,
            int num_vars,
            double max_complexity,
            int num_samples,
            bool marginals,
            double* log_pf,
            int** samples_data, int* samples_rows, int* samples_cols,
            double** single_mrg_data, int* single_mrg_len,
            double** pair_mrg_data, int* pair_mrg_rows, int* pair_mrg_cols,
            int** pair_data, int* pair_rows, int* pair_cols
) {
  // todo: isn't num_vars and task.numVars the same?
  VarVector var_order_vect = varOrderVec(num_vars, var_order, task.numVars());

  TreeDecomp decomp(task.graph(), var_order_vect, task.domSizes());
  if (!(decomp.complexity() <= max_complexity)) throw std::runtime_error("complexity exceeded");

  bool solvable = num_samples > 0;

  BucketTree<SampleTask> bucket_tree(task, decomp, DomIndexVector(task.numVars()), solvable, marginals);
  *log_pf = bucket_tree.problemValue();

  MallocPtr samples_mp;
  MallocPtr single_mrg_mp;
  MallocPtr pair_mrg_mp;
  MallocPtr pair_mp;

  if (solvable) {
    num_vars = static_cast<int>(task.numVars());
    *samples_rows = num_samples;
    *samples_cols = num_vars;
    if (numeric_limits<size_t>::max() / sizeof(**samples_data) / *samples_rows / *samples_cols == 0) {
      throw std::invalid_argument("samples size too large");
    }
    samples_mp.reset(mallocOrThrow(static_cast<size_t>(*samples_rows) * *samples_cols * sizeof(**samples_data)));

    int s[2] = {z, 1};
    int* samples_mp_data = static_cast<int*>(samples_mp.get());
    for (int i = 0; i < num_samples; ++i) {
      DomIndexVector samp = bucket_tree.solve();
      for (int j = 0; j < num_vars; ++j) {
        *samples_mp_data++ = s[samp[j]];
      }
    }

  } else {
    *samples_rows = 0;
    *samples_cols = 0;
    samples_mp.reset(mallocOrThrow(1));
  }

  if (marginals) {
    vector<double> single_mrg = singleMarginals(bucket_tree);
    *single_mrg_len = static_cast<int>(single_mrg.size());
    single_mrg_mp.reset(mallocOrThrow(single_mrg.size() * sizeof(**single_mrg_data)));
    copy(single_mrg.begin(), single_mrg.end(), static_cast<double*>(single_mrg_mp.get()));

    PairMrgMap pair_mrg = pairMarginals(bucket_tree);
    *pair_mrg_rows = static_cast<int>(pair_mrg.size());
    *pair_mrg_cols = 4;
    pair_mrg_mp.reset(mallocOrThrow(pair_mrg.size() * 4 * sizeof(**pair_mrg_data)));
    *pair_rows = static_cast<int>(pair_mrg.size());
    *pair_cols = 2;
    pair_mp.reset(mallocOrThrow(pair_mrg.size() * 2 * sizeof(**pair_data)));
    int* pair_mp_data = static_cast<int*>(pair_mp.get());
    double* pair_mrg_mp_data = static_cast<double*>(pair_mrg_mp.get());
    for (const auto &e: pair_mrg) {
      *pair_mp_data++ = static_cast<int>(e.first.first);
      *pair_mp_data++ = static_cast<int>(e.first.second);
      *pair_mrg_mp_data++ = e.second.values[0];
      *pair_mrg_mp_data++ = e.second.values[1];
      *pair_mrg_mp_data++ = e.second.values[2];
      *pair_mrg_mp_data++ = e.second.values[3];
    }

  } else {
    *single_mrg_len = 0;
    single_mrg_mp.reset(mallocOrThrow(1));

    *pair_mrg_rows = 0;
    *pair_mrg_cols = 4;
    pair_mrg_mp.reset(mallocOrThrow(1));

    *pair_rows = 0;
    *pair_cols = 2;
    pair_mp.reset(mallocOrThrow(1));
  }

  *samples_data = static_cast<int*>(samples_mp.release());
  *single_mrg_data = static_cast<double*>(single_mrg_mp.release());
  *pair_mrg_data = static_cast<double*>(pair_mrg_mp.release());
  *pair_data = static_cast<int*>(pair_mp.release());
}

} // namespace {anonymous}

template <class V, class B>
void sampleBQM(
  dimod::BinaryQuadraticModel<B, V> &bqm,
  int* var_order,
  double beta,
  int low,
  double max_complexity,
  int num_samples,
  bool marginals,
  int seed,
  double* log_pf,
  int** samples_data, int* samples_rows, int* samples_cols,
  double** single_mrg_data, int* single_mrg_len,
  double** pair_mrg_data, int* pair_mrg_rows, int* pair_mrg_cols,
  int** pair_data, int* pair_rows, int* pair_cols
) {
    std::mt19937 engine(randomSeed(seed));
    Rng rng(engine);

    vector<Table<double>::smartptr> tables = getTables(bqm, beta, low);

    int num_vars = bqm.num_variables();
    SampleTask task(tables.begin(), tables.end(), rng, num_vars);

    sample(task,
           var_order,
           low,   // -1 for SPIN, 0 for BINARY
           num_vars,
           max_complexity,
           num_samples,
           marginals,
           log_pf,
           samples_data, samples_rows, samples_cols,
           single_mrg_data, single_mrg_len,
           pair_mrg_data, pair_mrg_rows, pair_mrg_cols,
           pair_data, pair_rows, pair_cols);
}
