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

#include <vector>
#include <cstddef>

#include "dimod/quadratic_model.h"

using std::size_t;
using std::vector;

using orang::DomIndexVector;
using orang::VarVector;
using orang::Table;

class MallocPtr {
private:
  void* p_;
  MallocPtr(const MallocPtr&);
  MallocPtr& operator=(const MallocPtr&);
public:
  MallocPtr(void* p = 0) : p_(p) {}
  ~MallocPtr() { std::free(p_); }
  void* get() { return p_; }
  void reset(void* p = 0) {
    if (p != p_) {
      std::free(p_);
      p_ = p;
    }
  }
  void* release() {
    void* p = p_;
    p_ = 0;
    return p;
  }
};

inline void* mallocOrThrow(std::size_t sz) {
  void* p = std::malloc(sz);
  if (!p) throw std::bad_alloc();
  return p;
}

VarVector varOrderVec(int voLen, const int* voData, int numVars) {
  if (voLen < 0) throw std::invalid_argument("negative voLen");

  vector<char> seen(numVars, 0);
  VarVector varOrder;
  varOrder.reserve(voLen);

  for (int i = 0; i < voLen; ++i) {
    if (voData[i] < 0 || voData[i] >= numVars) throw std::invalid_argument("variable order out of range");
    if (seen[voData[i]]) throw std::invalid_argument("duplicate variable order entry");
    seen[voData[i]] = 1;
    varOrder.push_back(voData[i]);
  }

  return varOrder;
}

template <class V, class B>
vector<Table<double>::smartptr> getTables(
  dimod::BinaryQuadraticModel<B, V> &bqm,
  double beta,
  int low
) {
    static const DomIndexVector ds1(1, 2);  // [2]  index-linked to vars, specifies the domain size for each variable
    static const DomIndexVector ds2(2, 2);  // [2, 2]

    VarVector vars1(1);
    VarVector vars2(2);
    vector<Table<double>::smartptr> tables;

    auto num_vars = bqm.num_variables();

    for (int i = 0; i < num_vars; ++i) {
        auto var = i;
        if (bqm.linear(var) != 0.0) {
            vars1[0] = var;
            Table<double>::smartptr t(new Table<double>(vars1, ds1));
            (*t)[0] = -beta * bqm.linear(var) * low;
            (*t)[1] = -beta * bqm.linear(var);
            tables.push_back(t);
        }
    }

    for (int i = 0; i < num_vars; ++i) {
        auto var = i;
        auto span = bqm.neighborhood(var);
        while (span.first != span.second) {
            auto neighbor = span.first->v;
            auto bias = span.first->bias;

            if (neighbor > var && bias != 0.0) {
              vars2[0] = var;
              vars2[1] = neighbor;
              Table<double>::smartptr t(new Table<double>(vars2, ds2));

              (*t)[0] = -beta * bias * low * low;
              (*t)[1] = -beta * bias * low;
              (*t)[2] = -beta * bias * low;
              (*t)[3] = -beta * bias;

              tables.push_back(t);
            }
            ++span.first;
        }
    }
    return tables;
}
