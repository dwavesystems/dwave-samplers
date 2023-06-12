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
#ifndef INCLUDED_ORANG_OPERATIONS_LOGSUMPROD_H
#define INCLUDED_ORANG_OPERATIONS_LOGSUMPROD_H

#include <algorithm>
#include <vector>
#include <cmath>
#include <cstddef>

#include <base.h>
#include <table.h>
#include <combine.h>
#include <marginalizer.h>

namespace orang {

class LogSumMarginalizer : public Marginalizer<double> {
public:
  typedef Table<double> table_type;

private:
  virtual double marginalizeImpl(std::size_t, const table_type& mrgTable) {

    double yMax = *std::max_element(mrgTable.begin(), mrgTable.end());
    double fsum = 0.0;

    for (const auto &y: mrgTable) {
      fsum += exp(y - yMax);
    }

    return yMax + log(fsum);
  }
};



template<typename Rng>
class SolvableLogSumMarginalizer : public SolvableMarginalizer<double, DomIndexVector> {
public:
  typedef Table<double> table_type;

private:
  Rng& rng_;
  varstep_vector inVarsSteps_;
  Var outVar_;
  DomIndex outStepSize_;
  std::vector<double> cumProbs_;

  virtual double marginalizeImpl(std::size_t n, const table_type& mrgTable) {
    using std::transform;
    using std::make_pair;
    typedef std::vector<double>::iterator value_iterator;

    double yMax = *std::max_element(mrgTable.begin(), mrgTable.end());
    double fsum(0.0);

    value_iterator cpBegin = cumProbs_.begin() + n * outStepSize_;
    value_iterator cpIter = cpBegin;
    value_iterator cpEnd = cpBegin + outStepSize_;
    table_type::const_iterator tblIter = mrgTable.begin();
    while (cpIter != cpEnd) {
      fsum += exp(*tblIter - yMax);
      *cpIter = fsum;
      ++cpIter;
      ++tblIter;
    }
    fsum += exp(*tblIter - yMax);
    transform(cpBegin, cpEnd, cpBegin, [&](double a){ return a / fsum; });

    return yMax + log(fsum);
  }

  virtual void solveImpl(DomIndexVector& s) const {
    using std::find_if;
    using std::make_pair;
    typedef std::vector<double>::const_iterator value_const_iterator;

    size_t n = 0;
    for (const auto &vs: inVarsSteps_) {
      n += s[vs.first] * vs.second;
    }
    n *= outStepSize_;

    value_const_iterator cpBegin = cumProbs_.begin() + n;
    value_const_iterator cpEnd = cpBegin + outStepSize_;
    const double rng = rng_();
    value_const_iterator cpFound = find_if(cpBegin, cpEnd, [&](double a){ return rng < a; });
    DomIndex outI = static_cast<DomIndex>(cpFound - cpBegin);

    s[outVar_] = outI;
  }

public:

  SolvableLogSumMarginalizer(Rng& rng,
      const VarVector& inScope, const DomIndexVector& inDomSizes,
      Var outVar, DomIndex outDomSize) :
    rng_(rng),
    inVarsSteps_(),
    outVar_(outVar),
    outStepSize_(outDomSize - 1),
    cumProbs_() {

    size_t numInEntries = buildStepSizes(inScope, inDomSizes, inVarsSteps_);
    cumProbs_.assign(numInEntries * outStepSize_, 0.0);
  }

};

template<typename Rng>
class LogSumProductOperations : public Plus<double> {
public:
  typedef double value_type;
  typedef DomIndexVector solution_type;
  typedef MarginalizerTypes<value_type,solution_type> marginalizer_types;
  typedef typename marginalizer_types::marginalizer_type marginalizer_type;
  typedef typename marginalizer_types::solvablemarginalizer_type solvablemarginalizer_type;
  typedef typename marginalizer_types::marginalizer_smartptr marginalizer_smartptr;
  typedef typename marginalizer_types::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;

  struct CtorArgs {
    Rng& rng;
    CtorArgs(Rng& rng0) : rng(rng0) {}
  };

private:
  Rng& rng_;
  marginalizer_smartptr marginalizer_;

public:
  LogSumProductOperations(const CtorArgs& ca) :
    rng_(ca.rng), marginalizer_(new LogSumMarginalizer) {}

  marginalizer_smartptr marginalizer() const {
    return marginalizer_;
  }

  solvablemarginalizer_smartptr solvableMarginalizer(
      const VarVector& inScope, const DomIndexVector& inDomSizes,
      Var outVar, DomIndex outDomSize) const {
    return solvablemarginalizer_smartptr(
        new SolvableLogSumMarginalizer<Rng>(rng_, inScope, inDomSizes, outVar, outDomSize));
  }

  solution_type initSolution(const DomIndexVector& x0) const {
    return x0;
  }
};

} // namespace orang

#endif
