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
#ifndef INCLUDED_ORANG_OPERATIONS_MIN_H
#define INCLUDED_ORANG_OPERATIONS_MIN_H

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <set>
#include <vector>
#include <utility>
#include <functional>

#include <base.h>
#include <exception.h>
#include <table.h>
#include <marginalizer.h>

namespace orang {

template<typename Y, typename Compare=std::less<Y> >
class MinMarginalizer : public Marginalizer<Y> {
public:
  typedef Y value_type;
  typedef Table<value_type> table_type;
  typedef Compare value_compare;

private:
  virtual value_type marginalizeImpl(std::size_t, const table_type& mrgTable) {
    return *std::min_element(mrgTable.begin(), mrgTable.end(), value_compare());
  }
};


template<typename Y>
struct MinSolution {
  typedef Y value_type;
  value_type value;
  DomIndexVector solution;

  MinSolution(const value_type& value, const DomIndexVector& solution) : value(value), solution(solution) {}
};

template<typename Y, typename Compare=std::less<Y> >
struct MinSolutionCompare {
  bool operator()(const MinSolution<Y>& s1, const MinSolution<Y>& s2) const {
    return Compare()(s1.value, s2.value) || (s1.value == s2.value && s1.solution < s2.solution);
  }
};


template<typename Y, typename Compare=std::less<Y> >
class MinSolutionSet {
public:
  typedef Y value_type;
  typedef Compare value_compare;
  typedef MinSolution<value_type> solution_type;
  typedef std::set<solution_type, MinSolutionCompare<value_type, value_compare> > solution_set;

private:
  const std::size_t maxSolutions_;
  solution_set solSet_;

public:
  explicit MinSolutionSet(std::size_t maxSolutions) : maxSolutions_(maxSolutions), solSet_() {
    if (maxSolutions_ == 0) {
      throw InvalidArgumentException("maxSolutions must be positive");
    }
  }

  std::size_t maxSolutions() const { return maxSolutions_; }
  solution_set& solutions() { return solSet_; }
  const solution_set& solutions() const { return solSet_; }
};

template<typename Y, typename Cmb, typename Compare> class MinOperations;

template<typename Y, typename Cmb, typename Compare=std::less<Y> >
class SolvableMinMarginalizer : public SolvableMarginalizer<Y, MinSolutionSet<Y, Compare> > {
public:
  typedef Y value_type;
  typedef Compare value_compare;
  typedef Table<value_type> table_type;
  typedef MinSolutionSet<value_type, value_compare> solution_type;
  typedef MinSolutionCompare<value_type, value_compare> solution_compare;

private:
  typedef MinOperations<value_type, Cmb, value_compare> minproblem_type;
  typedef typename SolvableMarginalizer<value_type, solution_type>::varstep_pair varstep_pair;
  typedef typename SolvableMarginalizer<value_type, solution_type>::varstep_vector varstep_vector;
  typedef std::pair<value_type, DomIndex> valueindex_pair;
  typedef std::vector<valueindex_pair> valueindex_vector;
  typedef typename valueindex_vector::const_iterator valueindex_iterator;

  varstep_vector inVarsSteps_;
  Var outVar_;
  DomIndex outDomSize_;
  valueindex_vector solveVector_;

  virtual value_type marginalizeImpl(std::size_t outIndex, const table_type& mrgTable) {
    value_type minValue = *std::min_element(mrgTable.begin(), mrgTable.end(), value_compare());

    typename valueindex_vector::iterator svBegin = solveVector_.begin() + outIndex * outDomSize_;
    typename valueindex_vector::iterator svEnd = svBegin + outDomSize_;

    typename valueindex_vector::iterator svIter = svBegin;
    DomIndex index = 0;
    for (const auto &v: mrgTable) {
      svIter->first = minproblem_type::combineInverse(v, minValue);
      svIter->second = index++;
      ++svIter;
    }

    std::sort(svBegin, svEnd);
    return minValue;
  }

  virtual void solveImpl(solution_type& solSet) const {
    typedef typename solution_type::solution_type single_solution_type;
    typedef typename solution_type::solution_set solution_set;

    solution_set inSolSet;
    inSolSet.swap(solSet.solutions());
    solution_set& outSolSet = solSet.solutions();
    solution_compare solLess;

    for (auto sol: inSolSet) {
      bool added = false;
      value_type baseValue = sol.value;

      std::size_t svIndex = 0;
      for (const auto &vs: inVarsSteps_) {
        svIndex += sol.solution[vs.first] * vs.second;
      }
      svIndex *= outDomSize_;
      valueindex_iterator svBegin = solveVector_.begin() + svIndex;
      valueindex_iterator svEnd = svBegin + outDomSize_;

      for (auto it = svBegin; it != svEnd; ++it) {
        sol.solution[outVar_] = it->second;
        sol.value = minproblem_type::combine(baseValue, it->first);

        if (outSolSet.size() < solSet.maxSolutions()) {
          outSolSet.insert(sol);
        } else {
          typename solution_set::iterator last = outSolSet.end();
          --last;
          if (solLess(sol, *last)) {
            outSolSet.insert(sol);
            outSolSet.erase(last);
          } else {
            break;
          }
        }

        added = true;
      }

      if (!added) {
        break;
      }
    }
  }

public:
  SolvableMinMarginalizer(
      const VarVector& inScope, const DomIndexVector& inDomSizes,
      Var outVar, DomIndex outDomSize) :
    inVarsSteps_(),
    outVar_(outVar),
    outDomSize_(outDomSize),
    solveVector_() {

    if (inScope.size() != inDomSizes.size()) {
      throw InvalidArgumentException();
    }

    size_t numInEntries = SolvableMinMarginalizer::buildStepSizes(inScope, inDomSizes, inVarsSteps_);
    solveVector_.resize(numInEntries * outDomSize);
  }
};

template<typename Y, typename Cmb, typename Compare=std::less<Y> >
class MinOperations : public Cmb {
public:
  typedef Y value_type;
  typedef Compare value_compare;
  typedef MinSolutionSet<value_type, value_compare> solution_type;
  typedef MarginalizerTypes<value_type, solution_type> marginalizer_types;
  typedef typename marginalizer_types::marginalizer_type marginalizer_type;
  typedef typename marginalizer_types::solvablemarginalizer_type solvablemarginalizer_type;
  typedef typename marginalizer_types::marginalizer_smartptr marginalizer_smartptr;
  typedef typename marginalizer_types::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;

  struct CtorArgs {
    std::size_t maxSolutions;
    CtorArgs(std::size_t maxSolutions0 = 1) : maxSolutions(maxSolutions0) {}
  };

private:
  typedef Cmb combiner_type;

  marginalizer_smartptr marginalizer_;
  std::size_t maxSolutions_;

public:
  using combiner_type::combineIdentity;
  using combiner_type::combine;
  using combiner_type::combineInverse;

  explicit MinOperations(const CtorArgs& ctorArgs = CtorArgs()) :
    marginalizer_(new MinMarginalizer<value_type, value_compare>()),
    maxSolutions_(ctorArgs.maxSolutions) {}

  marginalizer_smartptr marginalizer() const {
    return marginalizer_;
  }

  solvablemarginalizer_smartptr solvableMarginalizer(
      const VarVector& inScope, const DomIndexVector& inDomSizes,
      Var outVar, DomIndex outDomSize) const {
    return solvablemarginalizer_smartptr(
        new SolvableMinMarginalizer<value_type, Cmb, value_compare>(inScope, inDomSizes, outVar, outDomSize));
  }

  solution_type initSolution(const DomIndexVector& x0) const {
    MinSolutionSet<value_type> initSol(maxSolutions_);
    initSol.solutions().insert(MinSolution<value_type>(combineIdentity(), x0));
    return initSol;
  }

  MinOperations& maxSolutions(std::size_t maxSols) {
    maxSolutions_ = maxSols;
    return *this;
  }

  std::size_t maxSolutions() const { return maxSolutions_; }
};

} // namespace orang

#endif
