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
#ifndef INCLUDED_ORANG_MERGER_H
#define INCLUDED_ORANG_MERGER_H

#include <cstddef>
#include <vector>
#include <memory>
#include <iterator>
#include <algorithm>
#include <functional>
#include <utility>
#include <memory>

#include <base.h>
#include <table.h>
#include <task.h>
#include <marginalizer.h>

namespace orang {

template<typename T>
class TableMerger {
public:
  typedef T task_type; ///< Task
  typedef typename task_type::table_type table_type;
  typedef typename task_type::table_smartptr table_smartptr;
  typedef typename task_type::table_vector table_vector;
  typedef typename task_type::marginalizer_type marginalizer_type;
private:
  const task_type& task_;

  static void constraints() {
    DomIndex (T::*ds)(Var) const = &T::domSize;
    typedef typename T::value_type value_type;
    value_type (*cmb)(const value_type&, const value_type&) = &T::combine;
    static_cast<void>(ds);
    static_cast<void>(cmb);
  }

public:
  TableMerger(const task_type& task) : task_(task) {
    void (*cs)() = &TableMerger::constraints;
    static_cast<void>(cs);
  }

  template<typename TblIter>
  table_smartptr operator()(
      const VarVector& outScope,
      const TblIter& tablesBegin,
      const TblIter& tablesEnd,
      marginalizer_type& marginalizer) const;
};



namespace internal {

template<typename Tbl>
struct TableVarIter {
  typedef typename std::vector<TableVar>::const_iterator var_iterator;
  typedef typename Tbl::const_iterator tbl_iterator;

  var_iterator varIter;
  var_iterator varEnd;
  tbl_iterator* tableIter;

  TableVarIter(var_iterator varIter0, var_iterator varEnd0, tbl_iterator* tableIter0) :
    varIter(varIter0), varEnd(varEnd0), tableIter(tableIter0) {}
};

template<typename Iter>
class StepIter {
private:
  Iter* iter_;
  std::size_t stepSize_;
public:
  StepIter(Iter* iter, std::size_t stepSize) : iter_(iter), stepSize_(stepSize) {}
  void operator+=(int n) {
    *iter_ += n * stepSize_;
  }
};

template<typename T>
class GrayVar {
public:
  typedef typename T::table_type table_type;
  typedef typename table_type::const_iterator table_const_iterator;
  typedef std::vector<StepIter<table_const_iterator> > tableiter_vector;

private:

  DomIndex domIndex_;
  int dir_;
  const DomIndex domSize_;

  tableiter_vector inIters_;
  StepIter<std::size_t> outIndex_;

public:
  GrayVar(
      DomIndex domSize,
      std::size_t* outIndex,
      std::size_t outDelta) :
        domIndex_(0),
        dir_(1),
        domSize_(domSize),
        inIters_(),
        outIndex_(outIndex, outDelta) {}

  void addInIter(table_const_iterator* inIter, std::size_t stepSize) {
    inIters_.push_back(StepIter<table_const_iterator>(inIter, stepSize));
  }

  bool advance() {
    DomIndex nextDomIndex = domIndex_ + dir_;

    // nextDomIndex is unsigned, so things >= domSize_
    // are precisely the invalid values, regardless of dir_
    if (nextDomIndex < domSize_) {
      for (auto &inIter: inIters_) {
        inIter += dir_;
      }
      outIndex_ += dir_;
      domIndex_ = nextDomIndex;
      return true;

    } else {
      dir_ = -dir_;
      return false;
    }
  }
};

}



template<typename T>
template<typename TblIter>
typename TableMerger<T>::table_smartptr
TableMerger<T>::operator ()(
    const VarVector& outScope,
    const TblIter& tablesBegin,
    const TblIter& tablesEnd,
    typename TableMerger<T>::marginalizer_type& marginalize) const {

  //===========================================================================================================
  // typedef and using declarations

  using std::size_t;
  using std::vector;
  using std::find_if;
  using std::distance;
  using std::make_pair;
  using std::unique_ptr;

  using internal::GrayVar;
  using internal::TableVarIter;

  typedef typename TableMerger<T>::table_type table_type;
  typedef typename table_type::const_iterator table_const_iterator;
  typedef typename table_type::smartptr table_smartptr;

  typedef typename vector<TableVar>::const_iterator tablevar_const_iterator;

  typedef TableVarIter<table_type> tablevariter_type;
  typedef vector<tablevariter_type> tablevariter_vector;

  typedef std::vector<unique_ptr<GrayVar<T>>> grayvar_ptrvector;

  //===========================================================================================================
  //
  //  Validate input
  //

  size_t numTables = distance(tablesBegin, tablesEnd);
  if (numTables == 0) {
    table_smartptr t( new table_type(VarVector(), DomIndexVector()) );
    (*t)[0] = task_type::combineIdentity();
    return t;
  }

  for (auto it = tablesBegin; it != tablesEnd; ++it) {
    for (const auto &v: (*it)->vars()) {
      if (v.domSize != task_.domSize(v.index)) {
        throw InvalidArgumentException("Table and Task domain sizes don't match");
      }
    }
  }

  //===========================================================================================================
  //
  //  Variables needed throughout this function
  //

  // inTable iterators.  For each inTable element, there will be an iterator pointer (with corresponding step
  // sizes) for each of its variables.  These pointers (for a single input table) all refer to the same
  // iterator which must exist somewhere.  That somewhere is here.
  vector<table_const_iterator> inTableIters;
  inTableIters.reserve(numTables);

  // output table, domain sizes, and index.
  DomIndexVector outDomSizes;
  outDomSizes.reserve(outScope.size());
  for (auto var: outScope) {
    outDomSizes.push_back(task_.domSize(var));
  }
  table_smartptr outTable( new table_type(outScope, outDomSizes));
  table_type& outTableRef = *outTable.get();
  size_t outIndex = 0;

  // marginalization table scope, domain sizes and index.  These are built during Phase 1
  VarVector mrgScope;
  DomIndexVector mrgDomSizes;
  size_t mrgIndex = 0;

  // Gray-counting vectors.  One for variables in outScope, one for all other variables appearing in inTables.
  grayvar_ptrvector outGrayVars;
  grayvar_ptrvector mrgGrayVars;
  outGrayVars.reserve(outScope.size());

  //===========================================================================================================
  //
  //  Phase 1: construct the GrayVar vectors for iterating over output variables and
  //           marginalized-out input variables.
  //

  // First pass: initialize TableVarIter vector and determine first (ie. lowest-index) variable to process.
  Var nextVar = outScope.empty() ? Var(-1) : outScope.front();
  tablevariter_vector tableVarIters;
  for (auto it = tablesBegin; it != tablesEnd; ++it) {
    inTableIters.push_back((*it)->begin());
    tableVarIters.push_back(tablevariter_type((*it)->vars().begin(), (*it)->vars().end(), &inTableIters.back()));

    Var tableVar0 = (*it)->vars().empty() ? Var(-1) : (*it)->var(0).index;
    nextVar = tableVar0 < nextVar ? tableVar0 : nextVar;
  }

  // Now iterate over all variables in all tables (input and output) and build GrayVar vectors
  tablevar_const_iterator outVarIter = outTableRef.vars().begin();
  tablevar_const_iterator outVarEnd = outTableRef.vars().end();
  size_t mrgDelta = 1;
  for (bool done = false; !done; ) {
    done = true;
    Var currentVar = nextVar;
    GrayVar<T>* grayVar;

    // Determine which GrayVar vector is to receive the new GrayVar.
    if (outVarIter != outVarEnd && currentVar == outVarIter->index) {
      outGrayVars.push_back(unique_ptr<GrayVar<T>>(new GrayVar<T>(outVarIter->domSize, &outIndex, outVarIter->stepSize)));
      grayVar = outGrayVars.back().get();
      ++outVarIter;
    } else {
      DomIndex domSize = task_.domSize(currentVar);
      mrgScope.push_back(currentVar);
      mrgDomSizes.push_back(domSize);
      mrgGrayVars.push_back(unique_ptr<GrayVar<T>>(new GrayVar<T>(domSize, &mrgIndex, mrgDelta)));
      grayVar = mrgGrayVars.back().get();
      mrgDelta *= domSize;
    }

    if (outVarIter != outVarEnd) {
      nextVar = outVarIter->index;
      done = false;
    } else {
      nextVar = Var(-1);
    }

    // Add all inTableIters (with appropriate step sizes) for tables containing currentVar to grayVar
    for (auto &tvi: tableVarIters) {
      if (tvi.varIter != tvi.varEnd && tvi.varIter->index == currentVar) {
        grayVar->addInIter(tvi.tableIter, tvi.varIter->stepSize);
        ++tvi.varIter;
      }

      if (tvi.varIter != tvi.varEnd) {
        nextVar = nextVar < tvi.varIter->index ? nextVar : tvi.varIter->index;
        done = false;
      }
    }
  }

  //===========================================================================================================
  //
  //  Phase 2: Loop through output table entries using outGrayVars.  At each step, loop through variables
  //           to be marginalized out using mrgGrayVars, combining entries from all input tables.
  //

  // construct table for combining input table values.  Scope is variables to be marginalized out.
  table_type mrgTable((mrgScope), mrgDomSizes);

  // now loop over output variables
  do {

    // build table of marginalized-out variable values for current setting of output variables
    do {
      mrgTable[mrgIndex] = *inTableIters.front();
      for (auto it = inTableIters.begin() + 1; it != inTableIters.end(); ++it) {
        mrgTable[mrgIndex] = task_type::combine(mrgTable[mrgIndex], **it);
      }
    } while (find_if(mrgGrayVars.begin(), mrgGrayVars.end(), std::mem_fn(&GrayVar<T>::advance)) != mrgGrayVars.end());

    outTableRef[outIndex] = marginalize(outIndex, mrgTable);

    // loop condition: GrayVar<T>::advance returns true iff its value changed.  Continue until no value changes
  } while (find_if(outGrayVars.begin(), outGrayVars.end(), std::mem_fn(&GrayVar<T>::advance)) != outGrayVars.end());

  //===========================================================================================================
  //
  //  Done!
  //
  return outTable;
}

}

#endif
