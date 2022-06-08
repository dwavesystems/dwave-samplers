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
#ifndef INCLUDED_ORANG_TASK_H
#define INCLUDED_ORANG_TASK_H

#include <cstddef>
#include <utility>
#include <vector>
#include <set>
#include <memory>

#include <base.h>
#include <graph.h>
#include <treedecomp.h>
#include <marginalizer.h>

namespace orang {

class TaskBase {
protected:
  DomIndexVector domSizes_;
  Graph graph_;

  TaskBase() : domSizes_(), graph_() {}

public:
  const DomIndexVector& domSizes() const { return domSizes_; }
  DomIndex domSize(Var v) const { return domSizes_[v]; }
  Var numVars() const { return static_cast<Var>(domSizes_.size()); }
  const Graph& graph() const { return graph_; }
};

template<typename Ops>
class Task : public TaskBase, public Ops {
public:
  typedef Ops ops_type;
  typedef typename ops_type::value_type value_type;
  typedef typename ops_type::solution_type solution_type;
  typedef typename ops_type::marginalizer_type marginalizer_type;
  typedef typename ops_type::solvablemarginalizer_type solvablemarginalizer_type;
  typedef typename ops_type::marginalizer_smartptr marginalizer_smartptr;
  typedef typename ops_type::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;
  typedef typename ops_type::CtorArgs problem_ctorargs;

  typedef Table<value_type> table_type;
  typedef std::shared_ptr<const table_type> const_table_smartptr;
  typedef std::shared_ptr<table_type> table_smartptr;
  typedef std::vector<const_table_smartptr> table_vector;

private:
  table_vector tables_;

public:
  using ops_type::combineIdentity;
  using ops_type::combine;
  using ops_type::marginalizer;
  using ops_type::solvableMarginalizer;
  using ops_type::initSolution;

  template<typename Iter>
  Task(
      Iter tablesBegin,
      Iter tablesEnd,
      const problem_ctorargs& mrgCtorArgs,
      Var minVars = 0
      ) :
        TaskBase(),
        ops_type(mrgCtorArgs),
        tables_() {

    using std::vector;
    using std::make_pair;
    using std::max;

    typedef vector<TableVar> tablevar_vector;
    typedef typename tablevar_vector::const_iterator tablevar_const_iterator;

    std::set<Graph::adj_pair> adjSet;
    Var numV = 0;

    for ( ; tablesBegin != tablesEnd; ++tablesBegin) {
      const_table_smartptr tableCopy( new table_type(**tablesBegin) );
      tables_.push_back(tableCopy);

      if (!tableCopy->vars().empty()) {
        numV = max(numV, tableCopy->vars().back().index + 1);
      }

      for (tablevar_const_iterator varsIter = tableCopy->vars().begin(), varsEnd = tableCopy->vars().end();
          varsIter != varsEnd; ++varsIter) {

        size_t dsSize = domSizes_.size();
        if (varsIter->index >= dsSize) {
          domSizes_.resize(varsIter->index + 1, 1);
          domSizes_.back() = varsIter->domSize;
        } else if (domSizes_[varsIter->index] == 1) {
          domSizes_[varsIter->index] = varsIter->domSize;
        } else if (domSizes_[varsIter->index] != varsIter->domSize) {
          throw InvalidArgumentException("Tables contain inconsistent domain sizes");
        }

        for (tablevar_const_iterator varsIter2 = varsIter + 1; varsIter2 != varsEnd; ++ varsIter2) {
          adjSet.insert(make_pair(varsIter->index, varsIter2->index));
        }
      }
    }

    if (domSizes_.size() < minVars) {
      domSizes_.resize(minVars, 1);
    }
    graph_.setAdjacencies(adjSet, max(numV, minVars));
  }

  const table_vector& tables() const { return tables_; }

  table_vector baseTables(const TreeDecompNode& dNode, const DomIndexVector& x0) const;

  value_type problemValue(
      const std::vector<value_type>& rootValues,
      const DomIndexVector& x0,
      const VarVector& clampedVars) const {

    value_type v = combineIdentity();
    for (const auto &rv: rootValues) {
      v = combine(v, rv);
    }

    VarVector::const_iterator cvEnd = clampedVars.end();
    for (const auto &tp: tables_) {
      bool allClamped = true;
      typename table_type::const_iterator tblIter = tp->begin();

      VarVector::const_iterator cvIter = clampedVars.begin();
      for (const auto &tv: tp->vars()) {
        while (cvIter != cvEnd && *cvIter < tv.index) {
          ++cvIter;
        }

        if (cvIter == cvEnd || *cvIter > tv.index) {
          allClamped = false;
          break;
        }

        tblIter += tv.stepSize * x0[tv.index];
      }

      if (allClamped) {
        v = combine(v, *tblIter);
      }
    }

    return v;
  }
};


template<typename Ops>
typename Task<Ops>::table_vector
Task<Ops>::baseTables(
    const TreeDecompNode& dNode,
    const DomIndexVector& x0) const {

  using std::size_t;
  using std::vector;

  typedef vector<TableVar> tablevar_vector;
  typedef typename table_type::const_iterator table_const_iterator;

  table_vector retTables;

  for (const auto &t: tables_) {

    // check if t's scope contains dNode.nodeVar and is otherwise contained in dNode.sepVars + dNodes.clampedVars
    bool goodScope = false;
    VarVector newScope;
    DomIndexVector newDomSizes;
    SizeVector oldStepSizes;
    newScope.reserve(t->vars().size());
    table_const_iterator tIter = t->begin();

    typename tablevar_vector::const_iterator varsIter = t->vars().begin();
    typename tablevar_vector::const_iterator varsEnd = t->vars().end();
    VarVector::const_iterator sepVarsIter = dNode.sepVars().begin();
    VarVector::const_iterator sepVarsEnd = dNode.sepVars().end();
    VarVector::const_iterator clampedVarsIter = dNode.clampedVars().begin();
    VarVector::const_iterator clampedVarsEnd = dNode.clampedVars().end();
    while (varsIter != varsEnd) {
      bool sepEq = sepVarsIter != sepVarsEnd && *sepVarsIter == varsIter->index;
      bool sepLessEq = sepVarsIter != sepVarsEnd && *sepVarsIter <= varsIter->index;
      bool clampedEq = clampedVarsIter != clampedVarsEnd && *clampedVarsIter == varsIter->index;
      bool clampedLessEq = clampedVarsIter != clampedVarsEnd && *clampedVarsIter <= varsIter->index;

      if (sepEq) {
        newScope.push_back(varsIter->index);
        newDomSizes.push_back(varsIter->domSize);
        oldStepSizes.push_back(varsIter->stepSize);
        ++varsIter;
      } else if (clampedEq) {
        tIter += x0[varsIter->index] * varsIter->stepSize;
        ++varsIter;
      }

      if (sepLessEq) {
        ++sepVarsIter;
      }

      if (clampedLessEq) {
        ++clampedVarsIter;
      }

      if (!sepLessEq && !clampedLessEq) {
        if (varsIter->index == dNode.nodeVar()) {
          goodScope = true;
          newScope.push_back(varsIter->index);
          newDomSizes.push_back(varsIter->domSize);
          oldStepSizes.push_back(varsIter->stepSize);
          ++varsIter;
        } else {
          goodScope = false;
          break;
        }
      }
    }

    if (!goodScope) {
      continue;
    }

    // table scope is good.  either add it directly to the return list or remove clamped vars
    Var newScopeSize = static_cast<Var>(newScope.size());
    if (newScopeSize == t->vars().size()) {
      retTables.push_back(t);

    } else {
      std::shared_ptr<table_type> newTable( new table_type(newScope, newDomSizes) );
      typename table_type::iterator newTableIter = newTable->begin();
      vector<int> dirs(newScopeSize, 1);
      SizeVector indices(newScopeSize, 0);

      for (bool done = false; !done; ) {
        *newTableIter = *tIter;

        done = true;
        for (Var i = 0; i < newScopeSize; ++i) {
          size_t newIndex = indices[i] + dirs[i];
          if (newIndex >= newDomSizes[i]) {
            dirs[i] = -dirs[i];
          } else {
            indices[i] = newIndex;
            newTableIter += dirs[i] * newTable->var(i).stepSize;
            tIter += dirs[i] * oldStepSizes[i];
            done = false;
            break;
          }
        }
      }

      retTables.push_back( newTable );
    }
  }

  return retTables;
}


} // end of namespace orang

#endif
