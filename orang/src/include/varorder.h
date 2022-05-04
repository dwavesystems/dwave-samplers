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
#ifndef INCLUDED_ORANG_VARORDER_H
#define INCLUDED_ORANG_VARORDER_H

#include <cstddef>
#include <cmath>
#include <algorithm>
#include <set>
#include <vector>
#include <utility>
#include <iterator>
#include <limits>
#include <memory>
#include <cassert>

#include <base.h>
#include <exception.h>
#include <graph.h>
#include <task.h>

namespace orang {

namespace greedyvarorder {
namespace internal {

struct Variable;
typedef std::shared_ptr<Variable> var_ptr;


struct Variable {
  const Var index;
  const double domSize;
  bool processed;
  int clampRank;
  double clampValue;
  double cost;
  double complexity;
  VarSet adjList;

  Variable(Var index0, const TaskBase& task, const std::vector<int>& clampRanks) :
    index(index0),
    domSize(task.domSize(index0)),
    processed(clampRanks[index0] < 0),
    clampRank(clampRanks[index0]),
    clampValue(),
    cost(),
    complexity(),
    adjList() {

    const Graph& g = task.graph();

    for (auto it = g.adjacencyBegin(index); it != g.adjacencyEnd(index); ++it) {
      if (clampRanks[*it] >= 0) {
        adjList.insert(*it);
      }
    }
  }

  Variable(Var index0, double domSize0, bool processed0, int clampRank0, double clampValue0,
      double cost0, double complexity0) :
        index(index0), domSize(domSize0), processed(processed0), clampRank(clampRank0), clampValue(clampValue0),
        cost(cost0), complexity(complexity0), adjList() {}

  static var_ptr upperBound(const Variable& var) {
    return std::make_shared<Variable>(std::numeric_limits<Var>::max(), var.domSize, var.processed, var.clampRank, var.clampValue,
                                      var.cost, var.complexity);
  }

  static var_ptr complexityUpperBound(double maxComplexity) {
    return std::make_shared<Variable>(std::numeric_limits<Var>::max(), 0.0, false, 0, 0.0,
                                      std::numeric_limits<double>::infinity(), maxComplexity);
  }

  static var_ptr clampRankUpperBound(int rank) {
    return std::make_shared<Variable>(std::numeric_limits<Var>::max(), 0.0, false, rank,
                                      -std::numeric_limits<double>::infinity(), 0.0, 0.0);
  }
};



//===================================================================================================================
//
//   C O M P A R I S O N   F U N C T O R S
//
//===================================================================================================================

/*
 * Var objects are sorted as follows:
 * 1. Processed variables appear last and are not sorted further.
 * 2. Within unprocessed variables, those whose complexity exceeds the maximum appear last and are not sorted further.
 * 3. Below-complexity-limit variables are sorted by increasing cost, with ties broken by variable index.
 */
class CostCmp {
private:
  double maxComplexity_;

public:
  CostCmp(double maxComplexity) : maxComplexity_(maxComplexity) {}

  bool operator()(const var_ptr& v1, const var_ptr& v2) const {
    return !v1->processed
        && (v2->processed
            || (v1->complexity <= maxComplexity_
                && (v2->complexity > maxComplexity_
                    || v1->cost < v2->cost
                    || (v1->cost == v2->cost && v1->index < v2->index))));
  }
};

struct ClampCmp {
  bool operator()(const var_ptr& v1, const var_ptr& v2) const {
    return !v1->processed
        && (v2->processed
            || v1->clampRank < v2->clampRank
            || (v1->clampRank == v2->clampRank
                && (v1->clampValue > v2->clampValue
                    || (v1->clampValue == v2->clampValue && v1->index < v2->index))));
  }
};



//===================================================================================================================
//
//   M U L T I - I N D E X   V A R I A B L E   C O N T A I N E R
//
//===================================================================================================================


/*
 // * Container of variables that provides multiple modes of access.  Random
 // * access is provided through a vector, while two sortings are also
 // * maintained through multisets.  Each internal container stores smart
 // * pointers to Variable instances.  Functions are provided for modifying an
 // * element specified by an iterator of any particular container.  Internally,
 // * modification is done by removing and adding back elements to the
 // * multisets, to preserve the order.
 */
class VarContainer {
public:
  std::vector<var_ptr> byIndex;
  std::multiset<var_ptr, CostCmp> byCost;
  std::multiset<var_ptr, ClampCmp> byClamp;

  VarContainer(const TaskBase& task, double maxComplexity, const std::vector<int>& clampRank)
    : byCost(CostCmp(maxComplexity)) {
    const Graph& g = task.graph();
    Var numVertices = g.numVertices();

    for (Var v = 0; v < numVertices; ++v) {
      add(std::make_shared<Variable>(v, task, clampRank));
    }
  }

  void add(var_ptr var) {
    byIndex.push_back(var);
    byCost.insert(var);
    byClamp.insert(var);
  }

  template<typename F>
  void modifyByIndex(std::vector<var_ptr>::iterator it, F& func) {
    func(**it);

    auto pos = find(byCost.begin(), byCost.end(), *it);
    assert(pos != byCost.end());
    byCost.erase(pos);
    byCost.insert(*it);

    pos = find(byClamp.begin(), byClamp.end(), *it);
    assert(pos != byClamp.end());
    byClamp.erase(pos);
    byClamp.insert(*it);
  }

  template<typename F>
  void modifyByCost(std::multiset<var_ptr, CostCmp>::iterator it, F func) {
    auto pos_index = find(byIndex.begin(), byIndex.end(), *it);
    func(**pos_index);

    byCost.erase(it);
    byCost.insert(*pos_index);

    auto pos = find(byClamp.begin(), byClamp.end(), *it);
    assert(pos != byClamp.end());
    byClamp.erase(pos);
    byClamp.insert(*pos_index);
  }

  template<typename F>
  void modifyByClamp(std::multiset<var_ptr, ClampCmp>::iterator it, F func) {
    auto pos_index = find(byIndex.begin(), byIndex.end(), *it);
    func(**pos_index);

    byClamp.erase(it);
    byClamp.insert(*pos_index);

    auto pos = find(byCost.begin(), byCost.end(), *it);
    assert(pos != byCost.end());
    byCost.erase(pos);
    byCost.insert(*pos_index);
  }

};



//===================================================================================================================
//
//   M O D I F I E R   F U N C T O R S
//
//===================================================================================================================


struct MarkAsProcessed {
  void operator()(Variable& var) const {
    var.processed = true;
  }
};

struct DecrementClampRank {
  void operator()(Variable& var) const {
    --var.clampRank;
  }
};

class ElimNeighbour {
private:
  const Var elimVar_;
  const VarSet& vars_;

public:
  ElimNeighbour(Var elimVar, const VarSet& vars = VarSet()) : elimVar_(elimVar), vars_(vars) {}
  void operator()(Variable& var) const {
    var.adjList.insert(vars_.begin(), vars_.end());
    var.adjList.erase(var.index);
    var.adjList.erase(elimVar_);
  }
};

class ClampNeighbour {
private:
  const Var clampVar_;

public:
  ClampNeighbour(Var clampVar) : clampVar_(clampVar) {}
  void operator()(Variable& var) const {
    var.adjList.erase(clampVar_);
  }
};



//===================================================================================================================
//
//   H E U R I S T I C - S P E C I F I C   S T U F F
//
//===================================================================================================================


//-------------------------------------------------------------------------------------------------------------------
// Var member data modifier functors
//-------------------------------------------------------------------------------------------------------------------

/*
 * These values are based on current contents of the variable's adjList and the adjList contents of its neighbours;
 * thus, this functor must be applied to all appropriate variables AFTER UpdateNeighbours has been applied to ALL
 * those same variables.
 *
 * Cost calculation is heuristic-dependent.  Derived classes exist for the different calculations.
 */
class UpdateVarData {
private:
  virtual void updateCost(Variable& var) const = 0;

protected:
  const std::vector<var_ptr>& varsByIndex_;

public:
  UpdateVarData(const std::vector<var_ptr>& varsByIndex) : varsByIndex_(varsByIndex) {}
  virtual ~UpdateVarData() {}

  void operator()(Variable& var) const {
    var.clampValue = static_cast<double>(var.domSize) * static_cast<double>(var.adjList.size());
    double p2Cplx = var.domSize;
    for (const auto &w: var.adjList) {
      p2Cplx *= varsByIndex_[w]->domSize;
    }
    static const double E_LOG2 = 1.4426950408889633;
    var.complexity = log(p2Cplx) * E_LOG2;
    updateCost(var);
  }
};

class UpdateMinDegreeVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = static_cast<double>(var.adjList.size());
  }
public:
  UpdateMinDegreeVarData(std::vector<var_ptr>& varsByIndex) : UpdateVarData(varsByIndex) {}
};

class UpdateWeightedMinDegreeVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = var.clampValue;
  }
public:
  UpdateWeightedMinDegreeVarData(const std::vector<var_ptr>& varsByIndex) : UpdateVarData(varsByIndex) {}
};

class UpdateMinFillVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = 0.0;
    for (VarSet::const_iterator vAdjIter = var.adjList.begin(), vAdjEnd = var.adjList.end();
        vAdjIter != vAdjEnd; ++vAdjIter) {
      const Var u = *vAdjIter;
      const Variable* uVar = varsByIndex_[u].get();
      VarSet::const_iterator uAdjIter = uVar->adjList.upper_bound(u);
      VarSet::const_iterator uAdjEnd = uVar->adjList.end();
      VarSet::const_iterator vAdjIter2 = vAdjIter;
      ++vAdjIter2;
      while (vAdjIter2 != vAdjEnd) {
        if (uAdjIter == uAdjEnd || *vAdjIter2 < *uAdjIter) {
          ++var.cost;
          ++vAdjIter2;
        } else if (*uAdjIter < *vAdjIter2) {
          ++uAdjIter;
        } else {
          ++vAdjIter2;
          ++uAdjIter;
        }
      }
    }
  }
public:
  UpdateMinFillVarData(const std::vector<var_ptr>& varsByIndex) : UpdateVarData(varsByIndex) {}
};

class UpdateWeightedMinFillVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = 0.0;
    for (VarSet::const_iterator vAdjIter = var.adjList.begin(), vAdjEnd = var.adjList.end();
        vAdjIter != vAdjEnd; ++vAdjIter) {
      const Var u = *vAdjIter;
      const Variable* uVar = varsByIndex_[u].get();
      VarSet::const_iterator uAdjIter = uVar->adjList.upper_bound(u);
      VarSet::const_iterator uAdjEnd = uVar->adjList.end();
      VarSet::const_iterator vAdjIter2 = vAdjIter;
      ++vAdjIter2;
      double cost = 0.0;
      while (vAdjIter2 != vAdjEnd) {
        if (uAdjIter == uAdjEnd || *vAdjIter2 < *uAdjIter) {
          cost += varsByIndex_[*vAdjIter2]->domSize;
          ++vAdjIter2;
        } else if (*uAdjIter < *vAdjIter2) {
          ++uAdjIter;
        } else {
          ++vAdjIter2;
          ++uAdjIter;
        }
      }
      var.cost += uVar->domSize * cost;
    }
  }
public:
  UpdateWeightedMinFillVarData(const std::vector<var_ptr>& varsByIndex) : UpdateVarData(varsByIndex) {}
};


//-------------------------------------------------------------------------------------------------------------------
// List-of-affected-variables functors
//-------------------------------------------------------------------------------------------------------------------

class AffectedVars {
private:
  virtual VarSet affectedVars(const Variable&) const = 0;
public:
  virtual ~AffectedVars() {}
  VarSet operator()(const Variable& var) const {
    return affectedVars(var);
  }
};

class MinDegreeAffectedVars : public AffectedVars {
private:
  virtual VarSet affectedVars(const Variable& var) const {
    return var.adjList;
  }
};

class MinFillAffectedVars : public AffectedVars {
private:
  const std::vector<var_ptr>& varsByIndex_;
  virtual VarSet affectedVars(const Variable& var) const {
    VarSet vars = var.adjList;
    for (const auto &u: var.adjList) {
      const Variable* uVar = varsByIndex_[u].get();
      vars.insert(uVar->adjList.begin(), uVar->adjList.end());
    }
    vars.erase(var.index);
    return vars;
  }

public:
  MinFillAffectedVars(const VarContainer& varContainer) : varsByIndex_(varContainer.byIndex) {}
};



//===================================================================================================================
//
//   R A N D O M   V A R I A B L E   S E L E C T O R
//
//===================================================================================================================

template<typename Iter, typename Rng>
Iter selectVar(Iter begin, Iter baseEnd, Iter finalEnd, Rng& rng, float selectionScale) {

  float baseRange = static_cast<float>(std::distance(begin, baseEnd));
  float totalRange = baseRange + static_cast<float>(std::distance(baseEnd, finalEnd));
  float selectionRange = std::min(baseRange * selectionScale, totalRange);
  float incr = std::floor(static_cast<float>(selectionRange * rng()));
  incr = std::max(incr, 0.0f);
  incr = std::min(incr, totalRange - 1);

  Iter ret = begin;
  std::advance(ret, incr);
  return ret;
}

} // namespace orang::greedyvarorder::internal




//===================================================================================================================
//
//   H E U R I S T I C   E N U M
//
//===================================================================================================================

enum Heuristics {
  MIN_DEGREE,
  WEIGHTED_MIN_DEGREE,
  MIN_FILL,
  WEIGHTED_MIN_FILL,
  NUM_HEURISTICS
};

} // namespace orang::greedyvarorder



//===================================================================================================================
//
//   T H E   F U N C T I O N
//
//===================================================================================================================


template<typename Rng>
VarVector greedyVarOrder(
    const TaskBase& task,
    double maxComplexity,
    const std::vector<int>& clampRank,
    greedyvarorder::Heuristics h,
    Rng& rng,
    float selectionScale = 1.0f) {

  using std::floor;
  using std::advance;
  using std::distance;
  using namespace greedyvarorder::internal;

  if (task.numVars() != clampRank.size()) {
    throw InvalidArgumentException("clampRank size must equal the number of variables in task");
  }

  if (task.numVars() == 0) {
    return VarVector();
  }

  VarContainer vars(task, maxComplexity, clampRank);

  std::unique_ptr<UpdateVarData> updateCostPtr;
  std::unique_ptr<AffectedVars> affectedVarsPtr;
  switch (h) {
    case greedyvarorder::MIN_DEGREE:
      updateCostPtr.reset( new UpdateMinDegreeVarData(vars.byIndex) );
      affectedVarsPtr.reset( new MinDegreeAffectedVars() );
      break;
    case greedyvarorder::WEIGHTED_MIN_DEGREE:
      updateCostPtr.reset( new UpdateWeightedMinDegreeVarData(vars.byIndex) );
      affectedVarsPtr.reset( new MinDegreeAffectedVars() );
      break;
    case greedyvarorder::MIN_FILL:
      updateCostPtr.reset( new UpdateMinFillVarData(vars.byIndex) );
      affectedVarsPtr.reset( new MinFillAffectedVars(vars) );
      break;
    case greedyvarorder::WEIGHTED_MIN_FILL:
      updateCostPtr.reset( new UpdateWeightedMinFillVarData(vars.byIndex) );
      affectedVarsPtr.reset( new MinFillAffectedVars(vars) );
      break;
    default:
      throw InvalidArgumentException("Invalid heuristic");
  }

  for (auto iter = vars.byIndex.begin(); iter != vars.byIndex.end(); ++iter) {
    vars.modifyByIndex(iter, *updateCostPtr);
  }

  VarVector varOrder;
  int lastClampRank = -1;
  const var_ptr complexityUpper = Variable::complexityUpperBound(maxComplexity);

  for (;;) {

    auto minCostLower = vars.byCost.begin();
    if ((*minCostLower)->processed) {
      break;
    }

    if ((*minCostLower)->complexity <= maxComplexity) {
      auto pickedIter = selectVar(minCostLower, vars.byCost.upper_bound( Variable::upperBound(**minCostLower)),
          vars.byCost.upper_bound(complexityUpper), rng, selectionScale);

      const Variable& v = **pickedIter;
      varOrder.push_back(v.index);
      VarSet affectedVars = (*affectedVarsPtr)(v);
      ElimNeighbour elimNeighbour(v.index, v.adjList);

      vars.modifyByCost(pickedIter, MarkAsProcessed());

      for (const auto &uIndex: v.adjList) {
        vars.modifyByIndex(vars.byIndex.begin() + uIndex, elimNeighbour);
      }

      for (const auto &uIndex: affectedVars) {
        vars.modifyByIndex(vars.byIndex.begin() + uIndex, *updateCostPtr);
      }

    } else {
      if (lastClampRank >= 0) {
        auto clampIter = vars.byClamp.upper_bound( Variable::clampRankUpperBound(lastClampRank) );
        auto clampEnd = vars.byClamp.end();
        while (clampIter != clampEnd && !(*clampIter)->processed) {
          auto here = clampIter++;
          vars.modifyByClamp(here, DecrementClampRank());
        }
      }

      auto clampLower = vars.byClamp.begin();
      auto pickedIter = selectVar(clampLower,
                                  vars.byClamp.upper_bound( Variable::upperBound(**clampLower) ),
                                  vars.byClamp.upper_bound( Variable::clampRankUpperBound((*clampLower)->clampRank) ),
                                  rng, selectionScale);

      const Variable& v = **pickedIter;
      lastClampRank = v.clampRank;
      vars.modifyByClamp(pickedIter, MarkAsProcessed());
      ClampNeighbour clampNeighbour(v.index);

      for (const auto &uIndex: v.adjList) {
        vars.modifyByIndex(vars.byIndex.begin() + uIndex, clampNeighbour);
      }

      for (const auto &uIndex: v.adjList) {
        vars.modifyByIndex(vars.byIndex.begin() + uIndex, *updateCostPtr);
      }
    }
  }

  return varOrder;
}

} // namespace orang

#endif
