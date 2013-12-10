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

#include <boost/foreach.hpp>
#include <boost/smart_ptr/scoped_ptr.hpp>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/random_access_index.hpp>
#include <boost/multi_index/tag.hpp>
#include <boost/tuple/tuple.hpp>

#include <orang/base.h>
#include <orang/exception.h>
#include <orang/graph.h>
#include <orang/task.h>

namespace orang {

namespace greedyvarorder {
namespace internal {

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

    BOOST_FOREACH( Var w, std::make_pair(g.adjacencyBegin(index), g.adjacencyEnd(index)) ) {
      if (clampRanks[w] >= 0) {
        adjList.insert(w);
      }
    }
  }

  Variable(Var index0, double domSize0, bool processed0, int clampRank0, double clampValue0,
      double cost0, double complexity0) :
        index(index0), domSize(domSize0), processed(processed0), clampRank(clampRank0), clampValue(clampValue0),
        cost(cost0), complexity(complexity0), adjList() {}

  static Variable upperBound(const Variable& var) {
    return Variable(std::numeric_limits<Var>::max(), var.domSize, var.processed, var.clampRank, var.clampValue,
        var.cost, var.complexity);
  }

  static Variable complexityUpperBound(double maxComplexity) {
    return Variable(std::numeric_limits<Var>::max(), 0.0, false, 0, 0.0,
        std::numeric_limits<double>::infinity(), maxComplexity);
  }

  static Variable clampRankUpperBound(int rank) {
    return Variable(std::numeric_limits<Var>::max(), 0.0, false, rank,
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

  bool operator()(const Variable& v1, const Variable& v2) const {
    return !v1.processed
        && (v2.processed
            || (v1.complexity <= maxComplexity_
                && (v2.complexity > maxComplexity_
                    || v1.cost < v2.cost
                    || (v1.cost == v2.cost && v1.index < v2.index))));
  }
};

struct ClampCmp {
  bool operator()(const Variable& v1, const Variable& v2) const {
    return !v1.processed
        && (v2.processed
            || v1.clampRank < v2.clampRank
            || (v1.clampRank == v2.clampRank
                && (v1.clampValue > v2.clampValue
                    || (v1.clampValue == v2.clampValue && v1.index < v2.index))));
  }
};



//===================================================================================================================
//
//   M U L T I - I N D E X   V A R I A B L E   C O N T A I N E R
//
//===================================================================================================================

struct Index {};
struct Cost {};
struct Clamp {};

struct VarIndices : public boost::multi_index::indexed_by<
  boost::multi_index::random_access<boost::multi_index::tag<Index> >,
  boost::multi_index::ordered_non_unique<boost::multi_index::tag<Cost>, boost::multi_index::identity<Variable>, CostCmp>,
  boost::multi_index::ordered_non_unique<boost::multi_index::tag<Clamp>, boost::multi_index::identity<Variable>, ClampCmp>
> {};

class VarContainer : public boost::multi_index::multi_index_container<Variable, VarIndices> {
private:
  typedef boost::multi_index::multi_index_container<Variable,VarIndices> base_type;

public:
  VarContainer(
      const TaskBase& task,
      double maxComplexity,
      const std::vector<int>& clampRank) :
          base_type(
              boost::make_tuple(
                  boost::multi_index::multi_index_container<Variable, VarIndices>::index<Index>::type::ctor_args(),
                  boost::make_tuple(boost::multi_index::identity<Variable>(), CostCmp(maxComplexity)),
                  boost::multi_index::multi_index_container<Variable, VarIndices>::index<Clamp>::type::ctor_args())) {

    const Graph& g = task.graph();
    Var numVertices = g.numVertices();

    index<Index>::type& thisByIndex = get<Index>();

    for (Var v = 0; v < numVertices; ++v) {
      Variable var(v, task, clampRank);
      thisByIndex.push_back(var);
    }
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
  const VarContainer::index<Index>::type& varsByIndex_;

public:
  UpdateVarData(const VarContainer::index<Index>::type& varsByIndex) : varsByIndex_(varsByIndex) {}
  virtual ~UpdateVarData() {}

  void operator()(Variable& var) const {
    var.clampValue = static_cast<double>(var.domSize) * static_cast<double>(var.adjList.size());
    double p2Cplx = var.domSize;
    BOOST_FOREACH( Var w, var.adjList ) {
      p2Cplx *= varsByIndex_[w].domSize;
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
  UpdateMinDegreeVarData(const VarContainer::index<Index>::type& varsByIndex) : UpdateVarData(varsByIndex) {}
};

class UpdateWeightedMinDegreeVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = var.clampValue;
  }
public:
  UpdateWeightedMinDegreeVarData(const VarContainer::index<Index>::type& varsByIndex) : UpdateVarData(varsByIndex) {}
};

class UpdateMinFillVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = 0.0;
    for (VarSet::const_iterator vAdjIter = var.adjList.begin(), vAdjEnd = var.adjList.end();
        vAdjIter != vAdjEnd; ++vAdjIter) {
      const Var u = *vAdjIter;
      const Variable& uVar = varsByIndex_[u];
      VarSet::const_iterator uAdjIter = uVar.adjList.upper_bound(u);
      VarSet::const_iterator uAdjEnd = uVar.adjList.end();
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
  UpdateMinFillVarData(const VarContainer::index<Index>::type& varsByIndex) : UpdateVarData(varsByIndex) {}
};

class UpdateWeightedMinFillVarData : public UpdateVarData {
private:
  virtual void updateCost(Variable& var) const {
    var.cost = 0.0;
    for (VarSet::const_iterator vAdjIter = var.adjList.begin(), vAdjEnd = var.adjList.end();
        vAdjIter != vAdjEnd; ++vAdjIter) {
      const Var u = *vAdjIter;
      const Variable& uVar = varsByIndex_[u];
      VarSet::const_iterator uAdjIter = uVar.adjList.upper_bound(u);
      VarSet::const_iterator uAdjEnd = uVar.adjList.end();
      VarSet::const_iterator vAdjIter2 = vAdjIter;
      ++vAdjIter2;
      double cost = 0.0;
      while (vAdjIter2 != vAdjEnd) {
        if (uAdjIter == uAdjEnd || *vAdjIter2 < *uAdjIter) {
          cost += varsByIndex_[*vAdjIter2].domSize;
          ++vAdjIter2;
        } else if (*uAdjIter < *vAdjIter2) {
          ++uAdjIter;
        } else {
          ++vAdjIter2;
          ++uAdjIter;
        }
      }
      var.cost += uVar.domSize * cost;
    }
  }
public:
  UpdateWeightedMinFillVarData(const VarContainer::index<Index>::type& varsByIndex) : UpdateVarData(varsByIndex) {}
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
  const VarContainer::index<Index>::type& varsByIndex_;
  virtual VarSet affectedVars(const Variable& var) const {
    VarSet vars = var.adjList;
    BOOST_FOREACH( Var u, var.adjList ) {
      const Variable& uVar = varsByIndex_[u];
      vars.insert(uVar.adjList.begin(), uVar.adjList.end());
    }
    vars.erase(var.index);
    return vars;
  }

public:
  MinFillAffectedVars(const VarContainer& varContainer) : varsByIndex_(varContainer.get<Index>()) {}
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
  using boost::scoped_ptr;
  using namespace greedyvarorder::internal;

  typedef VarContainer::index<Index>::type vars_by_index;
  typedef VarContainer::index<Cost>::type vars_by_cost;
  typedef VarContainer::index<Clamp>::type vars_by_clamp;
  typedef vars_by_cost::iterator cost_iterator;
  typedef vars_by_clamp::iterator clamp_iterator;

  if (task.numVars() != clampRank.size()) {
    throw InvalidArgumentException("clampRank size must equal the number of variables in task");
  }

  if (task.numVars() == 0) {
    return VarVector();
  }

  VarContainer vars(task, maxComplexity, clampRank);
  vars_by_index& varsByIndex = vars.get<Index>();
  vars_by_cost& varsByCost = vars.get<Cost>();
  vars_by_clamp& varsByClamp = vars.get<Clamp>();

  scoped_ptr<UpdateVarData> updateCostPtr;
  scoped_ptr<AffectedVars> affectedVarsPtr;
  switch (h) {
    case greedyvarorder::MIN_DEGREE:
      updateCostPtr.reset( new UpdateMinDegreeVarData(varsByIndex) );
      affectedVarsPtr.reset( new MinDegreeAffectedVars() );
      break;
    case greedyvarorder::WEIGHTED_MIN_DEGREE:
      updateCostPtr.reset( new UpdateWeightedMinDegreeVarData(varsByIndex) );
      affectedVarsPtr.reset( new MinDegreeAffectedVars() );
      break;
    case greedyvarorder::MIN_FILL:
      updateCostPtr.reset( new UpdateMinFillVarData(varsByIndex) );
      affectedVarsPtr.reset( new MinFillAffectedVars(vars) );
      break;
    case greedyvarorder::WEIGHTED_MIN_FILL:
      updateCostPtr.reset( new UpdateWeightedMinFillVarData(varsByIndex) );
      affectedVarsPtr.reset( new MinFillAffectedVars(vars) );
      break;
    default:
      throw InvalidArgumentException("Invalid heuristic");
  }

  for (vars_by_index::iterator iter = varsByIndex.begin(), end = varsByIndex.end(); iter != end; ++iter) {
    varsByIndex.modify<const UpdateVarData&>(iter, *updateCostPtr);
  }

  VarVector varOrder;
  int lastClampRank = -1;
  const Variable complexityUpper = Variable::complexityUpperBound(maxComplexity);

  for (;;) {

    cost_iterator minCostLower = varsByCost.begin();
    if (minCostLower->processed) {
      break;
    }

    if (minCostLower->complexity <= maxComplexity) {
      cost_iterator pickedIter = selectVar(minCostLower, varsByCost.upper_bound( Variable::upperBound(*minCostLower)),
          varsByCost.upper_bound(complexityUpper), rng, selectionScale);

      const Variable& v = *pickedIter;
      varOrder.push_back(v.index);
      VarSet affectedVars = (*affectedVarsPtr)(v);
      ElimNeighbour elimNeighbour(v.index, v.adjList);

      varsByCost.modify(pickedIter, MarkAsProcessed());

      BOOST_FOREACH( Var uIndex, v.adjList ) {
        varsByIndex.modify(varsByIndex.begin() + uIndex, elimNeighbour);
      }

      BOOST_FOREACH( Var uIndex, affectedVars ) {
        varsByIndex.modify<const UpdateVarData&>(varsByIndex.begin() + uIndex, *updateCostPtr);
      }

    } else {
      if (lastClampRank >= 0) {
        clamp_iterator clampIter = varsByClamp.upper_bound( Variable::clampRankUpperBound(lastClampRank) );
        clamp_iterator clampEnd = varsByClamp.end();
        while (clampIter != clampEnd && !clampIter->processed) {
          clamp_iterator here = clampIter++;
          varsByClamp.modify(here, DecrementClampRank());
        }
      }

      clamp_iterator clampLower = varsByClamp.begin();
      clamp_iterator pickedIter = selectVar(clampLower, varsByClamp.upper_bound( Variable::upperBound(*clampLower) ),
          varsByClamp.upper_bound( Variable::clampRankUpperBound(clampLower->clampRank) ), rng, selectionScale);

      const Variable& v = *pickedIter;
      lastClampRank = v.clampRank;
      varsByClamp.modify(pickedIter, MarkAsProcessed());
      ClampNeighbour clampNeighbour(v.index);

      BOOST_FOREACH( Var uIndex, v.adjList ) {
        varsByIndex.modify(varsByIndex.begin() + uIndex, clampNeighbour);
      }

      BOOST_FOREACH( Var uIndex, v.adjList ) {
        varsByIndex.modify<const UpdateVarData&>(varsByIndex.begin() + uIndex, *updateCostPtr);
      }
    }
  }

  return varOrder;
}

} // namespace orang

#endif
