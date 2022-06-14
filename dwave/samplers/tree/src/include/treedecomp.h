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
#ifndef INCLUDED_ORANG_TREEDECOMP_H
#define INCLUDED_ORANG_TREEDECOMP_H

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>
#include <sstream>
#include <cmath>
#include <memory>
#include <utility>

#include <base.h>
#include <exception.h>
#include <graph.h>

namespace orang {

class TreeDecompNode {
public:
  typedef std::vector<std::unique_ptr<TreeDecompNode>> node_vector;

private:
  TreeDecompNode* parent_;
  node_vector children_;
  Var nodeVar_;
  VarVector sepVars_;
  VarVector clampedVars_;

public:
  TreeDecompNode(Var nodeVar) :
    parent_(),
    children_(),
    nodeVar_(nodeVar),
    sepVars_(),
    clampedVars_() {}


  const TreeDecompNode* parent() const {
    return parent_;
  }

  TreeDecompNode*& parent() {
    return parent_;
  }

  const node_vector& children() const {
    return children_;
  }

  node_vector& children() {
    return children_;
  }

  Var nodeVar() const {
    return nodeVar_;
  }

  const VarVector& sepVars() const {
    return sepVars_;
  }

  VarVector& sepVars() {
    return sepVars_;
  }

  const VarVector& clampedVars() const {
    return clampedVars_;
  }

  VarVector& clampedVars() {
    return clampedVars_;
  }
};

class TreeDecomp {
private:
  std::size_t size_;
  std::size_t numVars_;
  double complexity_;
  VarVector clampedVars_;

  TreeDecompNode::node_vector roots_;

  void buildAdjSets(
      const Graph& g,
      const VarVector& varOrder,
      const VarVector& varRank,
      std::vector<VarSet>& adjSets,
      VarVector& parents) {

    // Copy partial graph structure in adjSets: for every edge {v,u} in g,
    // create an arc (v,u) iff both v and u appear in varOrder and u appears later than v.
    for (auto v: varOrder) {
      Var vRank = varRank[v];
      for (Graph::iterator adjIter = g.adjacencyBegin(v), adjEnd = g.adjacencyEnd(v); adjIter != adjEnd; ++adjIter) {
        Var u = *adjIter;
        Var uRank = varRank[u];
        if (uRank > vRank) {
          adjSets[v].insert(u);
        }
      }
    }

    // Now add extra arcs: fully connect (still respecting arc-rank relationship) all higher rank neighbours
    // of each node v in varOrder
    // Also record parent information: parent is lowest rank node u over all arcs (v,u)
    VarVector pRank(parents.size(), Var(-1));
    for (auto v: varOrder) {
      for (VarSet::const_iterator adjIter = adjSets[v].begin(), adjEnd = adjSets[v].end(); adjIter != adjEnd; ++adjIter) {
        Var u = *adjIter;
        Var uRank = varRank[u];

        if (uRank < pRank[v]) {
          pRank[v] = uRank;
          parents[v] = u;
        }

        VarSet::const_iterator adjIter2 = adjIter;
        for (++adjIter2; adjIter2 != adjEnd; ++adjIter2) {
          if (uRank > varRank[*adjIter2]) {
            adjSets[*adjIter2].insert(u);
          } else {
            adjSets[u].insert(*adjIter2);
          }
        }
      }
    }

  }

public:
  TreeDecomp(const Graph& g, const VarVector& varOrder, const DomIndexVector& domSizes) :
    size_(varOrder.size()),
    numVars_(g.numVertices()),
    complexity_(0.0),
    clampedVars_(),
    roots_() {

    using std::vector;
    using std::unique_ptr;
    using std::set_intersection;
    using std::back_inserter;
    using std::ostringstream;
    using std::max;
    using std::log;

    std::size_t orderSize = varOrder.size();

    if (domSizes.size() != numVars_) {
      ostringstream msg;
      msg << "Graph has " << numVars_ << " vertices but domSizes has " << domSizes.size() << " elements.";
      throw InvalidArgumentException(msg.str());
    }

    VarVector varRank(numVars_, 0);
    for (Var i = 0; i < orderSize; ++i) {
      Var voi = varOrder[i];

      if (voi >= numVars_) {
        ostringstream msg;
        msg << "Graph has " << numVars_ << " vertices (0-" << (numVars_ - 1)
          << ") but variable order contains " << voi;
        throw InvalidArgumentException(msg.str());
      }

      if (varRank[voi] != 0) {
        throw InvalidArgumentException("Variable order has duplicate entries");
      }

      if (domSizes[voi] == 0) {
        throw InvalidArgumentException("Domain size of zero given.");
      }

      varRank[varOrder[i]] = i + 1;
    }

    clampedVars_.reserve(numVars_ - orderSize);
    for (Var i = 0; i < numVars_; ++i) {
      if (varRank[i] == 0) {
        clampedVars_.push_back(i);
      }
    }

    vector<VarSet> adjSets(+numVars_);
    VarVector varParents(numVars_, Var(-1));
    buildAdjSets(g, varOrder, varRank, adjSets, varParents);

    vector<TreeDecompNode*> nodes(numVars_, 0);
    for (std::size_t i = orderSize; i-- > 0; ) {
      Var v = varOrder[i];
      Var vParent = varParents[v];

      unique_ptr<TreeDecompNode> vNode( new TreeDecompNode(v) );
      double nodeWidth = log(static_cast<double>(domSizes[v]));
      for (auto u: adjSets[v]) {
        vNode->sepVars().push_back(u);
        nodeWidth += log(static_cast<double>(domSizes[u]));
      }
      complexity_ = max(complexity_, nodeWidth);

      set_intersection(
          g.adjacencyBegin(v), g.adjacencyEnd(v),
          clampedVars_.begin(), clampedVars_.end(),
          back_inserter(vNode->clampedVars()));

      nodes[v] = vNode.get();
      if (vParent == Var(-1)) {
        roots_.push_back(std::move(vNode));
      } else {
        vNode->parent() = nodes[vParent];
        nodes[vParent]->children().push_back(std::move(vNode));
      }
    }

    static const double LOG2E = 1.4426950408889634074;
    complexity_ *= LOG2E;
  }

  std::size_t size() const {
    return size_;
  }

  std::size_t numVars() const {
    return numVars_;
  }

  double complexity() const {
    return complexity_;
  }

  const VarVector& clampedVars() const {
    return clampedVars_;
  }

  const TreeDecompNode::node_vector& roots() const {
    return roots_;
  }
};

} // end namespace orang

#endif
