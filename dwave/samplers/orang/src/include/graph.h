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
#ifndef INCLUDED_ORANG_GRAPH_H
#define INCLUDED_ORANG_GRAPH_H

#include <algorithm>
#include <set>
#include <utility>

#include <base.h>

namespace orang {

class Graph {
private:
  SizeVector adjOffsets_;
  VarVector adj_;

public:
  typedef VarVector::const_iterator iterator;
  typedef std::pair<Var, Var> adj_pair;

  Graph() : adjOffsets_(1), adj_() {}

  template<typename Adj>
  explicit Graph(const Adj& adjSet, Var minVars = 0) : adjOffsets_(), adj_() {
    setAdjacencies(adjSet, minVars);
  }

  template<typename Adj>
  void setAdjacencies(const Adj& adjSet, Var minVars = 0) {
    typedef std::set<adj_pair> adjacency_set;

    adjOffsets_.clear();
    adj_.clear();

    Var numVars = minVars;
    adjacency_set symAdjSet;
    for (const auto &adjPair: adjSet) {
      if (adjPair.first != adjPair.second) {
        symAdjSet.insert(adjPair);
        symAdjSet.insert(std::make_pair(adjPair.second, adjPair.first));
      }
      numVars = std::max(numVars, 1 + std::max(adjPair.first, adjPair.second));
    }

    adjOffsets_.reserve(numVars + 1);
    adj_.reserve(symAdjSet.size());

    Var lastI = 0;
    adjOffsets_.push_back(lastI);

    for (const auto &adjPair: symAdjSet) {
      while (lastI <= adjPair.first) {
        ++lastI;
        adjOffsets_.push_back(adjOffsets_.back());
      }

      ++adjOffsets_.back();
      adj_.push_back(adjPair.second);
    }

    adjOffsets_.resize(numVars + 1, adjOffsets_.back());
  }

  Var numVertices() const {
    return static_cast<Var>(adjOffsets_.size() - 1);
  }

  Var degree(Var v) const {
    return static_cast<Var>(adjOffsets_.at(v + 1) - adjOffsets_.at(v));
  }

  iterator adjacencyBegin(Var v) const {
    return adj_.begin() + adjOffsets_.at(v);
  }

  iterator adjacencyEnd(Var v) const {
    return adj_.begin() + adjOffsets_.at(v + 1);
  }
};

} // namespace orang

#endif
