#ifndef INCLUDED_ORANG_GRAPH_H
#define INCLUDED_ORANG_GRAPH_H

#include <algorithm>
#include <set>
#include <utility>

#include <boost/foreach.hpp>

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
    BOOST_FOREACH( const adj_pair& adjPair, adjSet ) {
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

    BOOST_FOREACH( const adj_pair& adjPair, symAdjSet ) {
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
