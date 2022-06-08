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
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <iterator>
#include <ostream>
#include <set>
#include <utility>
#include <vector>
using std::inserter;
using std::copy;
using std::equal;
using std::make_pair;
using std::ostream;
using std::pair;
using std::set;
using std::transform;
using std::vector;

#include <boost/lambda/bind.hpp>
#include <boost/assign/list_of.hpp>
using boost::lambda::bind;
using boost::lambda::_1;
using boost::assign::list_of;

#include <base.h>
#include <graph.h>
using orang::Var;
using orang::VarVector;
using orang::Graph;

#include "test.h"

namespace {
namespace graphData {
const vector<Graph::adj_pair> adjList =
    list_of<Graph::adj_pair> (0,1) (2,0) (1,3) (2,4) (1,3) (3,4) (3,2) (1,0) (2,2);

const VarVector degrees = list_of(2)(2)(3)(3)(2);

const VarVector allAdjIters =
    list_of (1)(2)(999) (0)(3)(999) (0)(3)(4)(999) (1)(2)(4)(999) (2)(3)(999);
}
}

BOOST_AUTO_TEST_SUITE( graph )

BOOST_AUTO_TEST_CASE( graph_constructors )
{
  Graph emptyGraph;
  BOOST_CHECK_EQUAL(emptyGraph.numVertices(), 0);

  Graph graph1(graphData::adjList);
  BOOST_CHECK_EQUAL(graph1.numVertices(), graphData::degrees.size());
  VarVector degreeVector;
  VarVector allGraph1Adj;
  for (Var v = 0; v < graph1.numVertices(); ++v) {
    degreeVector.push_back(graph1.degree(v));
    for (Graph::iterator adjIter = graph1.adjacencyBegin(v), adjEnd = graph1.adjacencyEnd(v);
        adjIter != adjEnd; ++adjIter) {
      allGraph1Adj.push_back(*adjIter);
    }
    allGraph1Adj.push_back(999);
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(degreeVector.begin(), degreeVector.end(),
      graphData::degrees.begin(), graphData::degrees.end());
  BOOST_CHECK_EQUAL_COLLECTIONS(allGraph1Adj.begin(), allGraph1Adj.end(),
      graphData::allAdjIters.begin(), graphData::allAdjIters.end());


  Graph graph2(graph1);
  BOOST_CHECK_EQUAL(graph1, graph2);

  Graph graph3;
  graph3 = graph1;
  BOOST_CHECK_EQUAL(graph1, graph3);

  Graph graph4;
  set<pair<Var, Var> > adjSet;
  copy(graphData::adjList.begin(), graphData::adjList.end(), inserter(adjSet, adjSet.begin()));
  graph4.setAdjacencies(adjSet);
  BOOST_CHECK_EQUAL(graph1, graph4);

  Graph graph5(graphData::adjList, 100);
  BOOST_CHECK_EQUAL(graph5.numVertices(), 100);

  Graph graph6;
  graph6.setAdjacencies(graphData::adjList, 100);
  BOOST_CHECK_EQUAL(graph5, graph6);
}

BOOST_AUTO_TEST_SUITE_END()
