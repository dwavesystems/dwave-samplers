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
#include <boost/test/tools/floating_point_comparison.hpp>

#include <cmath>
#include <functional>
#include <utility>
#include <vector>
using std::mem_fn;
using std::pair;
using std::vector;

#include <boost/assign/list_of.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/lambda/bind.hpp>
using boost::assign::list_of;
using boost::make_transform_iterator;
using boost::lambda::bind;
using boost::lambda::_1;

#include <base.h>
#include <exception.h>
#include <graph.h>
#include <treedecomp.h>
using orang::Var;
using orang::VarVector;
using orang::DomIndexVector;
using orang::Graph;
using orang::TreeDecomp;
using orang::TreeDecompNode;
using orang::Exception;
using orang::InvalidArgumentException;

namespace {
namespace decompData {
const vector<Graph::adj_pair> adjList = list_of<Graph::adj_pair> (0,1) (0,4) (1,2) (1,5) (2,6) (3,4) (3,8)
    (4,5) (4,9) (5,6) (5,10) (6,7) (6,11) (7,12) (8,9) (8,13) (9,10) (9,14) (10,11) (10,15) (11,12) (11,16) (12,17)
    (13,14) (14,15) (14,18) (15,16) (15,19) (16,17) (16,20) (18,19) (19,20);
const Graph graph(adjList);

const VarVector varOrder1 = (list_of(0), 1, 2, 5, 6, 7, 11, 12, 17, 3, 8, 13, 9, 20, 19, 18, 15, 14);
const DomIndexVector domSizes1 = DomIndexVector(graph.numVertices(), 2);
const VarVector expectedPreorder1 = (list_of(14), 15, 10, 16, 18, 19, 20, 16, 9, 4, 10, 13, 8, 3, 4,
    17, 16, 12, 11, 10, 16, 7, 6, 5, 4, 10, 2, 1, 0, 4);
const VarVector expectedPostorder1 = (list_of(20), 19, 19, 15, 18, 18, 14, 15, 15, 14,
    3, 8, 8, 9, 13, 13, 9, 14, 9, 14, 14,
    0, 1, 1, 2, 5, 2, 5, 6, 5, 6, 6, 7, 11, 7, 11, 12, 11, 12, 12, 17, 17);
const double expectedComplexity1 = 3.0;
const VarVector expectedRoots1 = list_of(14)(17);
const VarVector expectedClamped1 = list_of(4)(10)(16);

const VarVector varOrder2 = (list_of(13), 18, 14, 15, 20, 16, 17, 11, 12, 7, 3, 0, 4, 1, 5, 2, 6);
const DomIndexVector domSizes2 = (list_of(2), 3, 2, 2, 4, 2, 2, 3, 100, 100, 100, 2, 2, 5, 2, 2, 3, 2, 2, 100, 4);
const VarVector expectedPreorder2 = (list_of(6), 2, 5, 10, 1, 4, 9, 0, 3, 8,
    7, 12, 11, 10, 17, 16, 20, 19, 15, 10, 19, 14, 9, 18, 19, 13, 8);
const VarVector expectedPostorder2 = (list_of(0), 1, 4, 3, 4, 4, 1, 5, 1, 2, 5, 5, 2, 6, 2, 6,
    20, 16, 18, 14, 13, 14, 14, 15, 15, 16, 16, 11, 17, 17, 11, 12, 11, 6, 12, 12, 6, 7, 7, 6,
    6);
const double expectedComplexity2 = std::log(24.0) / std::log(2.0);
const VarVector expectedRoots2 = list_of(6);
const VarVector expectedClamped2 = list_of(8)(9)(10)(19);

const VarVector badVarOrder1 = (list_of(0), 1, 2, 100);
const VarVector badVarOrder2 = (list_of(0), 1, 2, 3, 2);
const DomIndexVector shortDomSizes = DomIndexVector(3, 2);
const DomIndexVector zeroDomSizes = list_of(2).repeat(19,2)(0);

}
}

void preorderRecurse(const TreeDecompNode& n, VarVector& order) {
  order.push_back(n.nodeVar());
  for (auto cv: n.clampedVars()) {
    order.push_back(cv);
  }
  for (const auto &cn: n.children()) {
    preorderRecurse(*cn, order);
  }
}

VarVector preorderPlusClamped(const TreeDecomp& decomp) {
  VarVector order;
  for (const auto &r: decomp.roots()) {
    preorderRecurse(*r, order);
  }
  return order;
}

void postorderRecurse(const TreeDecompNode& n, VarVector& order) {
  for (const auto &cn: n.children()) {
    postorderRecurse(*cn, order);
  }
  order.push_back(n.nodeVar());
  for (auto sv: n.sepVars()) {
    order.push_back(sv);
  }
}

VarVector postorderPlusClamped(const TreeDecomp& decomp) {
  VarVector order;
  for (const auto &r: decomp.roots()) {
    postorderRecurse(*r, order);
  }
  return order;
}

BOOST_AUTO_TEST_SUITE( treedecomp )

BOOST_AUTO_TEST_CASE( treedecomp1 )
{
  try {
    TreeDecomp decomp(decompData::graph, decompData::varOrder1, decompData::domSizes1);
    BOOST_CHECK_EQUAL(decomp.numVars(), decompData::graph.numVertices());
    BOOST_CHECK_EQUAL(decomp.size(), decompData::varOrder1.size());
    BOOST_CHECK_CLOSE(decomp.complexity(), decompData::expectedComplexity1, 1e-4);
    BOOST_CHECK_EQUAL_COLLECTIONS(decomp.clampedVars().begin(), decomp.clampedVars().end(),
        decompData::expectedClamped1.begin(), decompData::expectedClamped1.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
        make_transform_iterator(decomp.roots().begin(), mem_fn(&TreeDecompNode::nodeVar)),
        make_transform_iterator(decomp.roots().end(), mem_fn(&TreeDecompNode::nodeVar)),
        decompData::expectedRoots1.begin(), decompData::expectedRoots1.end());

    VarVector preorder = preorderPlusClamped(decomp);
    VarVector postorder = postorderPlusClamped(decomp);
    BOOST_CHECK_EQUAL_COLLECTIONS(preorder.begin(), preorder.end(),
        decompData::expectedPreorder1.begin(), decompData::expectedPreorder1.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(postorder.begin(), postorder.end(),
        decompData::expectedPostorder1.begin(), decompData::expectedPostorder1.end());
  } catch (const Exception& e) {
    BOOST_ERROR(e.what());
  }
}

BOOST_AUTO_TEST_CASE( treedecomp2 )
{
  try {
    TreeDecomp decomp(decompData::graph, decompData::varOrder2, decompData::domSizes2);
    BOOST_CHECK_EQUAL(decomp.numVars(), decompData::graph.numVertices());
    BOOST_CHECK_EQUAL(decomp.size(), decompData::varOrder2.size());
    BOOST_CHECK_CLOSE(decomp.complexity(), decompData::expectedComplexity2, 1e-4);
    BOOST_CHECK_EQUAL_COLLECTIONS(decomp.clampedVars().begin(), decomp.clampedVars().end(),
        decompData::expectedClamped2.begin(), decompData::expectedClamped2.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(
        make_transform_iterator(decomp.roots().begin(), mem_fn(&TreeDecompNode::nodeVar)),
        make_transform_iterator(decomp.roots().end(), mem_fn(&TreeDecompNode::nodeVar)),
        decompData::expectedRoots2.begin(), decompData::expectedRoots2.end());

    VarVector preorder = preorderPlusClamped(decomp);
    VarVector postorder = postorderPlusClamped(decomp);
    BOOST_CHECK_EQUAL_COLLECTIONS(preorder.begin(), preorder.end(),
        decompData::expectedPreorder2.begin(), decompData::expectedPreorder2.end());
    BOOST_CHECK_EQUAL_COLLECTIONS(postorder.begin(), postorder.end(),
        decompData::expectedPostorder2.begin(), decompData::expectedPostorder2.end());
  } catch (const Exception& e) {
    BOOST_ERROR(e.what());
  }
}

BOOST_AUTO_TEST_CASE( treedecomp_exceptions )
{
  BOOST_CHECK_THROW(TreeDecomp decomp1(decompData::graph, decompData::badVarOrder1, decompData::domSizes1), InvalidArgumentException);
  BOOST_CHECK_THROW(TreeDecomp decomp2(decompData::graph, decompData::badVarOrder2, decompData::domSizes1), InvalidArgumentException);
  BOOST_CHECK_THROW(TreeDecomp decomp1(decompData::graph, decompData::varOrder1, decompData::shortDomSizes), InvalidArgumentException);
  BOOST_CHECK_THROW(TreeDecomp decomp1(decompData::graph, decompData::varOrder1, decompData::zeroDomSizes), InvalidArgumentException);
}


BOOST_AUTO_TEST_SUITE_END()
