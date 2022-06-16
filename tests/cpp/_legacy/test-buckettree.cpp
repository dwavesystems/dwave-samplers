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
#include <limits>
#include <ostream>
#include <set>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/iterator/indirect_iterator.hpp>

#include <base.h>
#include <exception.h>
#include <table.h>
#include <treedecomp.h>
#include <combine.h>
#include <operations/min.h>
#include <task.h>
#include <buckettree.h>

#include "test.h"

using std::copy;
using std::inserter;
using std::numeric_limits;
using std::ostream;
using std::ostream_iterator;
using std::set;
using std::vector;

using boost::assign::list_of;
using boost::assign::map_list_of;
using boost::make_indirect_iterator;

using orang::Var;
using orang::VarVector;
using orang::DomIndexVector;
using orang::OperationUnavailable;
using orang::Table;
using orang::TreeDecomp;
using orang::Plus;
using orang::MinSolutionSet;
using orang::MinOperations;
using orang::Task;
using orang::BucketTree;
using orang::NodeTables;

using tableAssign::vars;
using tableAssign::domSizes;
using tableAssign::values;
using tableAssign::none;

using minsolsetAssign::MinSolSetBuilder;

namespace {

template<typename Y>
struct NodeTableSet {
  typedef Table<Y> table_type;
  Var nodeVar;
  VarVector sepVars;
  set<table_type> tables;
  NodeTableSet() :
    nodeVar(numeric_limits<Var>::max()),
    sepVars(),
    tables() {}
  NodeTableSet(Var nodeVar0, const VarVector& sepVars0, const set<table_type>& tables0) :
    nodeVar(nodeVar0),
    sepVars(sepVars0),
    tables(tables0) {}
  NodeTableSet(const NodeTables<Y>& nt) :
    nodeVar(nt.nodeVar),
    sepVars(nt.sepVars),
    tables(make_indirect_iterator(nt.tables.begin()), make_indirect_iterator(nt.tables.end())) {}
};

template<typename Y>
bool operator==(const NodeTableSet<Y>& ns1, const NodeTableSet<Y>& ns2) {
  return ns1.nodeVar == ns2.nodeVar && ns1.sepVars == ns2.sepVars && ns1.tables == ns2.tables;
}

template<typename Y>
bool operator!=(const NodeTableSet<Y>& ns1, const NodeTableSet<Y>& ns2) {
  return !(ns1 == ns2);
}

template<typename Y>
bool operator<(const NodeTableSet<Y>& ns1, const NodeTableSet<Y>& ns2) {
  return ns1.nodeVar < ns2.nodeVar
      || (ns1.nodeVar == ns2.nodeVar
          && (ns1.sepVars < ns2.sepVars
              || (ns1.sepVars == ns2.sepVars && ns1.tables < ns2.tables)));
}

template<typename Y>
ostream& operator<<(ostream& out, const NodeTableSet<Y>& ns) {
  out << "NodeTablesSet(nodeVar: " << ns.nodeVar << ", sepVars: [";
  copy(ns.sepVars.begin(), ns.sepVars.end(), ostream_iterator<Var>(out, ","));
  out << "], tables(" << ns.tables.size() << "): [";
  copy(ns.tables.begin(), ns.tables.end(), ostream_iterator<Table<int> >(out, ", "));
  out << "])";
  return out;
}


namespace testData {

vector<Table<int> > tables= list_of<Table<int> >
((vars =  0,  1,  2, domSizes = 2, 2, 2, values =  6,  8, -7,  8,  3, -8, -4,  1))
((vars =  0,  1,  3, domSizes = 2, 2, 2, values =  9,  9, -7,  9,  9,  0,  6, -7))
((vars =  1,  2,  4, domSizes = 2, 2, 2, values = -1,  8,  6,  9,  3, -9,  7,  8))
((vars =  3,  4,     domSizes = 2, 2,    values =  3,  5,  5, -2))
((vars =  3,  5,     domSizes = 2, 3,    values =  3, -6,  4, -9, -4, -9))
((vars =  4,  6,     domSizes = 2, 3,    values = -8,  6,  4, -3,  9, -9))
((vars =  4,  7,     domSizes = 2, 3,    values = -1, -2,  5,  6, -6,  0))
((vars =  5,  8,     domSizes = 3, 2,    values = -1,  3,  4,  5, -4,  3))
((vars =  6,         domSizes = 3,       values =  3, -6, -7))
((vars =  6,  8,     domSizes = 3, 2,    values =  0,  9, -3,  2, -5,  5))
((vars =  7,  9,     domSizes = 3, 2,    values = -5,  0,  4,  7,  9,  1))
((vars =  8,  9,     domSizes = 2, 2,    values = -7, -7, -5,  6))
((vars =  8, 10, 11, domSizes = 2, 2, 2, values = -5,  6, -5,  8, -3, -6, -5,  2))
((vars =  9, 11, 12, domSizes = 2, 2, 2, values = -1, -3,  6,  2,  1,  8, -4,  5))
((vars = 10, 11, 12, domSizes = 2, 2, 2, values =  5, -2,  1, -8, -8,  1,  5,  8));

const auto tablesPtr = []{
  vector<Table<int>*> tmp;
  for (auto &table: tables) {
    tmp.push_back(&table);
  }
  return tmp;
}();

const VarVector varOrderAllClamped;
const VarVector varOrderNoClamped = list_of(0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10)(11)(12);
const VarVector varOrderTwoRoots = list_of(2)(1)(0)(5)(3)(10)(11)(12)(9)(7);

const DomIndexVector x0AllClamped = list_of(0)(1)(0)(1)(0)(2)(1)(2)(0)(1)(0)(1)(0);
const DomIndexVector x0NoClamped(13, 0);
const DomIndexVector x0TwoRoots = list_of(0)(0)(0)(0)(0)(0)(2)(0)(1)(0)(0)(0)(0);

const int expectedProblemValueAllClamped = 4;
const int expectedProblemValueNoClamped = -64;
const int expectedProblemValueTwoRoots = -21;

const MinSolutionSet<int> expectedSolutionAllClamped = MinSolSetBuilder<int>(1) (0, x0AllClamped);
const MinSolutionSet<int> expectedSolutionNoClamped = MinSolSetBuilder<int>(1)
    (0, list_of(0)(1)(0)(1)(1)(0)(2)(0)(0)(0)(0)(0)(1));
const MinSolutionSet<int> expectedSolutionTwoRoots = MinSolSetBuilder<int>(1)
    (0, list_of(1)(0)(1)(1)(0)(1)(2)(0)(1)(0)(0)(1)(1));

const set<NodeTableSet<int> > expectedNodeTablesAllClamped;
const set<NodeTableSet<int> > expectedNodeTablesNoClamped = list_of<NodeTableSet<int> >
(NodeTableSet<int>( 0, list_of(1)(2)(3),    list_of<Table<int> >
  // base tables
  (tables[0])
  (tables[1])
  // pi table
  ((vars = 1, 2, 3, domSizes = 2, 2, 2, values = -37, -49, -33, -32, -51, -63, -47, -46))
))
(NodeTableSet<int>( 1, list_of(2)(3)(4),    list_of<Table<int> >
  // base tables
  (tables[2])
  // lambda tables
  ((vars = 1, 2, 3, domSizes = 2, 2, 2, values = 15, -14, 1, -11, 8, -1, -8, -6))
  // pi table
  ((vars = 2, 3, 4, domSizes = 2, 2, 2, values = -27, -27, -33, -33, -40, -40, -54, -54))
))
(NodeTableSet<int>( 2, list_of(3)(4),       list_of<Table<int> >
  // lambda tables
  ((vars = 2, 3, 4, domSizes = 2, 2, 2, values = -6, -2, 7, -2, -23, -3, -10, -1))
  // pi table
  ((vars = 3, 4, domSizes = 2, 2, values = -27, -33, -40, -54))
))
(NodeTableSet<int>( 3, list_of(4)(5),       list_of<Table<int> >
  // base tables
  (tables[3])
  (tables[4])
  // lambda tables
  ((vars = 3, 4, domSizes = 2, 2, values = -6, -2, -23, -10))
  // pi table
  ((vars = 4, 5, domSizes = 2, 3, values = -31, -46, -29, -42, -26, -41))
))
(NodeTableSet<int>( 4, list_of(5)(6)(7),    list_of<Table<int> >
  // base tables
  (tables[5])
  (tables[6])
  // lambda tables
  ((vars = 4, 5, domSizes = 2, 3, values = -3, -18, -6, -21, -7, -22))
  // pi table
  ((vars = 5, 6, 7, domSizes = 3, 3, 3, values = -22, -18, -17, -23, -32, -25, -35, -31, -30,
      -17, -13, -12, -18, -27, -20, -30, -26, -25, -13, -9, -8, -14, -23, -16, -26, -22, -21))
))
(NodeTableSet<int>( 5, list_of(6)(7)(8),    list_of<Table<int> >
  // base tables
  (tables[7])
  // lambda tables
  ((vars = 5, 6, 7, domSizes = 3, 3, 3, values = -14, -17, -18, -23, -26, -27, -29, -32, -33,
      -6, -9, -10, -15, -18, -19, -21, -24, -25, -17, -20, -21, -21, -24, -25, -27, -30, -31))
  // pi table
  ((vars = 6, 7, 8, domSizes = 3, 3, 2, values = -21, -21, -34, -16, -16, -29, -12, -12, -25,
      -12, -28, -19, -7, -23, -14, -3, -19, -10))
))
(NodeTableSet<int>( 6, list_of(7)(8),       list_of<Table<int> >
  // base tables
  (tables[8])
  (tables[9])
  // lambda tables
  ((vars = 6, 7, 8, domSizes = 3, 3, 2, values = -15, -24, -30, -7, -16, -22, -18, -22, -28,
      -21, -30, -36, -13, -22, -28, -24, -28, -34))
  // pi table
  ((vars = 7, 8, domSizes = 3, 2, values = -24, -19, -15, -17, -12, -8))
))
(NodeTableSet<int>( 7, list_of(8)(9),       list_of<Table<int> >
  // base tables
  (tables[10])
  // lambda tables
  ((vars = 7, 8, domSizes = 3, 2, values = -40, -32, -38, -41, -33, -39))
  // pi table
  ((vars = 8, 9, domSizes = 2, 2, values = -19, -12, -16, 2))
))
(NodeTableSet<int>( 8, list_of(9)(10)(11),  list_of<Table<int> >
  // base tables
  (tables[11])
  (tables[12])
  // lambda tables
  ((vars = 8, 9, domSizes = 2, 2, values = -45, -46, -37, -38))
  // pi table
  ((vars = 9, 10, 11, domSizes = 2, 2, 2, values = -7, 0, -3, -5, 1, 3, -2, -6))
))
(NodeTableSet<int>( 9, list_of(10)(11)(12), list_of<Table<int> >
// base tables
  (tables[13])
  // lambda tables
  ((vars = 9, 10, 11, domSizes = 2, 2, 2, values = -57, -47, -57, -47, -59, -45, -57, -47))
  // pi table
  ((vars = 10, 11, 12, domSizes = 2, 2, 2, values = 5, -2, 1, -8, -8, 1, 5, 8))
))
(NodeTableSet<int>(10, list_of(11)(12),     list_of<Table<int> >
// base tables
  (tables[14])
  // lambda tables
  ((vars = 10, 11, 12, domSizes = 2, 2, 2, values = -58, -58, -53, -51, -56, -56, -63, -61))
))
(NodeTableSet<int>(11, list_of(12),         list_of<Table<int> >
// lambda tables
  ((vars = 11, 12, domSizes = 2, 2, values = -60, -59, -64, -58))
))
(NodeTableSet<int>(12, VarVector(),         list_of<Table<int> >
// lambda tables
  ((vars = 12, domSizes = 2, values = -60, -64))
));

const set<NodeTableSet<int> > expectedNodeTablesTwoRoots = list_of<NodeTableSet<int> >
(NodeTableSet<int>( 0, list_of(3),       list_of<Table<int> >
  // lambda tables
  ((vars = 0, 3, domSizes = 2, 2, values = -6, 7, 7, -2))
  // pi table
  ((vars = 3, domSizes = 2, values = 2, -8))
))
(NodeTableSet<int>( 1, list_of(0)(3),    list_of<Table<int> >
  // base tables
  (tables[1])
  //lambda tables
  ((vars = 0, 1, domSizes = 2, 2, values = 5, -2, 1, 10))
  // pi table
  ((vars = 0, 3, domSizes = 2, 2, values = 2, 2, -8, -8))
))
(NodeTableSet<int>( 2, list_of(0)(1),    list_of<Table<int> >
  // base tables
  (tables[0])
  ((vars = 1, 2, domSizes = 2, 2, values = -1, 8, 6, 9)) // tables[2] with x4=0
  // pi table
  ((vars = 0, 1, domSizes = 2, 2, values = 1, -8, -5, -15))
))
(NodeTableSet<int>( 3, VarVector(),      list_of<Table<int> >
  // base tables
  ((vars = 3, domSizes = 2, values = 3, 5)) // tables[3] with x4=0
  //lambda tables
  ((vars = 3, domSizes = 2, values = -6, -2))
  ((vars = 3, domSizes = 2, values = -1, -13))
))
(NodeTableSet<int>( 5, list_of(3),       list_of<Table<int> >
  // base tables
  (tables[4])
  ((vars = 5, domSizes = 3, values = 5, -4, 3)) // tables[7] with x8=1
  // pi table
  ((vars = 3, domSizes = 2, values = -3, 3))
))
(NodeTableSet<int>( 7, VarVector(),      list_of<Table<int> >
  // base tables
  ((vars = 7, domSizes = 3, values = -1, 5, -6)) // tables[6] with x4=0
  //lambda tables
  ((vars = 7, domSizes = 3, values = -17, -12, -8))
))
(NodeTableSet<int>( 9, list_of(7),       list_of<Table<int> >
  // base tables
  (tables[10])
  ((vars = 9, domSizes = 2, values = -7, 6)) // tables[11] with x8=1
  //lambda tables
  ((vars = 9, domSizes = 2, values = -5, -4))
  // pi table
  ((vars = 7, domSizes = 3, values = -1, 5, -6))
))
(NodeTableSet<int>( 10, list_of(11)(12), list_of<Table<int> >
  // base tables
  (tables[14])
  ((vars = 10, 11, domSizes = 2, 2, values = 6, 8, -6, 2)) // tables[12] with x8=1
  // pi table
  ((vars = 11, 12, domSizes = 2, 2, values = -14, -7, -12, -17))
))
(NodeTableSet<int>( 11, list_of(9)(12),  list_of<Table<int> >
  // base tables
  (tables[13])
  //lambda tables
  ((vars = 11, 12, domSizes = 2, 2, values = 6, -6, -2, -1))
  // pi table
  ((vars = 9, 12, domSizes = 2, 2, values = -13, 1, -13, 1))
))
(NodeTableSet<int>( 12, list_of(9),      list_of<Table<int> >
  //lambda tables
  ((vars = 9, 12, domSizes = 2, 2, values = 0, -4, -5, 4))
  // pi table
  ((vars = 9, domSizes = 2, values = -13, 1))
))
;


} // namespace (anonymous)::testData


typedef MinOperations<int, Plus<int> > ops_type;
typedef ops_type::CtorArgs ops_ctorargs;
typedef Task<ops_type> task_type;

} // anonymous namespace

BOOST_AUTO_TEST_SUITE( bucket_tree )


BOOST_AUTO_TEST_CASE( allclamped_nosolve_notables )
{
  const VarVector& varOrder = testData::varOrderAllClamped;
  const DomIndexVector& x0 = testData::x0AllClamped;
  const int expectedProblemValue = testData::expectedProblemValueAllClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, false, false);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  BOOST_CHECK_THROW(bucketTree.solve(), OperationUnavailable);
  BOOST_CHECK_THROW(bucketTree.nodeTables(), OperationUnavailable);
}

BOOST_AUTO_TEST_CASE( noclamped_nosolve_notables )
{
  const VarVector& varOrder = testData::varOrderNoClamped;
  const DomIndexVector& x0 = testData::x0NoClamped;
  const int expectedProblemValue = testData::expectedProblemValueNoClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, false, false);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  BOOST_CHECK_THROW(bucketTree.solve(), OperationUnavailable);
  BOOST_CHECK_THROW(bucketTree.nodeTables(), OperationUnavailable);
}

BOOST_AUTO_TEST_CASE( tworoots_nosolve_notables )
{
  const VarVector& varOrder = testData::varOrderTwoRoots;
  const DomIndexVector& x0 = testData::x0TwoRoots;
  const int expectedProblemValue = testData::expectedProblemValueTwoRoots;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, false, false);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  BOOST_CHECK_THROW(bucketTree.solve(), OperationUnavailable);
  BOOST_CHECK_THROW(bucketTree.nodeTables(), OperationUnavailable);
}


BOOST_AUTO_TEST_CASE( allclamped_solve_notables )
{
  const VarVector& varOrder = testData::varOrderAllClamped;
  const DomIndexVector& x0 = testData::x0AllClamped;
  const int expectedProblemValue = testData::expectedProblemValueAllClamped;
  const MinSolutionSet<int>& expectedSolution = testData::expectedSolutionAllClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 2);
  task.maxSolutions(1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, true, false);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  MinSolutionSet<int> solution = bucketTree.solve();
  BOOST_CHECK_EQUAL_COLLECTIONS(solution.solutions().begin(), solution.solutions().end(),
      expectedSolution.solutions().begin(), expectedSolution.solutions().end());
  BOOST_CHECK_THROW(bucketTree.nodeTables(), OperationUnavailable);
}

BOOST_AUTO_TEST_CASE( noclamped_solve_notables )
{
  const VarVector& varOrder = testData::varOrderNoClamped;
  const DomIndexVector& x0 = testData::x0NoClamped;
  const int expectedProblemValue = testData::expectedProblemValueNoClamped;
  const MinSolutionSet<int>& expectedSolution = testData::expectedSolutionNoClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 2);
  task.maxSolutions(1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, true, false);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  MinSolutionSet<int> solution = bucketTree.solve();
  BOOST_CHECK_EQUAL_COLLECTIONS(solution.solutions().begin(), solution.solutions().end(),
      expectedSolution.solutions().begin(), expectedSolution.solutions().end());
  BOOST_CHECK_THROW(bucketTree.nodeTables(), OperationUnavailable);
}

BOOST_AUTO_TEST_CASE( tworoots_solve_notables )
{
  const VarVector& varOrder = testData::varOrderTwoRoots;
  const DomIndexVector& x0 = testData::x0TwoRoots;
  const int expectedProblemValue = testData::expectedProblemValueTwoRoots;
  const MinSolutionSet<int>& expectedSolution = testData::expectedSolutionTwoRoots;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 2);
  task.maxSolutions(1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, true, false);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  MinSolutionSet<int> solution = bucketTree.solve();
  BOOST_CHECK_EQUAL_COLLECTIONS(solution.solutions().begin(), solution.solutions().end(),
      expectedSolution.solutions().begin(), expectedSolution.solutions().end());
  BOOST_CHECK_THROW(bucketTree.nodeTables(), OperationUnavailable);
}


BOOST_AUTO_TEST_CASE( allclamped_nosolve_tables )
{
  const VarVector& varOrder = testData::varOrderAllClamped;
  const DomIndexVector& x0 = testData::x0AllClamped;
  const int expectedProblemValue = testData::expectedProblemValueAllClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, false, true);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  BOOST_CHECK_THROW(bucketTree.solve(), OperationUnavailable);
  BOOST_CHECK(bucketTree.nodeTables().empty());
}

BOOST_AUTO_TEST_CASE( noclamped_nosolve_tables )
{
  const VarVector& varOrder = testData::varOrderNoClamped;
  const DomIndexVector& x0 = testData::x0NoClamped;
  const int expectedProblemValue = testData::expectedProblemValueNoClamped;
  const set<NodeTableSet<int> >& expectedNodeTablesSet = testData::expectedNodeTablesNoClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, false, true);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  BOOST_CHECK_THROW(bucketTree.solve(), OperationUnavailable);

  set<NodeTableSet<int> > nodeTablesSet;
  const vector<NodeTables<int> >& nodeTablesVector = bucketTree.nodeTables();
  copy(nodeTablesVector.begin(), nodeTablesVector.end(), inserter(nodeTablesSet, nodeTablesSet.begin()));
  BOOST_CHECK_EQUAL_COLLECTIONS(nodeTablesSet.begin(), nodeTablesSet.end(),
      expectedNodeTablesSet.begin(), expectedNodeTablesSet.end());
}

BOOST_AUTO_TEST_CASE( tworoots_nosolve_tables )
{
  const VarVector& varOrder = testData::varOrderTwoRoots;
  const DomIndexVector& x0 = testData::x0TwoRoots;
  const int expectedProblemValue = testData::expectedProblemValueTwoRoots;
  const set<NodeTableSet<int> >& expectedNodeTablesSet = testData::expectedNodeTablesTwoRoots;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, false, true);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  BOOST_CHECK_THROW(bucketTree.solve(), OperationUnavailable);

  set<NodeTableSet<int> > nodeTablesSet;
  const vector<NodeTables<int> >& nodeTablesVector = bucketTree.nodeTables();
  copy(nodeTablesVector.begin(), nodeTablesVector.end(), inserter(nodeTablesSet, nodeTablesSet.begin()));
  BOOST_CHECK_EQUAL_COLLECTIONS(nodeTablesSet.begin(), nodeTablesSet.end(),
      expectedNodeTablesSet.begin(), expectedNodeTablesSet.end());
}


BOOST_AUTO_TEST_CASE( allclamped_solve_tables )
{
  const VarVector& varOrder = testData::varOrderAllClamped;
  const DomIndexVector& x0 = testData::x0AllClamped;
  const MinSolutionSet<int>& expectedSolution = testData::expectedSolutionAllClamped;
  const int expectedProblemValue = testData::expectedProblemValueAllClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  task.maxSolutions(1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, true, true);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);
  MinSolutionSet<int> solution = bucketTree.solve();
  BOOST_CHECK_EQUAL_COLLECTIONS(solution.solutions().begin(), solution.solutions().end(),
      expectedSolution.solutions().begin(), expectedSolution.solutions().end());
  BOOST_CHECK(bucketTree.nodeTables().empty());
}

BOOST_AUTO_TEST_CASE( noclamped_solve_tables )
{
  const VarVector& varOrder = testData::varOrderNoClamped;
  const DomIndexVector& x0 = testData::x0NoClamped;
  const int expectedProblemValue = testData::expectedProblemValueNoClamped;
  const MinSolutionSet<int>& expectedSolution = testData::expectedSolutionNoClamped;
  const set<NodeTableSet<int> >& expectedNodeTablesSet = testData::expectedNodeTablesNoClamped;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  task.maxSolutions(1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, true, true);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);

  MinSolutionSet<int> solution = bucketTree.solve();
  BOOST_CHECK_EQUAL_COLLECTIONS(solution.solutions().begin(), solution.solutions().end(),
      expectedSolution.solutions().begin(), expectedSolution.solutions().end());

  set<NodeTableSet<int> > nodeTablesSet;
  const vector<NodeTables<int> >& nodeTablesVector = bucketTree.nodeTables();
  copy(nodeTablesVector.begin(), nodeTablesVector.end(), inserter(nodeTablesSet, nodeTablesSet.begin()));
  BOOST_CHECK_EQUAL_COLLECTIONS(nodeTablesSet.begin(), nodeTablesSet.end(),
      expectedNodeTablesSet.begin(), expectedNodeTablesSet.end());
}

BOOST_AUTO_TEST_CASE( tworoots_solve_tables )
{
  const VarVector& varOrder = testData::varOrderTwoRoots;
  const DomIndexVector& x0 = testData::x0TwoRoots;
  const int expectedProblemValue = testData::expectedProblemValueTwoRoots;
  const MinSolutionSet<int>& expectedSolution = testData::expectedSolutionTwoRoots;
  const set<NodeTableSet<int> >& expectedNodeTablesSet = testData::expectedNodeTablesTwoRoots;

  task_type task(testData::tablesPtr.begin(), testData::tablesPtr.end(), 1);
  task.maxSolutions(1);
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  BucketTree<task_type> bucketTree(task, decomp, x0, true, true);

  BOOST_CHECK_EQUAL(bucketTree.problemValue(), expectedProblemValue);

  MinSolutionSet<int> solution = bucketTree.solve();
  BOOST_CHECK_EQUAL_COLLECTIONS(solution.solutions().begin(), solution.solutions().end(),
      expectedSolution.solutions().begin(), expectedSolution.solutions().end());

  set<NodeTableSet<int> > nodeTablesSet;
  const vector<NodeTables<int> >& nodeTablesVector = bucketTree.nodeTables();
  copy(nodeTablesVector.begin(), nodeTablesVector.end(), inserter(nodeTablesSet, nodeTablesSet.begin()));
  BOOST_CHECK_EQUAL_COLLECTIONS(nodeTablesSet.begin(), nodeTablesSet.end(),
      expectedNodeTablesSet.begin(), expectedNodeTablesSet.end());
}

BOOST_AUTO_TEST_SUITE_END()
