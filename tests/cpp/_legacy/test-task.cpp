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

#include <cstddef>
#include <exception>
#include <set>
#include <sstream>
#include <vector>

#include <boost/assign/list_of.hpp>
#include <boost/iterator/indirect_iterator.hpp>

#include <base.h>
#include <exception.h>
#include <table.h>
#include <combine.h>
#include <marginalizer.h>
#include <operations/min.h>
#include <operations/logsumprod.h>
#include <task.h>

#include "test.h"

using std::size_t;
using std::set;
using std::vector;
using std::ostringstream;

using boost::assign::list_of;
using boost::assign::pair_list_of;
using boost::make_indirect_iterator;
using boost::indirect_iterator;

using orang::Var;
using orang::VarVector;
using orang::DomIndexVector;
using orang::SizeVector;
using orang::Exception;
using orang::Table;
using orang::TableVar;
using orang::Graph;
using orang::Task;
using orang::Plus;
using orang::MinOperations;
using orang::LogSumProductOperations;
using orang::InvalidArgumentException;
using orang::TreeDecompNode;

using tableAssign::vars;
using tableAssign::domSizes;
using tableAssign::values;
using tableAssign::none;

namespace {

class TableVarIndexIter {
private:
  vector<TableVar>::const_iterator iter_;
public:
  TableVarIndexIter(const vector<TableVar>::const_iterator& iter) : iter_(iter) {}
  TableVarIndexIter(const TableVarIndexIter& other) : iter_(other.iter_) {}
  Var operator*() { return iter_->index; }
  TableVarIndexIter& operator++() { ++iter_; return *this; }
  TableVarIndexIter operator++(int) { TableVarIndexIter r(iter_); ++iter_; return r; }
};

class TableVarDomSizeIter {
private:
  vector<TableVar>::const_iterator iter_;
public:
  TableVarDomSizeIter(const vector<TableVar>::const_iterator& iter) : iter_(iter) {}
  TableVarDomSizeIter(const TableVarDomSizeIter& other) : iter_(other.iter_) {}
  size_t operator*() { return iter_->domSize; }
  TableVarDomSizeIter& operator++() { ++iter_; return *this; }
  TableVarDomSizeIter operator++(int) { TableVarDomSizeIter r(iter_); ++iter_; return r; }
};

TreeDecompNode initDNode(Var nodeVar, const VarVector& sepVars, const VarVector& clampedVars) {
  TreeDecompNode dNode(nodeVar);
  dNode.sepVars() = sepVars;
  dNode.clampedVars() = clampedVars;
  return dNode;
}

struct DummyRng {
  double operator()() { return 0.0; }
};

namespace testData {

vector<Table<int> > goodTables = list_of<Table<int> >
((vars = none,    domSizes = none,    values = 9999))
((vars = 0,       domSizes = 2,       values = -1, 1))
((vars = 5,       domSizes = 2,       values = 1, 10))
((vars = 0, 1,    domSizes = 2, 2,    values = 0, 1, 2, -4 ))
((vars = 4, 5,    domSizes = 3, 2,    values = -1, -1, -2, -3, -5, -8 ))
((vars = 0, 1, 2, domSizes = 2, 2, 4, values = 2, 7, 1, 8, 2, 8, 1, 8, 2, 8, 4, 5, 9, 0, 4, 5))
((vars = 1, 4, 5, domSizes = 2, 3, 2, values = 3, 1, 4, 1, 5, 9, 2, 6, 5, 3, 5, 6))
((vars = 6,       domSizes = 5,       values = 0, 0, 1, 0, 0));

const auto goodTablesPtr = []{
  vector<Table<int>*> tmp;
  for (auto &table: goodTables) {
    tmp.push_back(&table);
  }
  return tmp;
}();


const SizeVector expectedDomSizes = list_of(2)(2)(4)(1)(3)(2)(5)(1)(1);

const Graph expectedGraph(list_of<Graph::adj_pair> (0,1) (0,2) (1,2) (1,4) (1,5) (4,5)
    .convert_to_container<vector<Graph::adj_pair> >(), 9);

vector<Table<int> > badTables = list_of<Table<int> >
((vars = 1, 2, 3, domSizes = 2, 2, 2, values = 0, 0, 0, 0, 0, 0, 0, 0))
((vars = 0, 3,    domSizes = 2, 3,    values = 1, 1, 1, 1, 1, 1));

const auto badTablesPtr = []{
  vector<Table<int>*> tmp;
  for (auto &table: badTables) {
    tmp.push_back(&table);
  }
  return tmp;
}();

const TreeDecompNode dNode = initDNode( 1, (list_of(2)(5)), (list_of(4)) );
const DomIndexVector baseTablesX0 = list_of(0)(0)(0)(0)(1)(0)(0);

const set<Table<int> > expectedTables = list_of<Table<int> >
((vars = 1, 5, domSizes = 2, 2, values = 4, 1, 5, 3));

vector<Table<double> > doubleTables = list_of<Table<double> >
((vars = 0, domSizes = 2, values = 1.0, 2.0));

const auto doubleTablesPtr = []{
  vector<Table<double>*> tmp;
  for (auto &table: doubleTables) {
    tmp.push_back(&table);
  }
  return tmp;
}();


const vector<double> rootValues = list_of(1.0)(2.0)(3.0);

const VarVector clampedVars = list_of(1)(2)(4)(5);

const DomIndexVector problemValueX0 = list_of(0)(1)(0)(0)(2)(1)(0);

const double expectedProblemValue = 10013.0;

DummyRng rng;

} // namespace (anonymous)::testData
} // anonymous namespace

BOOST_AUTO_TEST_SUITE( task )

typedef MinOperations<int, Plus<int> > min_ops_type;
typedef Task<min_ops_type> min_task_type;

typedef LogSumProductOperations<DummyRng> logsumprod_ops_type;
typedef logsumprod_ops_type::CtorArgs logsumprod_ops_ctorargs;
typedef Task<logsumprod_ops_type> logsumprod_task_type;

BOOST_AUTO_TEST_CASE( constructor )
{
  min_task_type task(testData::goodTablesPtr.begin(), testData::goodTablesPtr.end(), 1,
    static_cast<Var>(testData::expectedDomSizes.size()));

  BOOST_CHECK_EQUAL(task.numVars(), testData::expectedDomSizes.size());

  BOOST_CHECK_EQUAL_COLLECTIONS(task.domSizes().begin(), task.domSizes().end(),
      testData::expectedDomSizes.begin(), testData::expectedDomSizes.end());

  SizeVector singleDomSizes;
  for (Var i = 0; i < testData::expectedDomSizes.size(); ++i) {
    singleDomSizes.push_back(task.domSize(i));
  }
  BOOST_CHECK_EQUAL_COLLECTIONS(singleDomSizes.begin(), singleDomSizes.end(),
      testData::expectedDomSizes.begin(), testData::expectedDomSizes.end());

  BOOST_CHECK_EQUAL(task.graph(), testData::expectedGraph);
}

BOOST_AUTO_TEST_CASE( constructor_exception )
{
  BOOST_CHECK_THROW(
      min_task_type task(testData::badTablesPtr.begin(), testData::badTablesPtr.end(), 1),
      InvalidArgumentException);
}

BOOST_AUTO_TEST_CASE( tables )
{
  min_task_type task(testData::goodTablesPtr.begin(), testData::goodTablesPtr.end(), 1);
  BOOST_CHECK_EQUAL_COLLECTIONS(make_indirect_iterator(
      task.tables().begin()), make_indirect_iterator(task.tables().end()),
      testData::goodTables.begin(), testData::goodTables.end());
}

BOOST_AUTO_TEST_CASE( base_tables )
{
  try {
    min_task_type task(testData::goodTablesPtr.begin(), testData::goodTablesPtr.end(), 1);

    min_task_type::table_vector tables = task.baseTables(testData::dNode, testData::baseTablesX0);
    BOOST_CHECK_EQUAL_COLLECTIONS(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()),
        testData::expectedTables.begin(), testData::expectedTables.end());
  } catch (Exception& e) {
    ostringstream msg;
    msg << "Caught orang exception: " << e.what();
    BOOST_FAIL(msg.str());
  }

}

BOOST_AUTO_TEST_CASE( problem_value )
{
  logsumprod_task_type task(testData::goodTablesPtr.begin(), testData::goodTablesPtr.end(), testData::rng);

  double pv = task.problemValue(testData::rootValues, testData::problemValueX0, testData::clampedVars);
  BOOST_CHECK_EQUAL(pv, testData::expectedProblemValue);
}

BOOST_AUTO_TEST_CASE ( task_as_ops )
{
  min_task_type task(testData::doubleTablesPtr.begin(), testData::doubleTablesPtr.end(), 1);
  task.maxSolutions(2);
  BOOST_CHECK_EQUAL(task.maxSolutions(), 2);
}

BOOST_AUTO_TEST_SUITE_END()
