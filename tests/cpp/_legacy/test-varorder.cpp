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
#include <vector>
#include <limits>

#include <algorithm>
#include <iterator>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/assign/list_of.hpp>

#include <base.h>
#include <table.h>
#include <operations/dummy.h>
#include <task.h>
#include <varorder.h>

#include "test.h"

using std::vector;
using std::numeric_limits;

using boost::assign::list_of;

using orang::VarVector;
using orang::Table;
using orang::DummyOperations;
using orang::Task;
using orang::greedyVarOrder;

using tableAssign::vars;
using tableAssign::domSizes;

namespace {
namespace testData {

vector<Table<int> > tables = list_of<Table<int> >
((vars =  0,  9, domSizes =  2,  3))
((vars =  0, 18, domSizes =  2,  2))
((vars =  0, 28, domSizes =  2,  3))
((vars =  1, 12, domSizes =  2,  3))
((vars =  1, 17, domSizes =  2,  2))
((vars =  1, 20, domSizes =  2,  2))
((vars =  1, 21, domSizes =  2,  1))
((vars =  1, 23, domSizes =  2,  2))
((vars =  2, 15, domSizes =  2,  2))
((vars =  2, 17, domSizes =  2,  2))
((vars =  2, 24, domSizes =  2,  2))
((vars =  2, 28, domSizes =  2,  3))
((vars =  3, 23, domSizes =  2,  2))
((vars =  3, 28, domSizes =  2,  3))
((vars =  4, 10, domSizes =  2,  2))
((vars =  4, 13, domSizes =  2,  3))
((vars =  4, 26, domSizes =  2,  2))
((vars =  5, 15, domSizes =  3,  2))
((vars =  5, 24, domSizes =  3,  2))
((vars =  5, 26, domSizes =  3,  2))
((vars =  5, 30, domSizes =  3,  1))
((vars =  6, 12, domSizes =  2,  3))
((vars =  6, 14, domSizes =  2,  2))
((vars =  6, 18, domSizes =  2,  2))
((vars =  6, 19, domSizes =  2,  2))
((vars =  6, 26, domSizes =  2,  2))
((vars =  7, 11, domSizes =  2,  3))
((vars =  7, 16, domSizes =  2,  2))
((vars =  7, 21, domSizes =  2,  1))
((vars =  8, 16, domSizes =  2,  2))
((vars =  8, 26, domSizes =  2,  2))
((vars =  9, 12, domSizes =  3,  3))
((vars =  9, 16, domSizes =  3,  2))
((vars =  9, 17, domSizes =  3,  2))
((vars = 10, 14, domSizes =  2,  2))
((vars = 10, 16, domSizes =  2,  2))
((vars = 10, 20, domSizes =  2,  2))
((vars = 10, 24, domSizes =  2,  2))
((vars = 10, 31, domSizes =  2,  2))
((vars = 11, 15, domSizes =  3,  2))
((vars = 11, 19, domSizes =  3,  2))
((vars = 11, 27, domSizes =  3,  3))
((vars = 12, 13, domSizes =  3,  3))
((vars = 12, 14, domSizes =  3,  2))
((vars = 12, 15, domSizes =  3,  2))
((vars = 12, 16, domSizes =  3,  2))
((vars = 13, 17, domSizes =  3,  2))
((vars = 13, 19, domSizes =  3,  2))
((vars = 13, 23, domSizes =  3,  2))
((vars = 13, 26, domSizes =  3,  2))
((vars = 14, 30, domSizes =  2,  1))
((vars = 15, 18, domSizes =  2,  2))
((vars = 16, 23, domSizes =  2,  2))
((vars = 16, 24, domSizes =  2,  2))
((vars = 16, 26, domSizes =  2,  2))
((vars = 16, 27, domSizes =  2,  3))
((vars = 17, 20, domSizes =  2,  2))
((vars = 17, 24, domSizes =  2,  2))
((vars = 19, 21, domSizes =  2,  1))
((vars = 19, 24, domSizes =  2,  2))
((vars = 20, 29, domSizes =  2,  1))
((vars = 21, 22, domSizes =  1,  2))
((vars = 22, 27, domSizes =  2,  3))
((vars = 23, 31, domSizes =  2,  2))
((vars = 25, 26, domSizes =  3,  2))
((vars = 26, 29, domSizes =  2,  1))
((vars = 27, 31, domSizes =  3,  2));

const auto tablesPtr = []{
  vector<Table<int>*> tmp;
  for (auto &table: tables) {
    tmp.push_back(&table);
  }
  return tmp;
}();

const Task<DummyOperations> task(tablesPtr.begin(), tablesPtr.end(), DummyOperations::CtorArgs());

const vector<int> clampRanks = list_of<int>
(2)(0)(2)(-1)(0)(1)(2)(1)(1)(1)(0)(1)(3)(0)(1)(3)(1)(1)(0)(1)(-1)(0)(0)(0)(0)(0)(0)(0)(1)(3)(1)(-1);

const VarVector expectedMinDegVarOrder   = list_of
    (25)(8)(29)(22)(28)(7)(30)(18)(23)(4)(27)(21)(10)(26)(14)(9)(24)(0)(2)(17)(15)(11)(19);
const VarVector expectedWMinDegVarOrder  = list_of
    (29)(25)(30)(8)(21)(28)(18)(4)(23)(10)(26)(14)(9)(6)(2)(17)(24)(27)(15)(22)(7)(1)(19);
const VarVector expectedMinFillVarOrder  = list_of
    (8)(22)(29)(25)(28)(7)(30)(18)(23)(4)(27)(1)(11)(10)(19)(9)(0)(17)(21)(2)(15)(24)(5)(26);
const VarVector expectedWMinFillVarOrder = list_of
    (8)(22)(29)(25)(28)(7)(30)(18)(23)(4)(27)(21)(10)(26)(14)(9)(24)(0)(2)(17)(15)(11)(19);

const vector<double> fngNums = list_of(0.0)(0.8)(0.5)(0.1);
const double maxComplexity = 4.0;
const double selectionScale = 1.5;

} // namespace (anonymous)::testData
} // anonymous namespace

BOOST_AUTO_TEST_SUITE( varorder )

BOOST_AUTO_TEST_CASE( emptyProblem )
{
  FixedNumberGenerator fng(testData::fngNums);
  const vector<Table<int>*> tables;
  Task<DummyOperations> task(tables.begin(), tables.end(), DummyOperations::CtorArgs());
  VarVector varOrder = greedyVarOrder(task, 1.0, vector<int>(), orang::greedyvarorder::MIN_DEGREE, fng, 1.0);

  BOOST_CHECK(varOrder.empty());
}

BOOST_AUTO_TEST_CASE( minDegree )
{
  FixedNumberGenerator fng(testData::fngNums);
  VarVector varOrder = greedyVarOrder(testData::task, testData::maxComplexity, testData::clampRanks,
      orang::greedyvarorder::MIN_DEGREE, fng, testData::selectionScale);

  BOOST_CHECK_EQUAL_COLLECTIONS(varOrder.begin(), varOrder.end(),
      testData::expectedMinDegVarOrder.begin(), testData::expectedMinDegVarOrder.end());

}

BOOST_AUTO_TEST_CASE( weightedMinDeg )
{
  FixedNumberGenerator fng(testData::fngNums);
  VarVector varOrder = greedyVarOrder(testData::task, testData::maxComplexity, testData::clampRanks,
      orang::greedyvarorder::WEIGHTED_MIN_DEGREE, fng, testData::selectionScale);

  BOOST_CHECK_EQUAL_COLLECTIONS(varOrder.begin(), varOrder.end(),
      testData::expectedWMinDegVarOrder.begin(), testData::expectedWMinDegVarOrder.end());
}

BOOST_AUTO_TEST_CASE( minFill )
{
  FixedNumberGenerator fng(testData::fngNums);
  VarVector varOrder = greedyVarOrder(testData::task, testData::maxComplexity, testData::clampRanks,
      orang::greedyvarorder::MIN_FILL, fng, testData::selectionScale);

  BOOST_CHECK_EQUAL_COLLECTIONS(varOrder.begin(), varOrder.end(),
      testData::expectedMinFillVarOrder.begin(), testData::expectedMinFillVarOrder.end());
}

BOOST_AUTO_TEST_CASE( weightedMinFill )
{
  FixedNumberGenerator fng(testData::fngNums);
  VarVector varOrder = greedyVarOrder(testData::task, testData::maxComplexity, testData::clampRanks,
      orang::greedyvarorder::WEIGHTED_MIN_FILL, fng, testData::selectionScale);

  BOOST_CHECK_EQUAL_COLLECTIONS(varOrder.begin(), varOrder.end(),
      testData::expectedWMinFillVarOrder.begin(), testData::expectedWMinFillVarOrder.end());
}

BOOST_AUTO_TEST_SUITE_END()

