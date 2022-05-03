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

#include <vector>

#include <boost/assign/list_of.hpp>

#include <table.h>
#include <task.h>
#include <combine.h>
#include <operations/min.h>
#include <merger.h>

#include "test.h"

using std::vector;

using boost::assign::list_of;

using orang::InvalidArgumentException;
using orang::VarVector;
using orang::Table;
using orang::Plus;
using orang::MinOperations;
using orang::Task;
using orang::TableMerger;

using tableAssign::vars;
using tableAssign::domSizes;
using tableAssign::values;
using tableAssign::none;

namespace testData {

vector<Table<int> > inTables = list_of<Table<int> >
((vars = none, domSizes = none, values = 9))
((vars = 0, 1, 2, domSizes = 2, 2, 2, values =  6, 9, 3, -9, 7, 8, 3, 5))
((vars = 0, 4, 6, domSizes = 2, 3, 2, values = 5, -2, 3, -6, 4, -9, -4, -9, -8, 6, 4, -3))
((vars = 1, 2, 3, 5, domSizes = 2, 2, 2, 4, values = 9, -9, -1, -2, 5, 6, -6, 0, -1, 3, 4, 5, -4, 3, 3, -6,
    -7, 0, 9, -3, 2, -5, 5, -5, 0, 4, 7, 9, 1, -7, -7, -5))
((vars = 3, 4, domSizes = 2, 3, values = 6, -5, 6, -5, 8, -3))
((vars = 5, 6, domSizes = 4, 2, values = -6, -5, 2, -1, -3, 6, 2, 1));

const auto inTablesPtr = []{
  vector<Table<int>*> tmp;
  for (auto &table: inTables) {
    tmp.push_back(&table);
  }
  return tmp;
}();


VarVector outVars = list_of(0)(4)(6);

Table<int> expectedTable = (vars = 0, 4, 6, domSizes = 2, 3, 2,
    values = 1, -15, -1, -19, 2, -20, -3, -20, -7, -5, 7, -12);

Table<int> expectedTrivialTable = (vars = none, domSizes = none, values = -20);

vector<Table<int>*> emptyTables;
VarVector emptyOutVars;
Table<int> expectedEmptyTrivialTable = (vars = none, domSizes = none, values = 0);


}

BOOST_AUTO_TEST_SUITE( merger )

typedef MinOperations<int, Plus<int> > ops_type;
typedef Task<ops_type> task_type;
typedef task_type::const_table_smartptr const_table_smartptr;
typedef task_type::marginalizer_smartptr marginalizer_smartptr;

BOOST_AUTO_TEST_CASE( merger )
{
  task_type task(testData::inTablesPtr.begin(), testData::inTablesPtr.end(), 1);
  TableMerger<task_type> merge(task);
  marginalizer_smartptr marginalizer = task.marginalizer();
  const_table_smartptr outTable = merge(testData::outVars, testData::inTablesPtr.begin(), testData::inTablesPtr.end(), *marginalizer);

  BOOST_CHECK_EQUAL(*outTable, testData::expectedTable);
}

BOOST_AUTO_TEST_CASE( merge_to_nullscope )
{
  task_type task(testData::inTablesPtr.begin(), testData::inTablesPtr.end(), 1);
  TableMerger<task_type> merge(task);
  marginalizer_smartptr marginalizer = task.marginalizer();
  const_table_smartptr outTable = merge(VarVector(), testData::inTablesPtr.begin(), testData::inTablesPtr.end(), *marginalizer);

  BOOST_CHECK_EQUAL(*outTable, testData::expectedTrivialTable);
}

BOOST_AUTO_TEST_CASE( empty_merger )
{
  task_type task(testData::inTablesPtr.begin(), testData::inTablesPtr.end(), 1);
  TableMerger<task_type> merge(task);
  marginalizer_smartptr marginalizer = task.marginalizer();
  const_table_smartptr outTable = merge(testData::emptyOutVars,
      testData::emptyTables.begin(), testData::emptyTables.end(), *marginalizer);
  BOOST_CHECK_EQUAL(*outTable, testData::expectedEmptyTrivialTable);
}


BOOST_AUTO_TEST_SUITE_END()
