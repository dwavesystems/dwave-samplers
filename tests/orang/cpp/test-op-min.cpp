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
#include <algorithm>
#include <numeric>
#include <vector>
#include <functional>

#include <boost/assign/list_of.hpp>

#include <base.h>
#include <table.h>
#include <combine.h>
#include <marginalizer.h>
#include <operations/min.h>

#include "test.h"

using std::size_t;
using std::copy;
using std::inner_product;
using std::vector;
using std::greater;

using boost::assign::list_of;

using orang::Var;
using orang::VarVector;
using orang::DomIndexVector;
using orang::Table;
using orang::Plus;
using orang::MinSolutionSet;
using orang::MinMarginalizer;
using orang::SolvableMinMarginalizer;
using orang::MinOperations;

using minsolsetAssign::MinSolSetBuilder;

typedef MinOperations<int, Plus<int> > ops_type;
typedef ops_type::marginalizer_type marginalizer_type;
typedef ops_type::marginalizer_smartptr marginalizer_smartptr;
typedef ops_type::solvablemarginalizer_type solvablemarginalizer_type;
typedef ops_type::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;

typedef MinOperations<int, Plus<int>, greater<int> > maxops_type;

namespace {
namespace tableData {

const vector<int> values = list_of(-2)(5)(1)(-3)(-4)(-1)(-2)(6);
const Var outVar = 7;
const size_t outDomSize = 8;

const VarVector inScope = list_of(1)(4);
const DomIndexVector inDomSizes = list_of(3)(2);
const VarVector inVarIndices = list_of(1)(1);
const size_t inIndex = 4;
const int expectedMinMrgValue = -4;
const int expectedMaxMrgValue = 6;

const MinSolutionSet<int> inSol3 = MinSolSetBuilder<int>(3)
    (100, list_of(9)(1)(9)(9)(1)(9)(9)(9)(9)(9))
    (101, list_of(8)(1)(8)(8)(1)(8)(8)(8)(8)(8));
const MinSolutionSet<int> expectedOutSol3 = MinSolSetBuilder<int>(3)
    (100, list_of(9)(1)(9)(9)(1)(9)(9)(4)(9)(9))
    (101, list_of(8)(1)(8)(8)(1)(8)(8)(4)(8)(8))
    (101, list_of(9)(1)(9)(9)(1)(9)(9)(3)(9)(9));

const MinSolutionSet<int> inSol10 = MinSolSetBuilder<int>(10)
    (200, list_of(9)(1)(9)(9)(1)(9)(9)(9)(9)(9));
const MinSolutionSet<int> expectedOutSol10 = MinSolSetBuilder<int>(10)
    (200, list_of(9)(1)(9)(9)(1)(9)(9)(4)(9)(9))
    (201, list_of(9)(1)(9)(9)(1)(9)(9)(3)(9)(9))
    (202, list_of(9)(1)(9)(9)(1)(9)(9)(0)(9)(9))
    (202, list_of(9)(1)(9)(9)(1)(9)(9)(6)(9)(9))
    (203, list_of(9)(1)(9)(9)(1)(9)(9)(5)(9)(9))
    (205, list_of(9)(1)(9)(9)(1)(9)(9)(2)(9)(9))
    (209, list_of(9)(1)(9)(9)(1)(9)(9)(1)(9)(9))
    (210, list_of(9)(1)(9)(9)(1)(9)(9)(7)(9)(9));

const MinSolutionSet<int, greater<int> > inSolMax = MinSolSetBuilder<int, greater<int> >(5)
    (300, list_of(9)(1)(9)(9)(1)(9)(9)(9)(9)(9));
const MinSolutionSet<int, greater<int> > expectedOutSolMax = MinSolSetBuilder<int, greater<int> >(5)
    (300, list_of(9)(1)(9)(9)(1)(9)(9)(7)(9)(9))
    (299, list_of(9)(1)(9)(9)(1)(9)(9)(1)(9)(9))
    (295, list_of(9)(1)(9)(9)(1)(9)(9)(2)(9)(9))
    (293, list_of(9)(1)(9)(9)(1)(9)(9)(5)(9)(9))
    (292, list_of(9)(1)(9)(9)(1)(9)(9)(0)(9)(9));


} // namespace tableData
} // anonymous namespace


BOOST_AUTO_TEST_SUITE( min_marginalizer )

BOOST_AUTO_TEST_CASE( marginalizer )
{
  ops_type ops;
  marginalizer_smartptr mrgP = ops.marginalizer();
  marginalizer_type& mrg = *mrgP;

  Table<int> mrgTable(VarVector(1, tableData::outVar), DomIndexVector(1, tableData::outDomSize));
  BOOST_REQUIRE_EQUAL(tableData::values.size(), mrgTable.size());
  copy(tableData::values.begin(), tableData::values.end(), mrgTable.begin());

  BOOST_CHECK_EQUAL(mrg(tableData::inIndex, mrgTable), tableData::expectedMinMrgValue);
}

BOOST_AUTO_TEST_CASE( solvable_marginalizer_hitmaxsols )
{
  ops_type ops;
  solvablemarginalizer_smartptr mrgP = ops.solvableMarginalizer(tableData::inScope, tableData::inDomSizes,
      tableData::outVar, tableData::outDomSize);
  solvablemarginalizer_type& mrg = *mrgP;

  Table<int> mrgTable(VarVector(1, tableData::outVar), DomIndexVector(1, tableData::outDomSize));
  BOOST_REQUIRE_EQUAL(tableData::values.size(), mrgTable.size());
  copy(tableData::values.begin(), tableData::values.end(), mrgTable.begin());

  BOOST_CHECK_EQUAL(mrg(tableData::inIndex, mrgTable), tableData::expectedMinMrgValue);
  MinSolutionSet<int> outSol = tableData::inSol3;
  mrg.solve(outSol);
  BOOST_CHECK_EQUAL_COLLECTIONS(outSol.solutions().begin(), outSol.solutions().end(),
      tableData::expectedOutSol3.solutions().begin(), tableData::expectedOutSol3.solutions().end());
}

BOOST_AUTO_TEST_CASE( solvable_marginalizer_nohitmaxsols )
{
  ops_type ops;
  solvablemarginalizer_smartptr mrgP = ops.solvableMarginalizer(tableData::inScope, tableData::inDomSizes,
      tableData::outVar, tableData::outDomSize);
  solvablemarginalizer_type& mrg = *mrgP;

  Table<int> mrgTable(VarVector(1, tableData::outVar), DomIndexVector(1, tableData::outDomSize));
  BOOST_REQUIRE_EQUAL(tableData::values.size(), mrgTable.size());
  copy(tableData::values.begin(), tableData::values.end(), mrgTable.begin());

  BOOST_CHECK_EQUAL(mrg(tableData::inIndex, mrgTable), tableData::expectedMinMrgValue);
  MinSolutionSet<int> outSol = tableData::inSol10;
  mrg.solve(outSol);
  BOOST_CHECK_EQUAL_COLLECTIONS(outSol.solutions().begin(), outSol.solutions().end(),
      tableData::expectedOutSol10.solutions().begin(), tableData::expectedOutSol10.solutions().end());
}

BOOST_AUTO_TEST_CASE( solvable_marginalizer_maxcompare )
{
  maxops_type ops;
  maxops_type::solvablemarginalizer_smartptr mrgP = ops.solvableMarginalizer(tableData::inScope, tableData::inDomSizes,
      tableData::outVar, tableData::outDomSize);
  maxops_type::solvablemarginalizer_type& mrg = *mrgP;

  Table<int> mrgTable(VarVector(1, tableData::outVar), DomIndexVector(1, tableData::outDomSize));
  BOOST_REQUIRE_EQUAL(tableData::values.size(), mrgTable.size());
  copy(tableData::values.begin(), tableData::values.end(), mrgTable.begin());

  BOOST_CHECK_EQUAL(mrg(tableData::inIndex, mrgTable), tableData::expectedMaxMrgValue);
  MinSolutionSet<int, greater<int> > outSol = tableData::inSolMax;
  mrg.solve(outSol);
  BOOST_CHECK_EQUAL_COLLECTIONS(outSol.solutions().begin(), outSol.solutions().end(),
      tableData::expectedOutSolMax.solutions().begin(), tableData::expectedOutSolMax.solutions().end());
}

BOOST_AUTO_TEST_SUITE_END()
