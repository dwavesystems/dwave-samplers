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

#include <boost/assign/list_of.hpp>

#include <base.h>
#include <table.h>
#include <marginalizer.h>
#include <operations/logsumprod.h>

#include "test.h"

using std::size_t;
using std::copy;
using std::inner_product;
using std::vector;

using boost::assign::list_of;

using orang::Var;
using orang::VarVector;
using orang::DomIndex;
using orang::DomIndexVector;
using orang::Table;
using orang::LogSumMarginalizer;
using orang::SolvableLogSumMarginalizer;
using orang::LogSumProductOperations;

namespace {

typedef LogSumProductOperations<FixedNumberGenerator> ops_type;
typedef ops_type::marginalizer_type marginalizer_type;
typedef ops_type::marginalizer_smartptr marginalizer_smartptr;
typedef ops_type::solvablemarginalizer_type solvablemarginalizer_type;
typedef ops_type::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;
typedef ops_type::CtorArgs ops_ctorargs;

namespace tableData {

const vector<int> values = list_of(0.0)(1.0)(-2.0)(-1.0)(0.0)(0.0)(2.0)(1.0)(0.0);
const Var outVar = 2;
const DomIndex outDomSize = 9;

const VarVector inScope = list_of(0)(6);
const DomIndexVector inDomSizes = list_of(4)(2);
const VarVector inVarIndices = list_of(2)(0);
const size_t inIndex = 2;
const double expectedMinMrgValue = 2.85237185;

const DomIndexVector inSol = list_of(2)(9)(9)(9)(9)(9)(0);

const vector<double> fixedNums = list_of(0.22)(0.23)(0.359)(0.4);

const VarVector expectedOutSol1 = list_of(2)(9)(2)(9)(9)(9)(0);
const VarVector expectedOutSol2 = list_of(2)(9)(3)(9)(9)(9)(0);
const VarVector expectedOutSol3 = list_of(2)(9)(5)(9)(9)(9)(0);
const VarVector expectedOutSol4 = list_of(2)(9)(6)(9)(9)(9)(0);
const vector<const VarVector*> expectedOutSols =
    list_of(&expectedOutSol1)(&expectedOutSol2)(&expectedOutSol3)(&expectedOutSol4);

}
}


BOOST_AUTO_TEST_SUITE( logsumprod_marginalizer )

BOOST_AUTO_TEST_CASE( marginalizer )
{
  FixedNumberGenerator fng(tableData::fixedNums);
  ops_type ops(fng);
  marginalizer_smartptr mrgP = ops.marginalizer();
  marginalizer_type& mrg = *mrgP;

  Table<int> mrgTable(VarVector(1, tableData::outVar), DomIndexVector(1, tableData::outDomSize));
  BOOST_REQUIRE_EQUAL(tableData::values.size(), mrgTable.size());
  copy(tableData::values.begin(), tableData::values.end(), mrgTable.begin());

  BOOST_CHECK_CLOSE(mrg(tableData::inIndex, mrgTable), tableData::expectedMinMrgValue, 1e-6);
}

BOOST_AUTO_TEST_CASE( solvable_marginalizer )
{
  FixedNumberGenerator fng(tableData::fixedNums);
  ops_type ops(fng);
  solvablemarginalizer_smartptr mrgP = ops.solvableMarginalizer(tableData::inScope, tableData::inDomSizes,
      tableData::outVar, tableData::outDomSize);
  solvablemarginalizer_type& mrg = *mrgP;

  Table<int> mrgTable(VarVector(1, tableData::outVar), DomIndexVector(1, tableData::outDomSize));
  BOOST_REQUIRE_EQUAL(tableData::values.size(), mrgTable.size());
  copy(tableData::values.begin(), tableData::values.end(), mrgTable.begin());

  BOOST_CHECK_CLOSE(mrg(tableData::inIndex, mrgTable), tableData::expectedMinMrgValue, 1e-6);

  for (const auto &expectedOutSol: tableData::expectedOutSols) {
    DomIndexVector outSol = tableData::inSol;
    mrg.solve(outSol);
    BOOST_CHECK_EQUAL_COLLECTIONS(outSol.begin(), outSol.end(),
        expectedOutSol->begin(), expectedOutSol->end());
  }
}

BOOST_AUTO_TEST_SUITE_END()
