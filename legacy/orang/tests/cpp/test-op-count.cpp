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
#include <functional>
#include <ostream>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/assign/list_of.hpp>

#include <base.h>
#include <table.h>
#include <marginalizer.h>
#include <operations/count.h>

using std::less;
using std::ostream;
using std::vector;

using boost::assign::list_of;

using orang::ValueCount;
using orang::CountOperations;
using orang::Table;
using orang::VarVector;
using orang::DomIndexVector;
using orang::Var;

typedef CountOperations<int> ops_type;
typedef ops_type::marginalizer_type marginalizer_type;
typedef ops_type::marginalizer_smartptr marginalizer_smartptr;
typedef ops_type::solvablemarginalizer_type solvablemarginalizer_type;
typedef ops_type::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;

namespace orang {
template<typename Y> bool operator==(const orang::ValueCount<Y>& a, const orang::ValueCount<Y>& b) {
  return a.value() == b.value() && a.count() == b.count();
}

template<typename Y> bool operator!=(const orang::ValueCount<Y>& a, const orang::ValueCount<Y>& b) {
  return !(a == b);
}

template<typename Y>
ostream& operator<<(ostream& o, const ValueCount<Y>& v) {
  return o << "(value=" << v.value() << ",count=" << v.count() << ")";
}

}

BOOST_AUTO_TEST_SUITE( count_operations )

BOOST_AUTO_TEST_CASE( combine ) {
  ops_type::value_type v1 = ops_type::combineIdentity();
  ops_type::value_type v2(-10, 100.0);
  ops_type::value_type v3(2, 4.0);
  ops_type::value_type ev23(-8, 400.0);

  BOOST_CHECK_EQUAL(ops_type::combine(v1, v2), v2);
  BOOST_CHECK_EQUAL(ops_type::combine(v3, v2), ev23);
  BOOST_CHECK_EQUAL(ops_type::combineInverse(ev23, v2), v3);
  BOOST_CHECK_EQUAL(ops_type::combineInverse(v3, v1), v3);
}

BOOST_AUTO_TEST_CASE( marginalizer_first )
{
  ops_type ops(0.0);
  marginalizer_smartptr mrgP = ops.marginalizer();
  marginalizer_type& mrg = *mrgP;

  Table<ValueCount<int> > mrgTable(VarVector(1, 10000), DomIndexVector(1, 6));
  BOOST_REQUIRE_EQUAL(mrgTable.size(), 6);
  mrgTable[0] = ValueCount<int>(-1);
  mrgTable[1] = ValueCount<int>(2, 100);
  mrgTable[2] = ValueCount<int>(-1, 20);
  mrgTable[3] = ValueCount<int>(0, 100);
  mrgTable[4] = ValueCount<int>(10, 100);
  mrgTable[5] = ValueCount<int>(0, 100);

  BOOST_CHECK_EQUAL(mrg(2000, mrgTable), ValueCount<int>(-1, 21));
}

BOOST_AUTO_TEST_CASE( marginalizer_tail )
{
  ops_type ops(0.0);
  marginalizer_smartptr mrgP = ops.marginalizer();
  marginalizer_type& mrg = *mrgP;

  Table<ValueCount<int> > mrgTable(VarVector(1, 10000), DomIndexVector(1, 6));
  BOOST_REQUIRE_EQUAL(mrgTable.size(), 6);
  mrgTable[0] = ValueCount<int>(-1);
  mrgTable[1] = ValueCount<int>(2, 100);
  mrgTable[2] = ValueCount<int>(-1, 20);
  mrgTable[3] = ValueCount<int>(0, 100);
  mrgTable[4] = ValueCount<int>(-10, 50);
  mrgTable[5] = ValueCount<int>(-10, 5);

  BOOST_CHECK_EQUAL(mrg(2000, mrgTable), ValueCount<int>(-10, 55));
}

BOOST_AUTO_TEST_CASE( marginalizer_eps )
{
  ops_type ops(1e-3);
  marginalizer_smartptr mrgP = ops.marginalizer();
  marginalizer_type& mrg = *mrgP;

  Table<ValueCount<int> > mrgTable(VarVector(1, 10000), DomIndexVector(1, 6));
  BOOST_REQUIRE_EQUAL(mrgTable.size(), 6);
  mrgTable[0] = ValueCount<int>(-1);
  mrgTable[1] = ValueCount<int>(2, 100);
  mrgTable[2] = ValueCount<int>(-1, 20);
  mrgTable[3] = ValueCount<int>(-10.001, 50);
  mrgTable[4] = ValueCount<int>(-9.98, 100);
  mrgTable[5] = ValueCount<int>(-10, 5);

  BOOST_CHECK_EQUAL(mrg(2000, mrgTable), ValueCount<int>(-10, 55));

}

BOOST_AUTO_TEST_SUITE_END()
