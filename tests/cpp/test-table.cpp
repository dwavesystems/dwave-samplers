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
#include <ostream>
#include <limits>
#include <vector>
using std::size_t;
using std::copy;
using std::ostream;
using std::transform;
using std::vector;
using std::numeric_limits;

#include <boost/lambda/bind.hpp>
#include <boost/assign/list_of.hpp>
using boost::lambda::bind;
using boost::lambda::_1;
using boost::assign::list_of;

#include <base.h>
#include <table.h>
#include <exception.h>
using orang::Var;
using orang::DomIndex;
using orang::VarVector;
using orang::DomIndexVector;
using orang::SizeVector;
using orang::Table;
using orang::TableVar;
using orang::InvalidArgumentException;
using orang::LengthException;

namespace {
namespace goodTableData {
const size_t numVars = 5;
const size_t tableSize = 48;

const VarVector vars = list_of(0)(3)(7)(10)(11);
const DomIndexVector domSizes = list_of(2)(3)(2)(2)(2);
const SizeVector stepSizes = list_of(1)(2)(6)(12)(24);

const vector<int> ints = list_of(1)(2)(3)(4)(5)(6)(7)(8)(9)(10)(11)(12)
    (13)(14)(15)(16)(17)(18)(19)(20)(21)(22)(23)(24)
    (25)(26)(27)(28)(29)(30)(31)(32)(33)(34)(35)(36)
    (37)(38)(39)(40)(41)(42)(43)(44)(45)(46)(47)(48);
const vector<double> doubles = list_of(1.0)(2.0)(3.0)(4.0)(5.0)(6.0)(7.0)(8.0)(9.0)(10.0)(11.0)(12.0)
    (13.0)(14.0)(15.0)(16.0)(17.0)(18.0)(19.0)(20.0)(21.0)(22.0)(23.0)(24.0)
    (25.0)(26.0)(27.0)(28.0)(29.0)(30.0)(31.0)(32.0)(33.0)(34.0)(35.0)(36.0)
    (37.0)(38.0)(39.0)(40.0)(41.0)(42.0)(43.0)(44.0)(45.0)(46.0)(47.0)(48.0);
}

namespace badTableData {
const size_t numVars = 3;

const VarVector repeatedVars = list_of(1)(1)(4)(5)(6);
const VarVector unsortedVars = list_of(1)(6)(2)(3)(4);

const DomIndexVector shortDomSizes = list_of(2)(4);
const DomIndexVector zeroDomSize = list_of(2)(0)(2)(2)(2);
const DomIndexVector hugeDomSizes =
    list_of(numeric_limits<DomIndex>::max() - 1)
    (numeric_limits<DomIndex>::max() - 1)
    (numeric_limits<DomIndex>::max() - 1)
    (numeric_limits<DomIndex>::max() - 1)
    (numeric_limits<DomIndex>::max() - 1);
}
}

namespace orang {
ostream& operator<<(ostream& out, const TableVar& tv) {
  out << "TableVar(index=" << tv.index << ",domSize=" << tv.domSize << ",stepSize=" << tv.stepSize << ")";
  return out;
}
}


BOOST_AUTO_TEST_SUITE(table)

BOOST_AUTO_TEST_CASE( table_constructors )
{
  static const int zero = 0;

  // expected
  vector<TableVar> expectedVars;
  for (size_t i = 0; i < goodTableData::numVars; ++i) {
    expectedVars.push_back(TableVar(
        goodTableData::vars.at(i), goodTableData::domSizes.at(i), goodTableData::stepSizes.at(i)));
  }

  // empty table
  Table<int> emptyTable;
  BOOST_CHECK(emptyTable.vars().empty());
  BOOST_CHECK_EQUAL_COLLECTIONS(emptyTable.begin(), emptyTable.end(), &zero, &zero + 1);

  // non-empty table of ints
  Table<int> intTable(goodTableData::vars, goodTableData::domSizes);
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.vars().begin(), intTable.vars().end(),
      expectedVars.begin(), expectedVars.end());
  BOOST_CHECK_EQUAL(intTable.size(), goodTableData::tableSize);
  copy(goodTableData::ints.begin(), goodTableData::ints.end(), intTable.begin());
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.begin(), intTable.end(),
      goodTableData::ints.begin(), goodTableData::ints.end());

  // copy constructor
  Table<int> intTableCopy(intTable);
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.vars().begin(), intTable.vars().end(),
      intTableCopy.vars().begin(), intTableCopy.vars().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.begin(), intTable.end(), intTableCopy.begin(), intTableCopy.end());

  // copy constructor with implicit type conversion
  Table<double> doubleTable(intTable);
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.vars().begin(), intTable.vars().end(),
      doubleTable.vars().begin(), doubleTable.vars().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.begin(), intTable.end(), doubleTable.begin(), doubleTable.end());

  // assignment operator
  Table<int> intTableAssign;
  intTableAssign = intTable;
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.vars().begin(), intTable.vars().end(),
      intTableAssign.vars().begin(), intTableAssign.vars().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.begin(), intTable.end(), intTableAssign.begin(), intTableAssign.end());

  // assignment operator with implicit type conversion
  Table<double> doubleTableAssign;
  doubleTableAssign = intTable;
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.vars().begin(), intTable.vars().end(),
      doubleTableAssign.vars().begin(), doubleTableAssign.vars().end());
  BOOST_CHECK_EQUAL_COLLECTIONS(intTable.begin(), intTable.end(), doubleTableAssign.begin(), doubleTableAssign.end());
}

BOOST_AUTO_TEST_CASE( table_constructor_exceptions )
{
  BOOST_CHECK_THROW(Table<int> t(goodTableData::vars, badTableData::shortDomSizes), InvalidArgumentException);
  BOOST_CHECK_THROW(Table<int> t(badTableData::repeatedVars, goodTableData::domSizes), InvalidArgumentException);
  BOOST_CHECK_THROW(Table<int> t(badTableData::unsortedVars, goodTableData::domSizes), InvalidArgumentException);
  BOOST_CHECK_THROW(Table<int> t(goodTableData::vars, badTableData::zeroDomSize), InvalidArgumentException);
  BOOST_CHECK_THROW(Table<int> t(goodTableData::vars, badTableData::hugeDomSizes), LengthException);
}

BOOST_AUTO_TEST_SUITE_END()
