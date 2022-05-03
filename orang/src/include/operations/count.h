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
#ifndef INCLUDED_ORANG_OPERATIONS_COUNT_H
#define INCLUDED_ORANG_OPERATIONS COUNT_H

#include <cmath>

#include <base.h>
#include <exception.h>
#include <table.h>
#include <marginalizer.h>

namespace orang {

template<typename Y>
class ValueCount {
public:
  typedef Y value_type;
private:
  value_type value_;
  double count_;
public:
  ValueCount(value_type value = value_type(), double count = 1.0) : value_(value), count_(count) {}
  value_type value() const { return value_; }
  double count() const { return count_; }
};

template<typename Y>
class ValueCountEpsCompare {
private:
  double eps_;
public:
  ValueCountEpsCompare(double eps) : eps_(eps) {}
  bool operator()(const Y& a, const Y&b) { return a < b - eps_ * (b > Y(0) ? b : -b); }
};

template<typename Y>
class CountMarginalizer : public Marginalizer<ValueCount<Y> > {
private:
  typedef ValueCount<Y> value_type;
  typedef Table<value_type> table_type;

  ValueCountEpsCompare<Y> cmp_;

  virtual value_type marginalizeImpl(std::size_t, const table_type& mrgTable) {
    typename value_type::value_type bestValue = mrgTable.begin()->value();
    double totalCount = 0.0;
    for (const auto &v: mrgTable) {
      if (cmp_(v.value(), bestValue)) {
        bestValue = v.value();
        totalCount = v.count();
      } else if (!cmp_(bestValue, v.value())) {
        totalCount += v.count();
      }
    }
    return value_type(bestValue, totalCount);
  }
public:
  CountMarginalizer(ValueCountEpsCompare<Y> cmp) : cmp_(cmp) {}
};

template<typename Y>
class CountOperations {
public:
  typedef ValueCount<Y> value_type;
  typedef ValueCountEpsCompare<Y> value_compare;
  typedef int solution_type;
  typedef MarginalizerTypes<value_type, solution_type> marginalizer_types;
  typedef typename marginalizer_types::marginalizer_type marginalizer_type;
  typedef typename marginalizer_types::solvablemarginalizer_type solvablemarginalizer_type;
  typedef typename marginalizer_types::marginalizer_smartptr marginalizer_smartptr;
  typedef typename marginalizer_types::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;

  struct CtorArgs {
    double eps;
    CtorArgs(double eps0) : eps(eps0) {}
  };

private:
  marginalizer_smartptr marginalizer_;

public:
  static value_type combineIdentity() { return value_type(0, 1); }
  static value_type combine(const value_type& x, const value_type& y) {
    return value_type(x.value() + y.value(), x.count() * y.count());
  }
  static value_type combineInverse(const value_type& c, const value_type& x) {
    return value_type(c.value() - x.value(), c.count() / x.count());
  }

  explicit CountOperations(const CtorArgs& ctorArgs) :
    marginalizer_(new CountMarginalizer<Y>(value_compare(ctorArgs.eps))) {}

  marginalizer_smartptr marginalizer() const {
    return marginalizer_;
  }

  solvablemarginalizer_smartptr solvableMarginalizer(
      const VarVector&, const DomIndexVector&,
      Var, DomIndex) const {
    throw OperationUnavailable();
  }

  solution_type initSolution(const DomIndexVector&) const {
    throw OperationUnavailable();
  }
};

} // namespace orang

#endif
