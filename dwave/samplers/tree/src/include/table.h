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
#ifndef INCLUDED_ORANG_TABLE_H
#define INCLUDED_ORANG_TABLE_H

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <limits>
#include <memory>

#include <exception.h>
#include <base.h>

namespace orang {

struct TableVar {
  Var index;
  DomIndex domSize;
  std::size_t stepSize;

  TableVar(Var index0, DomIndex domSize0, std::size_t stepSize0) :
    index(index0), domSize(domSize0), stepSize(stepSize0) {}
};


inline bool operator==(const TableVar& x, const TableVar& y) {
  return x.index == y.index && x.domSize == y.domSize && x.stepSize == y.stepSize;
}

inline bool operator!=(const TableVar& x, const TableVar& y) {
  return !(x == y);
}

inline bool operator<(const TableVar& x, const TableVar& y) {
  return x.index == y.index ? (x.domSize == y.domSize ?
      x.stepSize < y.stepSize : x.domSize < y.domSize) : x.index < y.index;
}

template<typename Y>
class Table {
public:
  typedef Y value_type;
  typedef typename std::vector<Y>::iterator iterator;
  typedef typename std::vector<Y>::const_iterator const_iterator;

  typedef std::shared_ptr<Table> smartptr;
  typedef std::shared_ptr<const Table> const_smartptr;

private:
  std::vector<TableVar> vars_;
  std::vector<value_type> values_;

public:

  //===========================================================================================================
  //
  //  Constructors & assignment operators
  //
  //===========================================================================================================

  Table() : vars_(), values_(1) {}

  Table(const VarVector& scope, const DomIndexVector& domSizes, const value_type& initVal = value_type()) :
    vars_(),
    values_() {

    if (scope.size() != domSizes.size()) {
      throw InvalidArgumentException("scope and domSizes vectors are not the same size");
    }

    vars_.reserve(scope.size());

    std::size_t coeffProd = 1;
    SizeVector::size_type limit = values_.max_size();
    for (std::size_t i = 0; i < scope.size(); ++i) {
      if (domSizes[i] == 0) {
        throw InvalidArgumentException("Domain size of zero encountered");
      }

      if (i > 0 && scope[i] <= scope[i - 1]) {
        throw InvalidArgumentException("Variables not listed in (strictly) increasing order");
      }

      vars_.push_back(TableVar(scope[i], domSizes[i], coeffProd));
      coeffProd *= domSizes[i];
      limit /= domSizes[i];
    }

    if (limit == 0) {
      throw LengthException();
    }

    values_.assign(coeffProd, initVal);
  }

  template<typename Z>
  Table(const Table<Z>& other) :
    vars_(other.vars()), values_(other.begin(), other.end()) {}

  template<typename Z>
  Table& operator=(const Table<Z>& other) {
    vars_ = other.vars();
    values_.reserve(other.size());
    values_.assign(other.begin(), other.end());
    return *this;
  }

  //===========================================================================================================
  //
  //  Comparison operators
  //
  //===========================================================================================================
  bool operator==(const Table& other) const {
    return vars_ == other.vars_ && values_ == other.values_;
  }

  bool operator!=(const Table& other) const {
    return !(*this == other);
  }

  bool operator<(const Table& other) const {
    return vars_ == other.vars_ ? values_ < other.values_ : vars_ < other.vars_;
  }

  const std::vector<TableVar>& vars() const {
    return vars_;
  }

  const TableVar& var(Var n) const {
    return vars_.at(n);
  }

  //===========================================================================================================
  //
  //  Iterators
  //
  //===========================================================================================================

  iterator begin() {
    return values_.begin();
  }

  const_iterator begin() const {
    return values_.begin();
  }

  iterator end() {
    return values_.end();
  }

  const_iterator end() const {
    return values_.end();
  }

  //===========================================================================================================
  //
  //  Value element access
  //
  //===========================================================================================================

  const value_type& operator[](std::size_t n) const {
    return values_[n];
  }

  value_type& operator[](std::size_t n) {
    return values_[n];
  }

  //===========================================================================================================
  //
  //  Other vector-like functions
  //
  //===========================================================================================================

  std::size_t size() const {
    return values_.size();
  }
};

} // namespace orang

#endif
