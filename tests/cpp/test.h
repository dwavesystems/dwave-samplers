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
#ifndef INCLUDED_TEST_H
#define INCLUDED_TEST_H

#include <cstddef>
#include <algorithm>
#include <ostream>
#include <vector>
#include <functional>
#include <base.h>
#include <table.h>
#include <graph.h>
#include <operations/min.h>

class FixedNumberGenerator {
private:
  const std::vector<double> x_;
  std::vector<double>::const_iterator it_;
public:
  FixedNumberGenerator(const std::vector<double>& x) : x_(x), it_(x_.begin()) {}
  double operator()() { if (it_ == x_.end()) it_ = x_.begin(); return *it_++; }
};

namespace orang {
bool operator==(const Graph& g1, const Graph& g2);
std::ostream& operator<<(std::ostream& out, const Graph& g);
template<typename Y> std::ostream& operator<<(std::ostream& out, const Table<Y>& t);

template<typename Y> bool operator==(const orang::MinSolution<Y>& s1, const orang::MinSolution<Y>& s2) {
  return s1.value == s2.value && s1.solution == s2.solution;
}
template<typename Y> bool operator!=(const orang::MinSolution<Y>& s1, const orang::MinSolution<Y>& s2) {
  return !(s1 == s2);
}

template<typename Y> std::ostream& operator<<(std::ostream& out, const orang::MinSolution<Y>& s);
template<typename Y> std::ostream& operator<<(std::ostream& out, const orang::MinSolutionSet<Y>& s);
}

namespace tableAssign {

struct WrongNumberOfDomSizes {};
struct WrongNumberOfValues {};

template<typename Y>
class ValuesBuilder {
private:
  const orang::VarVector& vars_;
  const orang::DomIndexVector& domSizes_;
  std::vector<Y> values_;
public:
  ValuesBuilder(const orang::VarVector& vars, const orang::DomIndexVector& domSizes, Y y) :
    vars_(vars), domSizes_(domSizes), values_(1, y) {}
  ValuesBuilder<Y>& operator,(Y y) { values_.push_back(y); return *this; }
  template<typename Z>
  operator orang::Table<Z>() {
    orang::Table<Z> t(vars_, domSizes_);
    if (t.size() != values_.size()) { throw WrongNumberOfValues(); }
    std::copy(values_.begin(), values_.end(), t.begin());
    return t;
  }
};

template<typename Y>
class ValuesBuilderInit {
private:
  Y value_;
public:
  ValuesBuilderInit(Y y) : value_(y) {}
  Y value() const { return value_; }
};

struct ValuesSentinel {
  template<typename Y>
  ValuesBuilderInit<Y> operator=(Y y) const { return ValuesBuilderInit<Y>(y); }
};

class DomSizesBuilder {
private:
  const orang::VarVector& vars_;
  orang::DomIndexVector domSizes_;
public:
  DomSizesBuilder(const orang::VarVector& vars, const orang::DomIndexVector& domSizes = orang::DomIndexVector()) :
    vars_(vars), domSizes_(domSizes) {}

  DomSizesBuilder& operator,(orang::DomIndex s) { domSizes_.push_back(s); return *this; }

  template<typename Y>
  ValuesBuilder<Y> operator,(const ValuesBuilderInit<Y>& b) const {
    if (vars_.size() != domSizes_.size()) {
      throw WrongNumberOfDomSizes();
    }
    return ValuesBuilder<Y>(vars_, domSizes_, b.value());
  }

  template<typename Y>
  operator orang::Table<Y>() const {
    if (vars_.size() != domSizes_.size()) {
      throw WrongNumberOfDomSizes();
    }
    return orang::Table<Y>(vars_, domSizes_);
  }
};

class DomSizesBuilderInit {
private:
  orang::DomIndexVector domSizes_;
public:
  DomSizesBuilderInit() : domSizes_() {}
  DomSizesBuilderInit(orang::DomIndex s) : domSizes_(1, s) {}
  const orang::DomIndexVector& domSizes() const { return domSizes_; }
};

struct NullSentinel {};

struct DomSizesSentinel {
  DomSizesBuilderInit operator=(orang::DomIndex s) const { return DomSizesBuilderInit(s); }
  DomSizesBuilderInit operator=(const NullSentinel&) const { return DomSizesBuilderInit(); }
};

class VarsBuilder {
private:
  orang::VarVector vars_;
public:
  VarsBuilder() : vars_() {}
  VarsBuilder(orang::Var v) : vars_(1, v) {}
  VarsBuilder& operator,(orang::Var v) { vars_.push_back(v); return *this; }
  DomSizesBuilder operator,(const DomSizesBuilderInit& b) { return DomSizesBuilder(vars_, b.domSizes()); }
};

struct VarsSentinel {
  VarsBuilder operator=(orang::Var v) const { return VarsBuilder(v); }
  VarsBuilder operator=(const NullSentinel&) const { return VarsBuilder(); }
};

extern const VarsSentinel vars;
extern const DomSizesSentinel domSizes;
extern const ValuesSentinel values;
extern const NullSentinel none;

} // namespace tableAssign



namespace minsolsetAssign {

template<typename Y, typename Compare=std::less<Y> >
class MinSolSetBuilder {
private:
  orang::MinSolutionSet<Y, Compare> minSolSet_;

public:
  MinSolSetBuilder(std::size_t maxSols) : minSolSet_(maxSols) {}

  MinSolSetBuilder& operator()(Y value, orang::DomIndexVector sol) {
    minSolSet_.solutions().insert(orang::MinSolution<Y>(value, sol));
    return *this;
  }

  operator orang::MinSolutionSet<Y, Compare>() const {
    return minSolSet_;
  }
};

}

#endif
