#ifndef INCLUDED_ORANG_MARGINALIZER_H
#define INCLUDED_ORANG_MARGINALIZER_H

#include <cstddef>

#include <boost/shared_ptr.hpp>

#include <orang/base.h>
#include <orang/table.h>

namespace orang {

template<typename Y>
class Marginalizer {
public:
  typedef Y value_type;
  typedef Table<value_type> table_type;

private:
  virtual value_type marginalizeImpl(std::size_t outIndex, const table_type& mrgTable) = 0;

public:
  virtual ~Marginalizer() {}
  value_type operator()(std::size_t outIndex, const table_type& mrgTable) {
    return marginalizeImpl(outIndex, mrgTable);
  }
};

template<typename Y, typename S>
class SolvableMarginalizer : public Marginalizer<Y> {
public:
  typedef Y value_type;
  typedef S solution_type;

protected:
  typedef std::pair<Var, std::size_t> varstep_pair;
  typedef std::vector<varstep_pair> varstep_vector;

  static std::size_t buildStepSizes(const VarVector& scope, const DomIndexVector& domSizes, varstep_vector& varsSteps) {
    using std::make_pair;
    using std::size_t;

    size_t limit = std::numeric_limits<size_t>::max();
    size_t stepSize = 1;
    varsSteps.reserve(domSizes.size());
    VarVector::const_iterator scopeIter = scope.begin();
    VarVector::const_iterator scopeEnd = scope.end();
    DomIndexVector::const_iterator domSizesIter = domSizes.begin();
    while (scopeIter != scopeEnd) {
      limit /= *domSizesIter;
      if (limit == 0) {
        throw LengthException();
      }
      varsSteps.push_back(make_pair(*scopeIter, stepSize));
      stepSize *= *domSizesIter;

      ++scopeIter;
      ++domSizesIter;
    }

    return stepSize;
  }

private:
  virtual void solveImpl(solution_type& s) const = 0;

public:
  virtual ~SolvableMarginalizer() {}
  void solve(solution_type& s) const {
    solveImpl(s);
  }
};

template<typename Y, typename S>
struct MarginalizerTypes {
  typedef Y value_type;
  typedef S solution_type;
  typedef Marginalizer<Y> marginalizer_type;
  typedef SolvableMarginalizer<Y,S> solvablemarginalizer_type;
  typedef boost::shared_ptr<marginalizer_type> marginalizer_smartptr;
  typedef boost::shared_ptr<solvablemarginalizer_type> solvablemarginalizer_smartptr;
};

} // namespace orang

#endif
