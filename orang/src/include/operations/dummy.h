#ifndef INCLUDED_ORANG_OPERATIONS_DUMMY_H
#define INCLUDED_ORANG_OPERATIONS_DUMMY_H

#include <exception.h>
#include <marginalizer.h>

namespace orang {

class DummyOperations {
public:
  typedef int value_type;
  typedef int solution_type;
  typedef MarginalizerTypes<value_type,solution_type> marginalizer_types;
  typedef marginalizer_types::marginalizer_type marginalizer_type;
  typedef marginalizer_types::solvablemarginalizer_type solvablemarginalizer_type;
  typedef marginalizer_types::marginalizer_smartptr marginalizer_smartptr;
  typedef marginalizer_types::solvablemarginalizer_smartptr solvablemarginalizer_smartptr;

  struct CtorArgs {};

  DummyOperations(const CtorArgs&) {}
  value_type combine(const value_type&, const value_type&) const { throw OperationUnavailable(); }
  value_type combineIdentity() const { throw OperationUnavailable(); }
  marginalizer_smartptr marginalizer() const { throw OperationUnavailable(); }
  solvablemarginalizer_smartptr solvableMarginalizer(
      const VarVector&, const SizeVector&, Var, DomIndex) const {
    throw OperationUnavailable();
  }
  solution_type initSolution(const DomIndexVector&) const { throw OperationUnavailable(); }
};

} // namespace orang

#endif
