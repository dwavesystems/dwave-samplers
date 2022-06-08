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
