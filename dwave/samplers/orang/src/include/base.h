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
#ifndef INCLUDED_ORANG_BASE_H
#define INCLUDED_ORANG_BASE_H

#include <cstddef>
#include <set>
#include <vector>
#include <cstdint>

namespace orang {

typedef std::uint_least32_t Var;
typedef std::vector<Var> VarVector;
typedef std::set<Var> VarSet;

typedef std::uint_least16_t DomIndex;
typedef std::vector<DomIndex> DomIndexVector;

typedef std::vector<std::size_t> SizeVector;

} // namespace orang

#endif
