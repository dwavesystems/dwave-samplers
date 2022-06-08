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
#ifndef ORANG_COMBINE_H
#define ORANG_COMBINE_H

namespace orang {

template<typename Y>
struct Plus {
  typedef Y value_type;
  static value_type combineIdentity() { return value_type(0); }
  static value_type combine(const value_type& x, const value_type& y) { return x + y; }
  static value_type combineInverse(const value_type& c, const value_type& x) { return c - x; }
};

template<typename Y>
struct Multiply {
  typedef Y value_type;
  static value_type combineIdentity() { return value_type(1); }
  static value_type combine(const value_type& x, const value_type& y) { return x * y; }
  static value_type combineInverse(const value_type& c, const value_type& x) { return c / x; }
};

}

#endif
