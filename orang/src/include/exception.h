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
#ifndef INCLUDED_ORANG_EXCEPTION_H
#define INCLUDED_ORANG_EXCEPTION_H

#include <string>

namespace orang {

class Exception {
private:
  std::string msg_;
public:
  Exception(const std::string& msg = "orang::Exception") : msg_(msg) {}

  const std::string& what() const { return msg_; }
};

class LengthException : public Exception {
public:
  LengthException(const std::string& msg = "orang::LengthException") : Exception(msg) {}
};

class InvalidArgumentException : public Exception {
public:
  InvalidArgumentException(const std::string& msg = "orang::InvalidArgumentException") : Exception(msg) {}
};

class OperationUnavailable : public Exception {
public:
  OperationUnavailable(const std::string& msg = "orang::OperationUnavailable") : Exception(msg) {}
};

}

#endif
