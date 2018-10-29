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
