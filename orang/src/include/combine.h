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
