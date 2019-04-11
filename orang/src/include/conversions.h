#ifndef ORANG_PY_CONVERSIONS_H_INCLUDED
#define ORANG_PY_CONVERSIONS_H_INCLUDED

#include <new>
#include <vector>

#include <cstdlib>

#include <base.h>
#include <table.h>

class MallocPtr {
private:
  void* p_;
  MallocPtr(const MallocPtr&);
  MallocPtr& operator=(const MallocPtr&);
public:
  MallocPtr(void* p = 0) : p_(p) {}
  ~MallocPtr() { std::free(p_); }
  void* get() { return p_; }
  void reset(void* p = 0) {
    if (p != p_) {
      std::free(p_);
      p_ = p;
    }
  }
  void* release() {
    void* p = p_;
    p_ = 0;
    return p;
  }
};

inline void* mallocOrThrow(std::size_t sz) {
  void* p = std::malloc(sz);
  if (!p) throw std::bad_alloc();
  return p;
}

std::vector<orang::Table<double>::smartptr> isingTables(
  int hLen, const double* hData,
  int jRows, int jCols, const double* jData,
  double beta
);

std::vector<orang::Table<double>::smartptr> quboTables(
  int qRows, int qCols, const double* qData,
  double beta
);

std::vector<orang::Table<double>::smartptr> cooTables(
    size_t numLinear,
    const double* lVals,
    size_t numQuadratic,
    const unsigned int* iRow, const unsigned int* iCol, const double* qVals,
    double low,
    double beta
);

orang::VarVector varOrderVec(int voLen, const int* voData, int numVars);

#endif
