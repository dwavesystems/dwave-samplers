#include <algorithm>
#include <stdexcept>
#include <vector>

#include <orang/table.h>

#include "conversions.h"

using std::max;
using std::min;
using std::vector;

using orang::Table;
using orang::DomIndexVector;
using orang::VarVector;

vector<Table<double>::smartptr> isingTables(
  int hLen, const double* hVals,
  int jRows, int jCols, const double* jVals,
  double beta
) {

  static const DomIndexVector ds1(1, 2);
  static const DomIndexVector ds2(2, 2);

  VarVector vars1(1);
  VarVector vars2(2);
  vector<Table<double>::smartptr> tables;

  for (int i = 0; i < hLen; ++i) {
    if (hVals[i] != 0.0) {
      vars1[0] = i;
      Table<double>::smartptr t(new Table<double>(vars1, ds1));
      (*t)[0] = beta * hVals[i];
      (*t)[1] = -beta * hVals[i];
      tables.push_back(t);
    }
  }

  for (int i = 0; i < jRows; ++i) {
    for (int j = 0; j < jCols; ++j, ++jVals) {
      if (*jVals != 0.0) {
        if (i == j) throw std::invalid_argument("nonzero J entry");
        vars2[0] = min(i, j);
        vars2[1] = max(i, j);
        Table<double>::smartptr t(new Table<double>(vars2, ds2));
        (*t)[0] = -beta * *jVals;
        (*t)[1] = beta * *jVals;
        (*t)[2] = beta * *jVals;
        (*t)[3] = -beta * *jVals;
        tables.push_back(t);
      }
    }
  }

  return tables;
}

vector<Table<double>::smartptr> quboTables(
  int qRows, int qCols, const double* qVals,
  double beta
) {

  static const DomIndexVector ds1(1, 2);
  static const DomIndexVector ds2(2, 2);

  VarVector vars1(1);
  VarVector vars2(2);
  vector<Table<double>::smartptr> tables;

  for (int i = 0; i < qRows; ++i) {
    for (int j = 0; j < qCols; ++j, ++qVals) {
      if (*qVals != 0.0) {
        if (i == j) {
          vars1[0] = i;
          Table<double>::smartptr t(new Table<double>(vars1, ds1));
          (*t)[0] = 0.0;
          (*t)[1] = -beta * *qVals;
          tables.push_back(t);
        } else {
          vars2[0] = min(i, j);
          vars2[1] = max(i, j);
          Table<double>::smartptr t(new Table<double>(vars2, ds2));
          (*t)[0] = 0.0;
          (*t)[1] = 0.0;
          (*t)[2] = 0.0;
          (*t)[3] = -beta * *qVals;
          tables.push_back(t);
        }
      }
    }
  }

  return tables;
}

VarVector varOrderVec(int voLen, const int* voData, int numVars) {
  if (voLen < 0) throw std::invalid_argument("negative voLen");

  vector<char> seen(numVars, 0);
  VarVector varOrder;
  varOrder.reserve(voLen);

  for (int i = 0; i < voLen; ++i) {
    if (voData[i] < 0 || voData[i] >= numVars) throw std::invalid_argument("variable order out of range");
    if (seen[voData[i]]) throw std::invalid_argument("duplicate variable order entry");
    seen[voData[i]] = 1;
    varOrder.push_back(voData[i]);
  }

  return varOrder;
}
