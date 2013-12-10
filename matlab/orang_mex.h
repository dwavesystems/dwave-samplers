#ifndef INCLUDED_ORANG_MEX_ORANGMEX_H
#define INCLUDED_ORANG_MEX_ORANGMEX_H

#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <vector>

#include <orang/orang.h>

// mexErrMsgAndId id values
extern const char* errIdInvalidArgument;
extern const char* errIdInternalError;
extern const char* errIdOutOfMemory;
extern const char* errIdExcessiveComplexity;

// matlab table structure field names
enum {
  TABLE_FIELD_VARS,
  TABLE_FIELD_DOMSIZES,
  TABLE_FIELD_VALUES,
  NUM_TABLE_FIELDS
};
extern const char* tableFieldNames[NUM_TABLE_FIELDS];

// matlab nodetables structure field names
enum {
  NODETABLES_FIELD_NODEVAR,
  NODETABLES_FIELD_SEPVARS,
  NODETABLES_FIELD_TABLES,
  NUM_NODETABLES_FIELDS
};
extern const char* nodeTablesFieldNames[NUM_NODETABLES_FIELDS];

// useful transform operators
struct Identity {
  template<typename T>
  T operator()(const T& t) const { return t; }
  template<typename T>
  T inv(const T& t) const { return t; }
};

struct AddOne {
  template<typename T>
  T operator()(const T& t) const { return t + T(1); }
  template<typename T>
  T inv(const T& t) const { return t - T(1); }
};

struct SubtractOne {
  template<typename T>
  T operator()(const T& t) const { return t - T(1); }
  template<typename T>
  T inv(const T& t) const { return t + T(1); }
};

class BadArray {
private:
  std::string what_;
public:
  BadArray(const std::string& what) : what_(what) {}
  const std::string& what() { return what_; }
};

// extract integral type value from scalar array
template<typename T>
T scalarArrayToIntegral(const mxArray* array) {
  if (mxIsComplex(array) || mxIsSparse(array) || !mxIsNumeric(array) || mxGetNumberOfElements(array) != 1) {
    throw BadArray("it must be a real scalar");
  }
  double value = mxGetScalar(array);
  double dMin = static_cast<double>(std::numeric_limits<T>::min());
  double dMax = static_cast<double>(std::numeric_limits<T>::max());
  if (value < dMin || value > dMax || value != floor(value)) {
    std::ostringstream msg;
    msg << "it must be an integer in the range [" << std::numeric_limits<T>::min()
        << "," << std::numeric_limits<T>::max() << "]";
    throw BadArray(msg.str());
  }
  return static_cast<T>(value);
}

// copy double mxArray to vector of some integral type
template<typename T, typename Op>
std::vector<T> doubleArrayToIntegralVector(const mxArray* array) {
  using std::vector;
  using std::size_t;
  using std::ostringstream;
  using std::numeric_limits;

  if (mxIsComplex(array) || !mxIsDouble(array) || mxIsSparse(array)) {
    throw BadArray("it must be a full, real vector of doubles");
  }

  const double* data = mxGetPr(array);
  size_t arraySize = mxGetNumberOfElements(array);
  vector<T> vec;
  vec.reserve(arraySize);

  Op op;
  double dMin = static_cast<double>(numeric_limits<T>::min());
  double dMax = static_cast<double>(numeric_limits<T>::max());

  for (const double* p = data; p != data + arraySize; ++p) {
    double x = op(*p);
    if (x < dMin || x > dMax || mxIsNaN(x) || x != floor(x)) {
      std::ostringstream msg;
      msg << "values must be integral and in the range [" << op.inv(dMin) << "," << op.inv(dMax) << "]";
      throw BadArray(msg.str());
    }
    vec.push_back(static_cast<T>(x));
  }

  return vec;
}

// initialize vector of tables.  no values are copied
template<typename Y>
std::vector<typename orang::Table<Y>::smartptr> initTableVector(const mxArray* tablesArray) {

  using std::floor;
  using std::numeric_limits;
  using std::ostringstream;
  using std::size_t;
  using std::vector;
  using orang::Var;
  using orang::VarVector;
  using orang::DomIndex;
  using orang::DomIndexVector;
  using orang::Table;

  // validate input
  if (!mxIsStruct(tablesArray)) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'tables' parameter is not a structure array");
  }

  int varsField = mxGetFieldNumber(tablesArray, tableFieldNames[TABLE_FIELD_VARS]);
  if (varsField < 0) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'tables' structure has no 'vars' field");
  }

  int domSizesField = mxGetFieldNumber(tablesArray, tableFieldNames[TABLE_FIELD_DOMSIZES]);
  if (domSizesField < 0) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'tables' structure has no 'domSizes' field");
  }

  size_t numTables = mxGetNumberOfElements(tablesArray);
  vector<typename Table<Y>::smartptr> tables;
  tables.reserve(numTables);

  for (size_t i = 0; i < numTables; ++i) {
    const mxArray* varsArray = mxGetFieldByNumber(tablesArray, i, varsField);
    const mxArray* domSizesArray = mxGetFieldByNumber(tablesArray, i, domSizesField);
    int fieldIndex = TABLE_FIELD_VARS;

    try {
      VarVector vars = doubleArrayToIntegralVector<Var,SubtractOne>(varsArray);
      fieldIndex = TABLE_FIELD_DOMSIZES;
      DomIndexVector domSizes = doubleArrayToIntegralVector<DomIndex, Identity>(domSizesArray);
      typename Table<Y>::smartptr tableptr( new Table<Y>(vars, domSizes) );
      tables.push_back(tableptr);
    } catch (BadArray& e) {
      mexErrMsgIdAndTxt(errIdInvalidArgument,
          "tables(%zu) has an invalid '%s' field: %s", i + 1, tableFieldNames[fieldIndex], e.what().c_str());
    } catch (orang::InvalidArgumentException& e) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "tables(%zu) is invalid: %s", i + 1, e.what().c_str());
    }
  }

  return tables;
}

std::vector<orang::Table<double>::smartptr> doubleTableVector(const mxArray* tablesArray);

template<typename Y>
mxArray* initTableMatlabArray(const std::vector<typename orang::Table<Y>::const_smartptr>& tables) {

  using std::size_t;
  using std::vector;
  using orang::TableVar;

  size_t numTables = tables.size();
  mxArray* tablesArray = mxCreateStructMatrix(1, numTables, NUM_TABLE_FIELDS, tableFieldNames);
  int varsField = mxGetFieldNumber(tablesArray, tableFieldNames[TABLE_FIELD_VARS]);
  int domSizesField = mxGetFieldNumber(tablesArray, tableFieldNames[TABLE_FIELD_DOMSIZES]);

  for (size_t i = 0; i < numTables; ++i) {
    const vector<TableVar>& vars = tables[i]->vars();
    size_t numVars = vars.size();
    mxArray* varsArray = mxCreateDoubleMatrix(1, numVars, mxREAL);
    mxArray* domSizesArray = mxCreateDoubleMatrix(1, numVars, mxREAL);

    double* varsData = mxGetPr(varsArray);
    double* domSizesData = mxGetPr(domSizesArray);
    for (size_t j = 0; j < numVars; ++j) {
      varsData[j] = vars[j].index + 1;
      domSizesData[j] = vars[j].domSize;
    }

    mxSetFieldByNumber(tablesArray, i, varsField, varsArray);
    mxSetFieldByNumber(tablesArray, i, domSizesField, domSizesArray);
  }

  return tablesArray;
}

template<typename T>
std::vector<typename orang::Table<T>::smartptr> doubleTableVector(const mxArray* tablesArray) {
  std::vector<typename orang::Table<T>::smartptr> tables = initTableVector<T>(tablesArray);

  int valuesField = mxGetFieldNumber(tablesArray, "values");
  if (valuesField < 0) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'tables' structure has no 'values' field");
  }

  size_t numTables = tables.size();
  for (size_t i = 0; i < numTables; ++i) {
    const mxArray* valuesArray = mxGetFieldByNumber(tablesArray, i, valuesField);
    if (mxIsComplex(valuesArray) || !mxIsDouble(valuesArray) || mxIsSparse(valuesArray)) {
      mexErrMsgIdAndTxt(errIdInvalidArgument,
          "tables(%zu) has an invalid 'values' field: it must be a full, real vector of doubles", i + 1);
    }

    size_t numValues = mxGetNumberOfElements(valuesArray);
    if (numValues != tables[i]->size()) {
      mexErrMsgIdAndTxt(errIdInvalidArgument,
          "tables(%zu).values size is wrong: it should be %zu", i + 1, tables[i]->size());
    }

    const double* valuesData = mxGetPr(valuesArray);
    copy(valuesData, valuesData + numValues, tables[i]->begin());
  }

  return tables;
}

template<typename Y>
mxArray* initNodeTablesMatlabArray(const std::vector<orang::NodeTables<Y> >& nodeTables) {

  using std::size_t;
  using std::transform;
  using std::vector;
  using orang::TableVar;
  using orang::NodeTables;

  size_t numNodeTables = nodeTables.size();
  mxArray* nodeTablesArray = mxCreateStructMatrix(1, numNodeTables, NUM_NODETABLES_FIELDS, nodeTablesFieldNames);
  int nodeVarField = mxGetFieldNumber(nodeTablesArray, nodeTablesFieldNames[NODETABLES_FIELD_NODEVAR]);
  int sepVarsField = mxGetFieldNumber(nodeTablesArray, nodeTablesFieldNames[NODETABLES_FIELD_SEPVARS]);
  int tablesField = mxGetFieldNumber(nodeTablesArray, nodeTablesFieldNames[NODETABLES_FIELD_TABLES]);

  for (size_t i = 0; i < numNodeTables; ++i) {
    mxArray* sepVarsArray = mxCreateDoubleMatrix(1, nodeTables[i].sepVars.size(), mxREAL);
    double* sepVarsData = mxGetPr(sepVarsArray);
    transform(nodeTables[i].sepVars.begin(), nodeTables[i].sepVars.end(), sepVarsData, AddOne());

    mxSetFieldByNumber(nodeTablesArray, i, nodeVarField, mxCreateDoubleScalar(nodeTables[i].nodeVar + 1));
    mxSetFieldByNumber(nodeTablesArray, i, sepVarsField, sepVarsArray);
    mxSetFieldByNumber(nodeTablesArray, i, tablesField, initTableMatlabArray<Y>(nodeTables[i].tables));
  }

  return nodeTablesArray;
}

mxArray* doubleNodeTablesMatlabArray(const std::vector<orang::NodeTables<double> >& nodeTables);

void validateVarOrder(const orang::VarVector& varOrder, orang::Var numVars);

#endif
