#include <cstddef>
#include <algorithm>
#include <vector>

#include <boost/foreach.hpp>

#include <mex.h>

#include <orang/orang.h>

#include "orang_mex.h"

using std::size_t;
using std::copy;
using std::vector;

using orang::Var;
using orang::VarVector;
using orang::Table;
using orang::NodeTables;

const char* errIdInvalidArgument = "orang_mex:invalid_argument";
const char* errIdOutOfMemory = "orang_mex:out_of_memory";
const char* errIdInternalError = "orang_mex:internal_error";
const char* errIdExcessiveComplexity = "orang_mex:excessive_complexity";

const char* tableFieldNames[NUM_TABLE_FIELDS] = {
    "vars",
    "domSizes",
    "values"
};

const char* nodeTablesFieldNames[NUM_NODETABLES_FIELDS] = {
    "nodeVar",
    "sepVars",
    "tables"
};

mxArray* doubleNodeTablesMatlabArray(const vector<NodeTables<double> >& nodeTables) {

  mxArray* nodeTablesArray = initNodeTablesMatlabArray(nodeTables);

  size_t numNodeTables = nodeTables.size();
  int tablesField = mxGetFieldNumber(nodeTablesArray, nodeTablesFieldNames[NODETABLES_FIELD_TABLES]);

  for (size_t i = 0; i < numNodeTables; ++i) {
    mxArray* tablesArray = mxGetFieldByNumber(nodeTablesArray, i, tablesField);
    int valuesField = mxGetFieldNumber(tablesArray, tableFieldNames[TABLE_FIELD_VALUES]);
    size_t numTables = nodeTables[i].tables.size();

    for (size_t j = 0; j < numTables; ++j) {
      mxArray* valuesArray = mxCreateDoubleMatrix(nodeTables[i].tables[j]->size(), 1, mxREAL);
      double* valuesData = mxGetPr(valuesArray);
      copy(nodeTables[i].tables[j]->begin(), nodeTables[i].tables[j]->end(), valuesData);
      mxSetFieldByNumber(tablesArray, j, valuesField, valuesArray);
    }
  }

  return nodeTablesArray;
}

void validateVarOrder(const VarVector& varOrder, Var numVars) {
  vector<char> seen(numVars, 0);
  BOOST_FOREACH( Var v, varOrder ) {
    if (v >= numVars) {
      mexErrMsgIdAndTxt(errIdInvalidArgument,
          "Invalid variable elimination order: it contains %u but there are only %u variables", v + 1, numVars);
    }

    if (seen[v]) {
      mexErrMsgIdAndTxt(errIdInvalidArgument,
          "Invalid variable elimination order: variable %u appears more than once", v + 1);
    }

    seen[v] = 1;
  }
}



