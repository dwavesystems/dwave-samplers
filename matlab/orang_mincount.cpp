#include <cmath>
#include <algorithm>
#include <functional>
#include <limits>
#include <new>
#include <vector>

#include <boost/iterator/indirect_iterator.hpp>

#include <mex.h>

#include <orang/orang.h>

#include "orang_mex.h"

using std::bad_alloc;
using std::plus;
using std::bind1st;
using std::transform;
using std::vector;
using std::numeric_limits;
using std::floor;

using boost::make_indirect_iterator;

using orang::Var;
using orang::VarVector;
using orang::DomIndex;
using orang::DomIndexVector;
using orang::TreeDecomp;
using orang::Table;
using orang::Task;
using orang::CountOperations;
using orang::ValueCount;
using orang::BucketTree;
using orang::Exception;

typedef ValueCount<double> value_type;
typedef CountOperations<double> ops_type;
typedef Task<ops_type> task_type;

enum {
  PARAM_TABLES,
  PARAM_VARORDER,
  PARAM_MAXCOMPLEXITY,
  PARAM_RELEPS,
  PARAM_X0,
  NUM_PARAMS,
  NUM_REQUIRED_PARAMS = PARAM_MAXCOMPLEXITY + 1,
};

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  if (nrhs < NUM_REQUIRED_PARAMS || nrhs > NUM_PARAMS) {
    mexErrMsgIdAndTxt(errIdInvalidArgument,
        "Between %d and %d parameters required", NUM_REQUIRED_PARAMS, NUM_PARAMS);
  }

  if (nlhs > 2) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "Too many (>%d) output arguments given", 2);
  }

  // max complexity (of tree decomp)
  const mxArray* mcArray = prhs[PARAM_MAXCOMPLEXITY];
  if (mxIsComplex(mcArray) || !mxIsNumeric(mcArray) || mxGetNumberOfElements(mcArray) != 1) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'maxComplexity' parameter must be a real scalar");
  }
  double maxComplexity = mxGetScalar(mcArray);
  if (mxIsNaN(maxComplexity)) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'maxComplexity' parameter must be a number");
  }

  // tolerable relative error for comparisons
  double relEps = 0.0;
  if (nrhs > PARAM_RELEPS && !mxIsEmpty(prhs[PARAM_RELEPS])) {
    const mxArray* relepsArray = prhs[PARAM_RELEPS];
    if (mxIsComplex(relepsArray) || !mxIsNumeric(relepsArray) || mxGetNumberOfElements(relepsArray) != 1) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "'relEps' parameter must be a real scalar");
    }
    relEps = mxGetScalar(relepsArray);
    if (mxIsNaN(relEps)) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "'relEps' parameter must be a number");
    }
  }

  const char* badParam = 0;
  try {
    vector<Table<value_type>::smartptr> tables = doubleTableVector<value_type>(prhs[PARAM_TABLES]);
    task_type task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), relEps );

    badParam = "varOrder";
    VarVector varOrder = doubleArrayToIntegralVector<Var, SubtractOne>(prhs[PARAM_VARORDER]);
    validateVarOrder(varOrder, task.numVars());

    badParam = "x0";
    DomIndexVector x0(task.numVars());
    if (nrhs > PARAM_X0 && !mxIsEmpty(prhs[PARAM_X0])) {
      x0 = doubleArrayToIntegralVector<DomIndex, SubtractOne>(prhs[PARAM_X0]);
      if (x0.size() != task.numVars()) {
        mexErrMsgIdAndTxt(errIdInvalidArgument, "'x0' parameter must have %zu variables", task.numVars());
      }
      for (Var i = 0; i < task.numVars(); ++i) {
        if (x0[i] > task.domSize(i)) {
          mexErrMsgIdAndTxt(errIdInvalidArgument,
              "x0(%u) is invalid: domain size of variable is %zu", i + 1, task.domSize(i));
        }
      }
    }

    TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
    if (decomp.complexity() > maxComplexity) {
      mexErrMsgIdAndTxt(errIdExcessiveComplexity,
          "Tree decomposition complexity is too high (%f)", decomp.complexity());
    }

    BucketTree<task_type> bucketTree(task, decomp, x0, false, false);
    ValueCount<double> baseValue = bucketTree.problemValue();
    plhs[0] = mxCreateDoubleScalar(baseValue.count());
    plhs[1] = mxCreateDoubleScalar(baseValue.value());

  } catch (BadArray& e) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'%s' parameter is invalid: %s", badParam, e.what().c_str());
  } catch (bad_alloc&) {
    mexErrMsgIdAndTxt(errIdOutOfMemory, "Out of memory");
  } catch (Exception& e) {
    mexErrMsgIdAndTxt(errIdInternalError, e.what().c_str());
  }


}
