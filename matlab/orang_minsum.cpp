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
using orang::Plus;
using orang::MinOperations;
using orang::MinSolution;
using orang::MinSolutionSet;
using orang::BucketTree;
using orang::Exception;

typedef Task<MinOperations<double, Plus<double> > > task_type;

enum {
  PARAM_TABLES,
  PARAM_VARORDER,
  PARAM_MAXCOMPLEXITY,
  PARAM_MAXSOLUTIONS,
  PARAM_X0,
  PARAM_MINVARS,
  NUM_PARAMS,
  NUM_REQUIRED_PARAMS = PARAM_MAXCOMPLEXITY + 1,
};

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  if (nrhs < NUM_REQUIRED_PARAMS || nrhs > NUM_PARAMS) {
    mexErrMsgIdAndTxt(errIdInvalidArgument,
        "Between %d and %d parameters required", NUM_REQUIRED_PARAMS, NUM_PARAMS);
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

  const char* badParam = 0;
  try {

    // min number of variables
    badParam = "minVars";
    Var minVars = 0;
    if (nrhs > PARAM_MINVARS && !mxIsEmpty(prhs[PARAM_MINVARS])) {
      minVars = scalarArrayToIntegral<Var>(prhs[PARAM_MINVARS]);
    }

    // maximum number of solutions
    badParam = "maxSolutions";
    size_t maxSolutions = 1;
    if (nrhs > PARAM_MAXSOLUTIONS && !mxIsEmpty(prhs[PARAM_MAXSOLUTIONS])) {
      maxSolutions = scalarArrayToIntegral<size_t>(prhs[PARAM_MAXSOLUTIONS]);
    }
    bool solvable = maxSolutions > 0;

    // return nodetables?
    bool hasNodeTables = nlhs > 2 || (nlhs > 1 && !solvable);

    // check number of outputs
    int maxOutputs = 1 + (solvable ? 1 : 0) + (hasNodeTables ? 1 : 0);
    if (nlhs > maxOutputs) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "Too many (>%d) output arguments given", maxOutputs);
    }

    vector<Table<double>::smartptr> tables = doubleTableVector<double>(prhs[PARAM_TABLES]);
    task_type task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), 1, minVars);

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

    BucketTree<task_type> bucketTree(task, decomp, x0, solvable, hasNodeTables);
    double baseValue = bucketTree.problemValue();

    if (!solvable) {
      plhs[0] = mxCreateDoubleScalar(baseValue);

    } else {
      task.maxSolutions(maxSolutions);
      MinSolutionSet<double> solutionSet = bucketTree.solve();

      size_t numSolutions = solutionSet.solutions().size();
      size_t numVars = task.numVars();

      plhs[0] = mxCreateDoubleMatrix(1, numSolutions, mxREAL);
      double* valueOffsetData = mxGetPr(plhs[0]);
      BOOST_FOREACH( const MinSolution<double>& s, solutionSet.solutions() ) {
        *valueOffsetData++ = baseValue + s.value;
      }

      if (nlhs > 1) {
        plhs[1] = mxCreateDoubleMatrix(numVars, numSolutions, mxREAL);
        double* solutionData = mxGetPr(plhs[1]);
        BOOST_FOREACH( const MinSolution<double>& s, solutionSet.solutions() ) {
          transform(s.solution.begin(), s.solution.end(), solutionData, AddOne());
          solutionData += numVars;
        }
      }
    }

    if (hasNodeTables) {
      mxArray*& ntRet = solvable ? plhs[2] : plhs[1];
      ntRet = doubleNodeTablesMatlabArray(bucketTree.nodeTables());
    }

  } catch (BadArray& e) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'%s' parameter is invalid: %s", badParam, e.what().c_str());
  } catch (bad_alloc&) {
    mexErrMsgIdAndTxt(errIdOutOfMemory, "Out of memory");
  } catch (Exception& e) {
    mexErrMsgIdAndTxt(errIdInternalError, e.what().c_str());
  }


}
