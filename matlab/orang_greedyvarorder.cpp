#include <cctype>
#include <ctime>
#include <algorithm>
#include <map>
#include <vector>
#include <new>

#include <mex.h>

#include <boost/random.hpp>
#include <boost/iterator/indirect_iterator.hpp>

#include <orang/orang.h>
#include "orang_mex.h"

using std::transform;
using std::map;
using std::vector;
using std::bad_alloc;

using boost::make_indirect_iterator;

using orang::Exception;
using orang::Var;
using orang::VarVector;
using orang::Table;
using orang::DummyOperations;
using orang::Task;
using orang::greedyvarorder::Heuristics;
using orang::greedyVarOrder;

typedef boost::uniform_01<boost::mt19937> Rng;
typedef Task<DummyOperations> task_type;
typedef DummyOperations::CtorArgs ctor_args;

enum {
  PARAM_TABLES,
  PARAM_MAXCOMPLEXITY,
  PARAM_CLAMPRANKS,
  PARAM_HEURISTIC,
  PARAM_SELECTIONSCALE,
  NUM_PARAMS,
  NUM_REQUIRED_PARAMS = PARAM_HEURISTIC + 1,
};

struct StrCmpi {
  bool operator()(const char* x, const char* y) const {
    int xl = std::tolower(*x);
    int yl = std::tolower(*y);
    while (*x != '\0' && *y != '\0') {
      if (xl != yl) {
        return xl < yl;
      }
      ++x;
      ++y;
      xl = std::tolower(*x);
      yl = std::tolower(*y);
    }
    return xl < yl;
  }
};

typedef map<const char*, Heuristics, StrCmpi> hmap_type;
const hmap_type& getHeuristicMap() {
  static hmap_type hmap;
  if (hmap.empty()) {
    hmap["mindeg"] = orang::greedyvarorder::MIN_DEGREE;
    hmap["wmindeg"] = orang::greedyvarorder::WEIGHTED_MIN_DEGREE;
    hmap["minfill"] = orang::greedyvarorder::MIN_FILL;
    hmap["wminfill"] = orang::greedyvarorder::WEIGHTED_MIN_FILL;
  }
  return hmap;
}



void mexFunction(int, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  static Rng rng(boost::mt19937(static_cast<unsigned int>(std::time(0))));

  if (nrhs < NUM_REQUIRED_PARAMS || nrhs > NUM_PARAMS) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "Between %d and %d arguments required", NUM_REQUIRED_PARAMS, NUM_PARAMS);
  }

  // max complexity
  const mxArray* mcArray = prhs[PARAM_MAXCOMPLEXITY];
  if (mxIsComplex(mcArray) || !mxIsNumeric(mcArray) || mxGetNumberOfElements(mcArray) != 1) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'maxComplexity' parameter must be a real scalar");
  }
  double maxComplexity = mxGetScalar(mcArray);
  if (mxIsNaN(maxComplexity)) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'maxComplexity' parameter must be a number");
  }

  // heuristic
  if (!mxIsChar(prhs[PARAM_HEURISTIC])) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'heuristic' parameter must be a string");
  }
  char* heurStr = mxArrayToString(prhs[PARAM_HEURISTIC]);
  if (!heurStr) {
    mexErrMsgIdAndTxt(errIdOutOfMemory, "Out of memory");
  }
  const hmap_type& hmap = getHeuristicMap();
  hmap_type::const_iterator hIter = hmap.find(heurStr);
  mxFree(heurStr);
  if (hIter == hmap.end()) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "Invalid heuristic");
  }
  Heuristics heuristic = hIter->second;

  // selection scale
  float selectionScale = 1.0f;
  if (nrhs > PARAM_SELECTIONSCALE) {
    const mxArray* ssArray = prhs[PARAM_SELECTIONSCALE];
    if (mxIsComplex(ssArray) || !mxIsNumeric(ssArray) || mxGetNumberOfElements(ssArray) != 1) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "'selectionScale' parameter must be a real scalar");
    }
    selectionScale = static_cast<float>(mxGetScalar(ssArray));
    if (selectionScale < 0.0 || !mxIsFinite(selectionScale) || mxIsNaN(selectionScale)) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "'selectionScale' parameter must be non-negative and finite");
    }
  }

  try {

    // clampranks
    vector<int> clampRanks;
    if (!mxIsEmpty(prhs[PARAM_CLAMPRANKS])) {
      clampRanks = doubleArrayToIntegralVector<int,Identity>(prhs[PARAM_CLAMPRANKS]);
    }

    // tables
    vector<Table<char>::smartptr> tables = initTableVector<char>(prhs[PARAM_TABLES]);
    task_type task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()),
        ctor_args(), static_cast<Var>(clampRanks.size()));

    // default clampranks
    if (mxIsEmpty(prhs[PARAM_CLAMPRANKS])) {
      clampRanks.assign(task.numVars(), 0);
    }
    if (clampRanks.size() != task.numVars()) {
      mexErrMsgIdAndTxt(errIdInvalidArgument,
          "'clampRanks' parameter must be empty or have size no less than the largest variable index");
    }

    // go!
    VarVector varOrder = greedyVarOrder(task, maxComplexity, clampRanks, heuristic, rng, selectionScale);

    // return
    plhs[0] = mxCreateDoubleMatrix(1, varOrder.size(), mxREAL);
    double* varOrderData = mxGetPr(plhs[0]);
    transform(varOrder.begin(), varOrder.end(), varOrderData, AddOne());

  } catch (BadArray& e) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'clampRanks' parameter is invalid: %s", e.what().c_str());
  } catch (bad_alloc&) {
    mexErrMsgIdAndTxt(errIdOutOfMemory, "Out of memory");
  } catch (Exception& e) {
    mexErrMsgIdAndTxt(errIdInternalError, e.what().c_str());
  }
}
