#include <cmath>
#include <ctime>
#include <algorithm>
#include <functional>
#include <limits>
#include <new>
#include <vector>

#include <boost/iterator/indirect_iterator.hpp>
#include <boost/random.hpp>
#include <boost/cstdint.hpp>

#include <mex.h>

#include <orang/orang.h>

#include "orang_mex.h"

using std::bad_alloc;
using std::transform;
using std::vector;
using std::numeric_limits;
using std::floor;
using std::exp;

using boost::make_indirect_iterator;
using boost::uint32_t;

using orang::Var;
using orang::VarVector;
using orang::DomIndex;
using orang::DomIndexVector;
using orang::TreeDecomp;
using orang::Table;
using orang::Task;
using orang::Plus;
using orang::LogSumProductOperations;
using orang::BucketTree;
using orang::TableMerger;
using orang::LogSumMarginalizer;
using orang::Exception;

typedef boost::variate_generator<boost::mt19937&, boost::uniform_01<> > Rng;
typedef Task<LogSumProductOperations<Rng> > task_type;
typedef LogSumProductOperations<Rng>::CtorArgs logsumprod_ops_ctorargs;

enum {
  PARAM_TABLES,
  PARAM_VARORDER,
  PARAM_MAXCOMPLEXITY,
  PARAM_NUMSAMPLES,
  PARAM_X0,
  PARAM_MINVARS,
  PARAM_RNGSEED,
  NUM_PARAMS,
  NUM_REQUIRED_PARAMS = PARAM_MAXCOMPLEXITY + 1,
};

class Normalizer {
private:
  double logPf_;

public:
  Normalizer(double logPf) : logPf_(logPf) {}
  double operator()(double x) const { return exp(x - logPf_); }
};

mxArray* createMarginals(const BucketTree<task_type>& bucketTree) {

  static const char* marginalsFields[] = { "vars", "values" };

  size_t numMarginals = 0;
  BOOST_FOREACH( const BucketTree<task_type>::nodetables_type& nt, bucketTree.nodeTables() ) {
    numMarginals += nt.sepVars.size() + 1;
  }

  mxArray* marginalsArray = mxCreateStructMatrix(1, numMarginals, 2, marginalsFields);
  int varsField = mxGetFieldNumber(marginalsArray, "vars");
  int valuesField = mxGetFieldNumber(marginalsArray, "values");

  size_t i = 0;
  VarVector vars1(1);
  VarVector vars2(2);
  TableMerger<task_type> mergeTables(bucketTree.task());
  task_type::marginalizer_smartptr marginalizer = bucketTree.task().marginalizer();
  BOOST_FOREACH( const BucketTree<task_type>::nodetables_type& nt, bucketTree.nodeTables() ) {
    // single-variable marginal
    {
      vars1[0] = nt.nodeVar;
      task_type::table_smartptr mTable = mergeTables(vars1, make_indirect_iterator(nt.tables.begin()),
          make_indirect_iterator(nt.tables.end()), *marginalizer);

      mxArray* varsArray = mxCreateDoubleScalar(nt.nodeVar + 1);
      mxSetFieldByNumber(marginalsArray, i, varsField, varsArray);
      mxArray* valuesArray = mxCreateDoubleMatrix(2, 1, mxREAL);
      Normalizer normalize((*marginalizer)(0, *mTable));
      transform(mTable->begin(), mTable->end(), mxGetPr(valuesArray), normalize);
      mxSetFieldByNumber(marginalsArray, i, valuesField, valuesArray);
      ++i;
    }

    // pair marginals
    BOOST_FOREACH( Var v, nt.sepVars ) {
      if (v < nt.nodeVar) {
        vars2[0] = v;
        vars2[1] = nt.nodeVar;
      } else {
        vars2[0] = nt.nodeVar;
        vars2[1] = v;
      }
      task_type::table_smartptr mTable = mergeTables(vars2, make_indirect_iterator(nt.tables.begin()),
          make_indirect_iterator(nt.tables.end()), *marginalizer);

      mxArray* varsArray = mxCreateDoubleMatrix(1, 2, mxREAL);
      double* varsData = mxGetPr(varsArray);
      varsData[0] = vars2[0] + 1;
      varsData[1] = vars2[1] + 1;
      mxSetFieldByNumber(marginalsArray, i, varsField, varsArray);
      mxArray* valuesArray = mxCreateDoubleMatrix(2, 2, mxREAL);
      Normalizer normalize((*marginalizer)(0, *mTable));
      transform(mTable->begin(), mTable->end(), mxGetPr(valuesArray), normalize);
      mxSetFieldByNumber(marginalsArray, i, valuesField, valuesArray);
      ++i;
    }
  }

  return marginalsArray;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  static boost::mt19937 defaultRngEngine(std::time(0));
  static Rng defaultRng(defaultRngEngine, boost::uniform_01<>());
  static boost::mt19937 seededRngEngine;
  static Rng seededRng(seededRngEngine, boost::uniform_01<>());


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

    // number of samples
    badParam = "numSamples";
    size_t numSamples = 1;
    if (nrhs > PARAM_NUMSAMPLES && !mxIsEmpty(prhs[PARAM_NUMSAMPLES])) {
      numSamples = scalarArrayToIntegral<size_t>(prhs[PARAM_NUMSAMPLES]);
    }
    bool solvable = numSamples > 0;

    // min number of variables
    badParam = "minVars";
    Var minVars = 0;
    if (nrhs > PARAM_MINVARS && !mxIsEmpty(prhs[PARAM_MINVARS])) {
      minVars = scalarArrayToIntegral<Var>(prhs[PARAM_MINVARS]);
    }

    // return nodetables?
    bool hasNodeTables = nlhs > 2 || (nlhs > 1 && !solvable);

    // check number of outputs
    int maxOutputs = 1 + (solvable ? 1 : 0) + (hasNodeTables ? 1 : 0);
    if (nlhs > maxOutputs) {
      mexErrMsgIdAndTxt(errIdInvalidArgument, "Too many (>%d) output arguments given", maxOutputs);
    }

    // given rng seed?
    Rng* rng = &defaultRng;
    if (nrhs > PARAM_RNGSEED && !mxIsEmpty(prhs[PARAM_RNGSEED])) {
      uint32_t seed = scalarArrayToIntegral<uint32_t>(prhs[PARAM_RNGSEED]);
      seededRngEngine.seed(seed);
      rng = &seededRng;
    }

    vector<Table<double>::smartptr> tables = doubleTableVector<double>(prhs[PARAM_TABLES]);
    task_type task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), *rng, minVars);

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

    plhs[0] = mxCreateDoubleScalar(bucketTree.problemValue());

    if (solvable) {
      size_t numVars = task.numVars();
      plhs[1] = mxCreateDoubleMatrix(numVars, numSamples, mxREAL);
      double* data = mxGetPr(plhs[1]);

      for (size_t i = 0; i < numSamples; ++i) {
        DomIndexVector s = bucketTree.solve();
        transform(s.begin(), s.end(), data, AddOne());
        data += numVars;
      }
    }

    if (hasNodeTables) {
      mxArray*& ntRet = solvable ? plhs[2] : plhs[1];
      ntRet = createMarginals(bucketTree);
    }

  } catch (BadArray& e) {
    mexErrMsgIdAndTxt(errIdInvalidArgument, "'%s' parameter is invalid: %s", badParam, e.what().c_str());
  } catch (bad_alloc&) {
    mexErrMsgIdAndTxt(errIdOutOfMemory, "Out of memory");
  } catch (Exception& e) {
    mexErrMsgIdAndTxt(errIdInternalError, e.what().c_str());
  }


}
