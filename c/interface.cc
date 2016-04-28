
#include <vector>
#include <algorithm>
#include <stdio.h>
#include <ctime>
#include <boost/random.hpp>
#include <boost/iterator/indirect_iterator.hpp>


#include "orang/orang.h"
#include "orang/varorder.h"
#include "interface.h"
#include "helper.h"
#include "errors.h"
#include "mempool.h"

using std::vector;
using std::transform;

using orang::Exception;
using orang::Table;
using orang::VarVector;
using orang::Task;
using orang::DummyOperations;
using orang::greedyvarorder::Heuristics;
using orang::greedyVarOrder;
using orang::BucketTree;
using orang::TreeDecomp;
using orang::Var;
using orang::MinOperations;
using orang::Plus;
using orang::MinSolution;
using orang::MinSolutionSet;
using orang::LogSumProductOperations;
using orang::TableMerger;
using orang::LogSumMarginalizer;

typedef boost::uniform_01<boost::mt19937> varOrder_Rng;
typedef boost::variate_generator<boost::mt19937&, boost::uniform_01<> > sample_Rng;

typedef Task<DummyOperations> varOrder_task_type;
typedef Task<MinOperations<double, Plus<double> > > optimize_task_type;
typedef Task<LogSumProductOperations<sample_Rng> > sample_task_type;

typedef DummyOperations::CtorArgs ctor_args;
typedef LogSumProductOperations<sample_Rng>::CtorArgs logsumprod_ops_ctorargs;

const int MAX_ERROR_LENGTH = 200;

void validateVarOrder(const VarVector& varOrder, Var numVars) {
    vector<char> seen(numVars, 0);
    BOOST_FOREACH( Var v, varOrder ) {
        if (v >= numVars) {
            char errMsg[MAX_ERROR_LENGTH];
            snprintf(errMsg, MAX_ERROR_LENGTH, "Invalid variable elimination order: it contains %u but there are only %u variables", v + 1, numVars);
            throw Errors(errMsg);
        }
        if (seen[v]) {
            char errMsg[MAX_ERROR_LENGTH];
            snprintf(errMsg, MAX_ERROR_LENGTH, "Invalid variable elimination order: variable %u appears more than once", v + 1);
            throw Errors(errMsg);
        }
        seen[v] = 1;
    }
}

class Normalizer {
private:
  double logPf_;

public:
  Normalizer(double logPf) : logPf_(logPf) {}
  double operator()(double x) const { return exp(x - logPf_); }
};

int createMarginals(const BucketTree<sample_task_type>& bucketTree, Marginal ** marginals,
                    MemPool & mempool){
    size_t numMarginals = 0;
    BOOST_FOREACH( const BucketTree<sample_task_type>::nodetables_type& nt, bucketTree.nodeTables() ) {
        numMarginals += nt.sepVars.size() + 1;
    }
    *marginals = (Marginal *)mempool.malloc(numMarginals * sizeof(Marginal));
    Marginal * m = *marginals;
    size_t i = 0;
    VarVector vars1(1);
    VarVector vars2(2);
    TableMerger<sample_task_type> mergeTables(bucketTree.task());
    sample_task_type::marginalizer_smartptr marginalizer = bucketTree.task().marginalizer();
    BOOST_FOREACH( const BucketTree<sample_task_type>::nodetables_type& nt, bucketTree.nodeTables() ) {
    // unary marginals
    {
      vars1[0] = nt.nodeVar;
      sample_task_type::table_smartptr mTable = mergeTables(vars1, make_indirect_iterator(nt.tables.begin()),
          make_indirect_iterator(nt.tables.end()), *marginalizer);
      m[i].vars_len = 1;
      m[i].vars = (int *)mempool.malloc(1 * sizeof(int));
      m[i].vars[0] = nt.nodeVar;
      m[i].values = (double *)mempool.malloc(mTable->size() * sizeof(double));
      Normalizer normalize((*marginalizer)(0, *mTable));
      transform(mTable->begin(), mTable->end(), m[i].values, normalize);
      ++i;
    }
    // pairwise m
    BOOST_FOREACH( Var v, nt.sepVars ) {
      if (v < nt.nodeVar) {
        vars2[0] = v;
        vars2[1] = nt.nodeVar;
      } else {
        vars2[0] = nt.nodeVar;
        vars2[1] = v;
      }
      sample_task_type::table_smartptr mTable = mergeTables(vars2, make_indirect_iterator(nt.tables.begin()),
          make_indirect_iterator(nt.tables.end()), *marginalizer);
      m[i].vars_len = 2;
      m[i].vars = (int *)mempool.malloc(2 * sizeof(int));
      m[i].vars[0] = vars2[0];
      m[i].vars[1] = vars2[1];
      m[i].values = (double *)mempool.malloc(mTable->size() * sizeof(double));

      Normalizer normalize((*marginalizer)(0, *mTable));
      transform(mTable->begin(), mTable->end(), m[i].values, normalize);
      ++i;
    }
  }
  return numMarginals;
}

extern "C"{


int greedyVarOrder(TableEntry * tables, int tables_len, int maxComplexity,
                   int * clampRanks, int clampRanks_len,
                   int heuristic, float selectionScale,
                   int ** variableOrder, int * variableOrder_len, char * errorMessage){
    try{
        if (variableOrder == NULL)
            throw Errors("variableOrder is NULL");
        if (variableOrder_len == NULL)
            throw Errors("variableOrder_len is NULL");
        *variableOrder = NULL;
        static varOrder_Rng rng(boost::mt19937(static_cast<unsigned int>(std::time(0))));
        vector<Table<char>::smartptr> tb = createTables<char>(tables, tables_len);
        vector<int> cr;
        cr.reserve(clampRanks_len);
        for (int i = 0; i < clampRanks_len; i++)
            cr.push_back(clampRanks[i]);
        varOrder_task_type task(make_indirect_iterator(tb.begin()), make_indirect_iterator(tb.end()),
            ctor_args(), static_cast<Var>(cr.size()));
        if (clampRanks_len == 0)
            cr.assign(task.numVars(), 0);
        if (cr.size() != task.numVars())
            throw Errors("'clampRanks' parameter must be empty or have size no less than the largest variable index");
        Heuristics h;
        switch (heuristic){
            case HEURISTIC_MIN_DEG:
                h = orang::greedyvarorder::MIN_DEGREE;
                break;
            case HEURISTIC_W_MIN_DEG:
                h = orang::greedyvarorder::WEIGHTED_MIN_DEGREE;
                break;
            case HEURISTIC_MIN_FILL:
                h = orang::greedyvarorder::MIN_FILL;
                break;
            case HEURISTIC_W_MIN_FILL:
                h = orang::greedyvarorder::WEIGHTED_MIN_FILL;
                break;
            default:
                throw Errors("Invalid heuristic");
        };
        VarVector varOrder = greedyVarOrder(task, maxComplexity, cr, h, rng, selectionScale);
        *variableOrder_len = varOrder.size();
        // No mempool call is required as there is no later code that may throw.
        *variableOrder = (int *)malloc(varOrder.size() * sizeof(int));
        if (*variableOrder == NULL)
            throw std::bad_alloc();
        for (unsigned int i = 0; i < varOrder.size(); i++)
            (*variableOrder)[i] = varOrder[i];
        return 0;
    }
    catch(Errors & e){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "%s", e.what());
        return 1;
    }
    catch (std::bad_alloc &){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "Out of memory");
        return 1;
    }
    catch (Exception & e){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "%s", e.what().c_str());
        return 1;
    }
}

int optimize(TableEntry * tables, int tables_len, int * variableOrder, int variableOrder_len,
             int maxComplexity, int maxSolutions, int * initState, int initState_len, int minVars,
             double ** energies, int * energies_len, int ** states, int * varNum, char * errorMessage){
    char errMsg[MAX_ERROR_LENGTH];
    MemPool mempool;
    try{
        if (energies == NULL)
            throw Errors("energies is NULL");
        if (energies_len == NULL)
            throw Errors("energies_len is NULL");
        if (states == NULL)
            throw Errors("states is NULL");
        if (varNum == NULL)
            throw Errors("varNum is NULL");
        *energies = NULL;
        *states = NULL;
        vector<Table<double>::smartptr> tb = createTables<double>(tables, tables_len);
        bool solvable = maxSolutions > 0;
        optimize_task_type task(make_indirect_iterator(tb.begin()), make_indirect_iterator(tb.end()), 1, (Var)minVars);
        VarVector varOrder;
        varOrder.reserve(variableOrder_len);
        for (int i = 0; i < variableOrder_len; i++){
            varOrder.push_back(variableOrder[i]);
        }
        validateVarOrder(varOrder, task.numVars());
        DomIndexVector x0(task.numVars());
        if (initState_len > 0 && initState_len != (int)task.numVars()){
            snprintf(errMsg, MAX_ERROR_LENGTH, "'x0' parameter must have %d variables", task.numVars());
            throw Errors(errMsg);
        }
        for (int i = 0; i < initState_len; i++){
            x0[i] = initState[i];
            if (x0[i] > task.domSize(i)) {
                snprintf(errMsg, MAX_ERROR_LENGTH, "x0(%u) is invalid: domain size of variable is %d", i + 1, task.domSize(i));
                throw Errors(errMsg);
            }
        }
        TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
        if (decomp.complexity() > maxComplexity) {
            snprintf(errMsg, MAX_ERROR_LENGTH, "Tree decomposition complexity is too high (%f)", decomp.complexity());
            throw Errors(errMsg);
        }
        BucketTree<optimize_task_type> bucketTree(task, decomp, x0, solvable, false);
        double baseValue = bucketTree.problemValue();
        if (!solvable){
            *energies_len = 1;
            *energies = (double *)malloc(1 * sizeof(double));
            (*energies)[0] = baseValue;
            return 0;
        }
        task.maxSolutions(maxSolutions);
        MinSolutionSet<double> solutionSet = bucketTree.solve();

        size_t numSolutions = solutionSet.solutions().size();
        size_t numVars = task.numVars();

        *energies_len = numSolutions;
        *energies = (double *)mempool.malloc(numSolutions * sizeof(double));

        double* valueOffsetData = *energies;
        BOOST_FOREACH( const MinSolution<double>& s, solutionSet.solutions() ) {
            *valueOffsetData++ = baseValue + s.value;
        }
        *varNum = numVars;
        *states = (int *)mempool.malloc(numVars * numSolutions * sizeof(int));
        int * st = *states;
        BOOST_FOREACH( const MinSolution<double>& s, solutionSet.solutions() ) {
            transform(s.solution.begin(), s.solution.end(), st, doNothing<double, int>);
            st += numVars;
        }
        mempool.release();
        return 0;
    }
    catch(Errors & e){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "%s", e.what());
        return 1;
    }
    catch (std::bad_alloc &){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "Out of memory");
        return 1;
    }
    catch (Exception & e){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "%s", e.what().c_str());
        return 1;
    }
}

int sample(TableEntry * tables, int tables_len, int * variableOrder, int variableOrder_len,
             int maxComplexity, int sampleNum, int * initState, int initState_len, int minVars,
             int seed, int returnMarginals,
             double * logZ, int ** samples, int * varNum,
             Marginal ** marginals, int * marginals_len,
             char * errorMessage){
    char errMsg[MAX_ERROR_LENGTH];
    MemPool mempool;
    static boost::mt19937 defaultRngEngine(std::time(0));
    static sample_Rng defaultRng(defaultRngEngine, boost::uniform_01<>());
    static boost::mt19937 seededRngEngine;
    static sample_Rng seededRng(seededRngEngine, boost::uniform_01<>());
    try{
        if (logZ == NULL)
            throw Errors("logZ is NULL");
        if (samples == NULL)
            throw Errors("samples is NULL");
        if (varNum == NULL)
            throw Errors("varNum is NULL");
        if (returnMarginals && marginals == NULL)
            throw Errors("marginals is NULL");
        if (returnMarginals && marginals_len == NULL)
            throw Errors("marginals_len is NULL");
        *samples = NULL;
        *marginals = NULL;
        sample_Rng * rng = &defaultRng;
        if (seed >= 0){
            seededRngEngine.seed((unsigned int)seed);
            rng = &seededRng;
        }
        vector<Table<double>::smartptr> tb = createTables<double>(tables, tables_len);
        bool solvable = sampleNum > 0;
        sample_task_type task(make_indirect_iterator(tb.begin()), make_indirect_iterator(tb.end()), *rng, (Var)minVars);
        VarVector varOrder;
        varOrder.reserve(variableOrder_len);
        for (int i = 0; i < variableOrder_len; i++){
            varOrder.push_back(variableOrder[i]);
        }
        validateVarOrder(varOrder, task.numVars());
        DomIndexVector x0(task.numVars());
        if (initState_len > 0 && initState_len != (int)task.numVars()){
            snprintf(errMsg, MAX_ERROR_LENGTH, "'x0' parameter must have %d variables", task.numVars());
            throw Errors(errMsg);
        }
        for (int i = 0; i < initState_len; i++){
            x0[i] = initState[i];
            if (x0[i] > task.domSize(i)) {
                snprintf(errMsg, MAX_ERROR_LENGTH, "x0(%u) is invalid: domain size of variable is %d", i + 1, task.domSize(i));
                throw Errors(errMsg);
            }
        }
        TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
        if (decomp.complexity() > maxComplexity) {
            snprintf(errMsg, MAX_ERROR_LENGTH, "Tree decomposition complexity is too high (%f)", decomp.complexity());
            throw Errors(errMsg);
        }
        BucketTree<sample_task_type> bucketTree(task, decomp, x0, solvable, returnMarginals);
        *logZ = bucketTree.problemValue();
        if (solvable){
            size_t numVars = task.numVars();
            *varNum = numVars;
            *samples = (int *)mempool.malloc(numVars * sampleNum * sizeof(int));
            int * st = *samples;
            for (int i = 0; i < sampleNum; ++i){
                DomIndexVector s = bucketTree.solve();
                transform(s.begin(), s.end(), st, doNothing<double, int>);
                st += numVars;
            }
        }
        if (returnMarginals){
            *marginals_len = createMarginals(bucketTree, marginals, mempool);
        }
        mempool.release();
        return 0;
    }
    catch(Errors & e){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "%s", e.what());
        return 1;
    }
    catch (std::bad_alloc &){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "Out of memory");
        return 1;
    }
    catch (Exception & e){
        snprintf(errorMessage, MAX_ERROR_LENGTH, "%s", e.what().c_str());
        return 1;
    }
}

void free_varOrder(int * varOrder){
    free(varOrder);
}

void free_energies(double * energies){
    free(energies);
}

void free_states(int * states){
    free(states);
}

void free_marginals(Marginal * marginals, int marginals_len){
    if (marginals == NULL || marginals_len == 0)
        return;
    for (int i = 0; i < marginals_len; i++){
        free(marginals[i].vars);
        free(marginals[i].values);
    }
    free(marginals);
}

} // extern "C"
