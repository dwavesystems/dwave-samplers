#include <algorithm>
#include <limits>
#include <map>
#include <utility>
#include <vector>

#include <cstddef>
#include <cstdlib>

#include <boost/iterator/indirect_iterator.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#include <boost/random.hpp>
#include <boost/foreach.hpp>

#include <base.h>
#include <combine.h>
#include <table.h>
#include <task.h>
#include <treedecomp.h>
#include <buckettree.h>
#include <operations/logsumprod.h>

#include "python-api.h"
#include "conversions.h"

using std::size_t;
using std::copy;
using std::map;
using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::vector;

using boost::make_indirect_iterator;
using boost::posix_time::ptime;
using boost::posix_time::microsec_clock;
using boost::posix_time::min_date_time;

using orang::Var;
using orang::DomIndexVector;
using orang::VarVector;
using orang::Table;
using orang::TreeDecomp;
using orang::BucketTree;
using orang::TableMerger;

namespace {

typedef boost::variate_generator<boost::mt19937&, boost::uniform_01<> > Rng;
typedef orang::Task<orang::LogSumProductOperations<Rng> > SampleTask;
typedef std::vector<orang::Table<double>::smartptr> Tables;

typedef pair<Var, Var> VarPair;

struct PairMrgVals {
  double values[4];
};

typedef map<VarPair, PairMrgVals> PairMrgMap;

class Normalizer {
private:
  double logPf_;

public:
  Normalizer(double logPf) : logPf_(logPf) {}
  double operator()(double x) const { return exp(x - logPf_); }
};

unsigned int randomSeed(int userSeed) {
  if (userSeed >= 0) {
    return userSeed;
  } else {
    return static_cast<unsigned>((microsec_clock::local_time() - ptime(min_date_time)).total_microseconds());
  }
}

vector<double> singleMarginals(const BucketTree<SampleTask>& bucketTree) {
  vector<double> mrg(bucketTree.task().numVars());

  VarVector vars1(1);
  TableMerger<SampleTask> mergeTables(bucketTree.task());
  SampleTask::marginalizer_smartptr marginalizer = bucketTree.task().marginalizer();
  BOOST_FOREACH( const BucketTree<SampleTask>::nodetables_type& nt, bucketTree.nodeTables() ) {
    vars1[0] = nt.nodeVar;
    SampleTask::table_smartptr mTable = mergeTables(vars1, make_indirect_iterator(nt.tables.begin()),
        make_indirect_iterator(nt.tables.end()), *marginalizer);

    Normalizer normalize((*marginalizer)(0, *mTable));
    mrg[nt.nodeVar] = normalize((*mTable)[1]);
  }

  return mrg;
}

PairMrgMap pairMarginals(const BucketTree<SampleTask>& bucketTree) {
  PairMrgMap mrg;

  BOOST_FOREACH( const SampleTask::const_table_smartptr& t, bucketTree.task().tables() ) {
    if (t->vars().size() == 2) {
      VarPair p(t->vars()[0].index, t->vars()[1].index);
      mrg[p];
    }
  }

  VarVector vars2(2);
  TableMerger<SampleTask> mergeTables(bucketTree.task());
  SampleTask::marginalizer_smartptr marginalizer = bucketTree.task().marginalizer();
  BOOST_FOREACH( const BucketTree<SampleTask>::nodetables_type& nt, bucketTree.nodeTables() ) {
    BOOST_FOREACH( Var v, nt.sepVars ) {
      VarPair p(min(nt.nodeVar, v), max(nt.nodeVar, v));
      if (mrg.find(p) == mrg.end()) continue;
      vars2[0] = p.first;
      vars2[1] = p.second;
      SampleTask::table_smartptr mTable = mergeTables(vars2, make_indirect_iterator(nt.tables.begin()),
          make_indirect_iterator(nt.tables.end()), *marginalizer);

      Normalizer normalize((*marginalizer)(0, *mTable));
      PairMrgVals& mv = mrg[p];
      mv.values[0] = normalize((*mTable)[0]);
      mv.values[1] = normalize((*mTable)[1]);
      mv.values[2] = normalize((*mTable)[2]);
      mv.values[3] = normalize((*mTable)[3]);
    }
  }

  return mrg;
}

void sample(SampleTask& task, int z,
  int* voData, int voLen,
  double maxComplexity, int numSamples, bool marginals,
  double* logPf,
  int** samplesData, int* samplesRows, int* samplesCols,
  double** singleMrgData, int* singleMrgLen,
  double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
  int** pairData, int* pairRows, int* pairCols) {

  VarVector varOrder = varOrderVec(voLen, voData, task.numVars());
  TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
  if (!(decomp.complexity() <= maxComplexity)) throw std::runtime_error("complexity exceeded");

  bool solvable = numSamples > 0;
  BucketTree<SampleTask> bucketTree(task, decomp, DomIndexVector(task.numVars()), solvable, marginals);
  *logPf = bucketTree.problemValue();

  MallocPtr samplesMp;
  MallocPtr singleMrgMp;
  MallocPtr pairMrgMp;
  MallocPtr pairMp;

  if (solvable) {
    int numVars = static_cast<int>(task.numVars());
    *samplesRows = numSamples;
    *samplesCols = numVars;
    if (numeric_limits<size_t>::max() / sizeof(**samplesData) / *samplesRows / *samplesCols == 0) {
      throw std::invalid_argument("samples size too large");
    }
    samplesMp.reset(mallocOrThrow(static_cast<size_t>(*samplesRows) * *samplesCols * sizeof(**samplesData)));

    int s[2] = {z, 1};
    int* samplesMpData = static_cast<int*>(samplesMp.get());
    for (int i = 0; i < numSamples; ++i) {
      DomIndexVector samp = bucketTree.solve();
      for (int j = 0; j < numVars; ++j) {
        *samplesMpData++ = s[samp[j]];
      }
    }

  } else {
    *samplesRows = 0;
    *samplesCols = 0;
    samplesMp.reset(mallocOrThrow(1));
  }

  if (marginals) {
    vector<double> singleMrg = singleMarginals(bucketTree);
    *singleMrgLen = static_cast<int>(singleMrg.size());
    singleMrgMp.reset(mallocOrThrow(singleMrg.size() * sizeof(**singleMrgData)));
    copy(singleMrg.begin(), singleMrg.end(), static_cast<double*>(singleMrgMp.get()));

    PairMrgMap pairMrg = pairMarginals(bucketTree);
    *pairMrgRows = static_cast<int>(pairMrg.size());
    *pairMrgCols = 4;
    pairMrgMp.reset(mallocOrThrow(pairMrg.size() * 4 * sizeof(**pairMrgData)));
    *pairRows = static_cast<int>(pairMrg.size());
    *pairCols = 2;
    pairMp.reset(mallocOrThrow(pairMrg.size() * 2 * sizeof(**pairData)));
    int* pairMpData = static_cast<int*>(pairMp.get());
    double* pairMrgMpData = static_cast<double*>(pairMrgMp.get());
    BOOST_FOREACH( const PairMrgMap::value_type& e, pairMrg ) {
      *pairMpData++ = static_cast<int>(e.first.first);
      *pairMpData++ = static_cast<int>(e.first.second);
      *pairMrgMpData++ = e.second.values[0];
      *pairMrgMpData++ = e.second.values[1];
      *pairMrgMpData++ = e.second.values[2];
      *pairMrgMpData++ = e.second.values[3];
    }

  } else {
    *singleMrgLen = 0;
    singleMrgMp.reset(mallocOrThrow(1));

    *pairMrgRows = 0;
    *pairMrgCols = 4;
    pairMrgMp.reset(mallocOrThrow(1));

    *pairRows = 0;
    *pairCols = 2;
    pairMp.reset(mallocOrThrow(1));
  }

  *samplesData = static_cast<int*>(samplesMp.release());
  *singleMrgData = static_cast<double*>(singleMrgMp.release());
  *pairMrgData = static_cast<double*>(pairMrgMp.release());
  *pairData = static_cast<int*>(pairMp.release());
}

} // namespace {anonymous}

void sample_ising(
  double* hData, int hLen,
  double* jData, int jRows, int jCols,
  int* voData, int voLen,
  double maxComplexity, int numSamples, bool marginals, double beta, int rngSeed,
  double* logPf,
  int** samplesData, int* samplesRows, int* samplesCols,
  double** singleMrgData, int* singleMrgLen,
  double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
  int** pairData, int* pairRows, int* pairCols) {

  boost::mt19937 rngEngine(randomSeed(rngSeed));
  Rng rng(rngEngine, boost::uniform_01<>());

  int minVars = max(hLen, max(jRows, jCols));
  vector<Table<double>::smartptr> tables = isingTables(hLen, hData, jRows, jCols, jData, beta);
  SampleTask task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), rng, minVars);

  sample(task, -1, voData, voLen, maxComplexity, numSamples, marginals,
    logPf, samplesData, samplesRows, samplesCols, singleMrgData, singleMrgLen,
    pairMrgData, pairMrgRows, pairMrgCols, pairData, pairRows, pairCols);
}

void sample_qubo(
  double* qData, int qRows, int qCols,
  int* voData, int voLen,
  double maxComplexity, int numSamples, bool marginals, double beta, int rngSeed,
  double* logPf,
  int** samplesData, int* samplesRows, int* samplesCols,
  double** singleMrgData, int* singleMrgLen,
  double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
  int** pairData, int* pairRows, int* pairCols) {

  boost::mt19937 rngEngine(randomSeed(rngSeed));
  Rng rng(rngEngine, boost::uniform_01<>());

  int minVars = max(qRows, qCols);
  vector<Table<double>::smartptr> tables = quboTables(qRows, qCols, qData, beta);
  SampleTask task(make_indirect_iterator(tables.begin()), make_indirect_iterator(tables.end()), rng, minVars);

  sample(task, 0, voData, voLen, maxComplexity, numSamples, marginals,
    logPf, samplesData, samplesRows, samplesCols, singleMrgData, singleMrgLen,
    pairMrgData, pairMrgRows, pairMrgCols, pairData, pairRows, pairCols);
}

void sampleTables(
    Tables tables, int numVars,
    int low,  // -1 for SPIN, 0 for BINARY
    int* voData, int voLen, double maxComplexity,  // elimination order
    int numSamples,
    bool marginals, // whether to compute the marginals or not
    int rngSeed,

    double* logPf,
    int** samplesData, int* samplesRows, int* samplesCols,
    double** singleMrgData, int* singleMrgLen,
    double** pairMrgData, int* pairMrgRows, int* pairMrgCols,
    int** pairData, int* pairRows, int* pairCols
    ) {

    boost::mt19937 rngEngine(randomSeed(rngSeed));
    Rng rng(rngEngine, boost::uniform_01<>());

    SampleTask task(make_indirect_iterator(tables.begin()),
                    make_indirect_iterator(tables.end()),
                    rng,
                    numVars);

    sample(
        task,
        low,  // -1 for SPIN, 0 for BINARY
        voData, voLen, maxComplexity,
        numSamples,
        marginals,
        logPf,
        samplesData, samplesRows, samplesCols,
        singleMrgData, singleMrgLen,
        pairMrgData, pairMrgRows, pairMrgCols,
        pairData, pairRows, pairCols);
}
