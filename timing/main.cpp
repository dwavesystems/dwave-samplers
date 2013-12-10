#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>

#include <boost/random.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/random.hpp>

#include <orang/orang.h>

#include "uai.h"
#include "settings.h"

using std::copy;
using std::ostream;
using std::ostream_iterator;
using std::cin;
using std::cout;
using std::cerr;
using std::vector;

using boost::uniform_01;
using boost::mt19937;
using boost::make_indirect_iterator;

using orang::TableVar;
using orang::Table;
using orang::Multiply;
using orang::MinOperations;
using orang::Task;
using orang::BucketTree;
using orang::greedyVarOrder;
using orang::VarVector;
using orang::varIndex;
using orang::TreeDecomp;

namespace orang {
ostream& operator<<(ostream& out, const TableVar& var) {
  out << "[index=" << var.index << ",domSize=" << var.domSize << ",stepSize=" << var.stepSize << "]";
  return out;
}

ostream& operator<<(ostream& out, const Table<double>& table) {
  out << "vars: ";
  copy(table.vars().begin(), table.vars().end(), ostream_iterator<TableVar>(out, " "));
  out << "values: ";
  copy(table.begin(), table.end(), ostream_iterator<double>(out, " "));
  return out;
}

ostream& operator<<(ostream& out, const Table<double>::smartptr& tablePtr) {
  out << *tablePtr;
  return out;
}

}

typedef Task<MinOperations<double, Multiply<double> > > TaskType;

typedef boost::uniform_01<boost::mt19937> Rng;
int seed = 0;

int main() {

  limitMemory();
  Rng rng((boost::mt19937(seed)));

  try {
    ParsedProblem parsed = parseUaiProblem(cin);
    ParsedEvidence evidence = parseUaiEvidence(cin, parsed);

    if (evidence.size() > 1) {
      cerr << "Too many evidence sets\n";
      return 1;
    }

//    cout << "numVars: " << parsed.domSizes.size() << "\n";
//    cout << "tables:\n";
//    copy(parsed.tables.begin(), parsed.tables.end(), ostream_iterator<Table<double>::smartptr>(cout, "\n"));
//
//    cout << "Evidence:\n";
//    BOOST_FOREACH( const EvidenceMap& e, evidence ) {
//      BOOST_FOREACH( const EvidenceMap::value_type& v, e ) {
//        cout << " x" << v.first << "=" << v.second;
//      }
//      cout << "\n";
//    }

    TaskType task(
        make_indirect_iterator(parsed.tables.begin()),
        make_indirect_iterator(parsed.tables.end()),
        1,
        parsed.domSizes.size());

    vector<int> clampRanks(task.numVars());
    VarVector x(task.numVars());
    size_t numClamped = 0;
    if (evidence.size() == 1) {
      BOOST_FOREACH( const Evidence::value_type& e, evidence[0] ) {
        if (clampRanks[e.first] >= 0) {
          ++numClamped;
          clampRanks[e.first] = -1;
          x[e.first] = e.second;
        }
      }
    }

    bool foundExact = false;
    VarVector varOrder;
    for (int i = 0; !foundExact && i < settings::maxExactAttempts; ++i) {
      varOrder = greedyVarOrder(task, settings::maxExactComplexity, clampRanks,
          orang::greedyvarorder::WEIGHTED_MIN_DEGREE, rng, settings::exactScaling);
      foundExact = varOrder.size() == task.numVars() - numClamped;
    }

    if (foundExact) {
      cout << "Found exact elimination order\n";
      TreeDecomp decomp(task.graph(), varOrder, task.domSizes());
      BucketTree<TaskType> bucketTree(task, decomp, VarVector(task.numVars()), true, false);
      TaskType::solution_type x = bucketTree.solve();
      cout << bucketTree.problemValue() << "   ";
      copy(x.solutions().begin()->solution.begin(), x.solutions().begin()->solution.end(), ostream_iterator<varIndex>(cout, ""));
      cout << "\n";
    }


  } catch (ParseFailure& e) {

    cerr << "Crap: " << e.what() << "\n";
  }

  return 0;
}



