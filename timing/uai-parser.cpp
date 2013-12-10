#include <cstddef>
#include <algorithm>
#include <ios>
#include <istream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <functional>
#include <utility>

#include <boost/foreach.hpp>
#include <boost/io/ios_state.hpp>

#include <orang/orang.h>

#include "uai.h"

using std::size_t;
using std::istream;
using std::ios_base;
using std::ios;
using std::reverse_copy;
using std::sort;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::greater;
using std::make_pair;

using boost::io::ios_exception_saver;

using orang::varIndex;
using orang::VarVector;
using orang::SizeVector;
using orang::TableVar;
using orang::Table;



class ReindexedTableIter {
private:
  SizeVector cur_;
  VarVector scopeIndices;
  const vector<TableVar>& vars_;
  Table<double>::iterator iter_;
public:
  ReindexedTableIter(Table<double>& table, const VarVector& inputScope) :
    cur_(inputScope.size(), 0),
    scopeIndices(),
    vars_(table.vars()),
    iter_(table.begin()) {

    typedef map<varIndex, varIndex, greater<varIndex> > simap_type;
    simap_type scopeIndexMap;
    varIndex index = 0;
    BOOST_FOREACH( varIndex v, inputScope ) {
      scopeIndexMap[v] = index++;
    }

    scopeIndices.reserve(inputScope.size());
    BOOST_FOREACH( const simap_type::value_type& entry, scopeIndexMap ) {
      scopeIndices.push_back(entry.second);
    }
  }

  ReindexedTableIter& operator++() {
    for (size_t i = 0, n = scopeIndices.size(); i < n; ++i) {
      const TableVar& var = vars_[scopeIndices[i]];
      if (cur_[i] + 1 < var.domSize) {
        iter_ += var.stepSize;
        ++cur_[i];
        break;
      } else {
        iter_ -= cur_[i] * var.stepSize;
        cur_[i] = 0;
      }
    }

    return *this;
  }

  double& operator*() {
    return *iter_;
  }
};


ParsedProblem parseUaiProblem(istream& in) {
  ios_exception_saver ies(in);
  in.exceptions(istream::badbit | istream::failbit);

  // validate problem type
  static const string MARKOV_PROBLEM_TYPE = "MARKOV";
  string problemType;
  in >> problemType;
  if (problemType != MARKOV_PROBLEM_TYPE) {
    throw ParseFailure("Unknown problem type: " + problemType);
  }

  // read number of variables and domain sizes
  ParsedProblem result;
  varIndex numVars;
  in >> numVars;
  result.domSizes.reserve(numVars);
  for (varIndex i = 0; i < numVars; ++i) {
    in >> result.domSizes[i];
  }

  // read tables scopes
  size_t numTables;
  in >> numTables;
  vector<VarVector> inputScopes(numTables);
  for (size_t i = 0; i < numTables; ++i) {
    varIndex scopeSize;
    in >> scopeSize;
    VarVector& inputScope = inputScopes[i];
    inputScope.resize(scopeSize);
    for (varIndex j = 0; j < scopeSize; ++j) {
      in >> inputScope[j];
      if (inputScope[j] >= numVars) {
        ostringstream msg;
        msg << "Invalid variable (" << inputScope[j] << ") in table " << i
            << ".  Maximum variable index is " << (numVars - 1) << ".";
        throw ParseFailure(msg.str());
      }
    }
  }

  // build tables
  result.tables.reserve(numTables);
  for (size_t i = 0; i < numTables; ++i) {
    VarVector tableScope = inputScopes[i];
    sort(tableScope.begin(), tableScope.end());
    SizeVector tableDomSizes;
    tableDomSizes.reserve(tableScope.size());
    for (size_t j = 0, n = tableScope.size(); j < n; ++j) {
      tableDomSizes.push_back(result.domSizes[tableScope[j]]);
    }
    Table<double>::smartptr tablePtr( new Table<double>(tableScope, tableDomSizes) );

    size_t numEntries;
    in >> numEntries;
    if (numEntries != tablePtr->size()) {
      ostringstream msg;
      msg << "Given number of entries (" << numEntries << ") for table " << i
          << " is wrong.  It should be " << tablePtr->size() << ".";
      throw ParseFailure(msg.str());
    }

    ReindexedTableIter entryIter(*tablePtr, inputScopes[i]);
    for (size_t k = 0; k < numEntries; ++k) {
      in >> *entryIter;
      ++entryIter;
    }

    result.tables.push_back(tablePtr);
  }

  return result;
}


ParsedEvidence parseUaiEvidence(istream& in, const ParsedProblem& pp) {
  ios_exception_saver ies(in);
  in.exceptions(istream::badbit | istream::failbit);

  in >> std::ws;
  if (in.eof()) {
    return ParsedEvidence(0);
  }

  size_t numEvidenceLines;
  in >> numEvidenceLines;

  ParsedEvidence evidence(numEvidenceLines);

  for (size_t i = 0; i < numEvidenceLines; ++i) {
    varIndex numClampedVars;
    in >> numClampedVars;
    for (size_t j = 0; j < numClampedVars; ++j) {
      varIndex var;
      varIndex val;
      in >> var >> val;
      if (var >= pp.domSizes.size()) { // || val >= pp.domSizes[var]) {
        ostringstream msg;
        msg << "Evidence line " << i << " lists variable " << var
            << " but the maximum valid variable index for this problem is " << (pp.domSizes.size() - 1) << ".";
        throw ParseFailure(msg.str());
      }
      if (val >= pp.domSizes[var]) {
        ostringstream msg;
        msg << "Observed value for variable " << var << " on evidence line " << i << " is " << val
            << " but its maximum valid domain index is " << pp.domSizes[var] << ".";
        throw ParseFailure(msg.str());
      }
      evidence[i].push_back(make_pair(var, val));
    }
  }

  return evidence;
}

