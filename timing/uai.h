#ifndef INCLUDED_UAI_H
#define INCLUDED_UAI_H

#include <istream>
#include <string>
#include <vector>
#include <map>
#include <utility>

#include <orang/orang.h>

class ParseFailure {
private:
  std::string what_;
public:
  ParseFailure(const std::string& what = "") : what_(what) {}
  const std::string what() const { return what_; }
};

struct ParsedProblem {
  orang::SizeVector domSizes;
  std::vector<orang::Table<double>::smartptr > tables;
};

typedef std::vector<std::pair<orang::varIndex, orang::varIndex> > Evidence;
typedef std::vector<Evidence> ParsedEvidence;

ParsedProblem parseUaiProblem(std::istream& in);

ParsedEvidence parseUaiEvidence(std::istream& in, const ParsedProblem& pp);

void limitMemory();

#endif
