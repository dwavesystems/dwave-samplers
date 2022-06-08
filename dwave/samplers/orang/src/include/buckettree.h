/**
# Copyright 2019 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.
#
# =============================================================================
*/
#ifndef INCLUDED_ORANG_BUCKETTREE_H
#define INCLUDED_ORANG_BUCKETTREE_H

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <vector>
#include <memory>

#include <base.h>
#include <exception.h>
#include <table.h>
#include <treedecomp.h>
#include <task.h>
#include <merger.h>

namespace orang {

template<typename Y>
struct NodeTables {
  typedef Table<Y> table_type;
  typedef typename table_type::const_smartptr const_table_smartptr;
  Var nodeVar;
  VarVector sepVars;
  std::vector<const_table_smartptr> tables;
};


template<typename T>
class BucketTree {
public:
  typedef T task_type;
  typedef typename task_type::value_type value_type;
  typedef typename task_type::solution_type solution_type;
  typedef NodeTables<value_type> nodetables_type;

private:
  typedef typename task_type::table_type table_type;
  typedef typename table_type::const_smartptr const_table_smartptr;
  typedef typename task_type::table_vector table_vector;
  typedef typename task_type::marginalizer_type marginalizer_type;
  typedef typename task_type::marginalizer_smartptr marginalizer_smartptr;
  typedef typename task_type::solvablemarginalizer_type solvablemarginalizer_type;

  struct Node {
    typedef std::shared_ptr<Node> smartptr;

    std::vector<smartptr> children;
    table_vector baseTables;
    table_vector lambdaTables;
    const_table_smartptr piTable;
    marginalizer_smartptr marginalizer;
  };

  typedef typename Node::smartptr node_smartptr;

  const bool solvable_;
  const bool hasNodeTables_;
  const task_type& task_;
  const DomIndexVector x0_;
  const TableMerger<T> mergeTables_;
  value_type problemValue_;
  std::vector<node_smartptr> roots_;
  std::size_t numNodes_;
  std::vector<NodeTables<value_type> > nodeTables_;

  typename Node::smartptr buildNode(
      const TreeDecompNode& dNode,
      table_vector* parentTables,
      std::vector<value_type>& rootValues,
      const DomIndexVector& x0);

  void buildNodeTables(node_smartptr& node, const TreeDecompNode& dNode);

  void solveRecursive(const node_smartptr& n, solution_type& s) const {
    dynamic_cast<solvablemarginalizer_type&>(*n->marginalizer).solve(s);
    for (const auto &c: n->children) {
      solveRecursive(c, s);
    }
  }

public:
  BucketTree(
      const task_type& task,
      const TreeDecomp& decomp,
      const DomIndexVector& x0,
      bool solvable,
      bool hasNodeTables) :
        solvable_(solvable),
        hasNodeTables_(hasNodeTables),
        task_(task),
        x0_(x0),
        mergeTables_(TableMerger<task_type>(task)),
        problemValue_(),
        roots_(),
        numNodes_(0),
        nodeTables_() {

    if (x0_.size() != task_.numVars()) {
      throw InvalidArgumentException("x0 has incorrect size");
    }

    std::vector<value_type> rootValues;
    rootValues.reserve(decomp.roots().size());
    roots_.reserve(decomp.roots().size());
    for (const auto &dNode: decomp.roots()) {
      roots_.push_back(buildNode(*dNode, 0, rootValues, x0));
      if (!hasNodeTables_) {
        roots_.back()->lambdaTables.clear();
        roots_.back()->baseTables.clear();
      }
    }

    problemValue_ = task_.problemValue(rootValues, x0, decomp.clampedVars());

    if (hasNodeTables_) {
      nodeTables_.reserve(numNodes_);
      TreeDecompNode::node_vector::const_iterator dRootsIter = decomp.roots().begin();
      for (auto &r: roots_) {
        buildNodeTables(r, **dRootsIter);
        ++dRootsIter;
      }
    }

    if (!solvable_) {
      roots_.clear();
    }
  }

  const task_type& task() const { return task_; }

  bool solvable() const { return solvable_; }

  bool hasNodeTables() const { return hasNodeTables_; }

  value_type problemValue() const { return problemValue_; }

  solution_type solve() const {
    if (solvable_) {
      solution_type s = task_.initSolution(x0_);
      for (const auto &r: roots_) {
        solveRecursive(r, s);
      }
      return s;
    } else {
      throw OperationUnavailable();
    }
  }

  const std::vector<nodetables_type>& nodeTables() const {
    if (hasNodeTables_) {
      return nodeTables_;
    } else {
      throw OperationUnavailable();
    }
  }
};


template<typename T>
typename BucketTree<T>::node_smartptr BucketTree<T>::buildNode(
    const TreeDecompNode& dNode,
    typename BucketTree<T>::table_vector* parentTables,
    std::vector<value_type>& rootValues,
    const DomIndexVector& x0) {

  ++numNodes_;

  node_smartptr n = node_smartptr( new Node );
  n->baseTables = task_.baseTables(dNode, x0);
  n->children.reserve(dNode.children().size());
  for (const auto &cdn: dNode.children()) {
    n->children.push_back( buildNode(*cdn, &n->lambdaTables, rootValues, x0) );
  }

  if (solvable_) {
    DomIndexVector inDomSizes;
    inDomSizes.reserve(dNode.sepVars().size());
    for (auto v: dNode.sepVars()) {
      inDomSizes.push_back(task_.domSize(v));
    }

    n->marginalizer = task_.solvableMarginalizer(dNode.sepVars(), inDomSizes,
        dNode.nodeVar(), task_.domSize(dNode.nodeVar()));
  } else {
    n->marginalizer = task_.marginalizer();
  }

  table_vector inTables;
  inTables.reserve(n->baseTables.size() + n->lambdaTables.size());
  copy(n->baseTables.begin(), n->baseTables.end(), back_inserter(inTables));
  copy(n->lambdaTables.begin(), n->lambdaTables.end(), back_inserter(inTables));
  const_table_smartptr pLambdaTable = mergeTables_(
      dNode.sepVars(), inTables.begin(), inTables.end(), *n->marginalizer);

  if (parentTables) {
    parentTables->push_back(pLambdaTable);
  } else {
    rootValues.push_back((*pLambdaTable)[0]);
  }

  if (!hasNodeTables_) {
    n->baseTables.clear();
    n->lambdaTables.clear();
  }

  return n;
}

template<typename T>
void BucketTree<T>::buildNodeTables(typename BucketTree<T>::Node::smartptr& node, const TreeDecompNode& dNode) {

  using std::size_t;
  using std::copy;
  using std::back_inserter;

  nodeTables_.push_back(NodeTables<value_type>());
  NodeTables<value_type>& nt = nodeTables_.back();

  size_t numTables = node->baseTables.size() + node->lambdaTables.size();

  nt.nodeVar = dNode.nodeVar();
  nt.sepVars = dNode.sepVars();

  nt.tables.reserve(numTables + 1);
  copy(node->baseTables.begin(), node->baseTables.end(), back_inserter(nt.tables));
  copy(node->lambdaTables.begin(), node->lambdaTables.end(), back_inserter(nt.tables));
  if (node->piTable) {
    nt.tables.push_back(node->piTable);
  }

  size_t numChildren = dNode.children().size();

  if (numChildren > 0) {
    table_vector inTables;
    inTables.reserve(numTables);
    copy(node->lambdaTables.begin() + 1, node->lambdaTables.end(), back_inserter(inTables));
    copy(node->baseTables.begin(), node->baseTables.end(), back_inserter(inTables));
    if (node->piTable) {
      inTables.push_back(node->piTable);
    }

    const TreeDecompNode::node_vector& dChildren = dNode.children();

    if (!inTables.empty()) {
      TableMerger<task_type> merge(task_);
      marginalizer_smartptr mrg = task_.marginalizer();


      for (size_t i = 0; i < numChildren; ++i) {
        if (i > 0) {
          inTables[i - 1] = node->lambdaTables[i - 1];
        }

        node->children[i]->piTable = merge(
            dChildren[i]->sepVars(), inTables.begin(), inTables.end(), *mrg);
      }
    }

    for (size_t i = 0; i < numChildren; ++i) {
      buildNodeTables(node->children[i], *dChildren[i]);
    }
  }

  node->baseTables.clear();
  node->lambdaTables.clear();
  node->piTable.reset();
}

} // namespace orang

#endif
