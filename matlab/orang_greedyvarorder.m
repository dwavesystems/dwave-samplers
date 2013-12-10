% orang_greedyvarorder generate a variable elimination order
%
% Usage:
%
%   varOrder = orang_greedyvarorder(
%       tables, maxComplexity, clampRanks, heuristic, selectionScale)
%
%
% Inputs:
%   tables: problem data formatted as a structure array with fields 'vars'
%     and 'domSizes'.  The 'vars' field must contain positive integral
%     values listed in increasing order.  Variable domain sizes must be
%     positive integers and consistent across all tables.  The dimensions
%     of each field are irrelevant; only the number of elements matters.
%     Note that other orang_* functions also require a 'values' field; this
%     field is not needed here but there is no harm in including it.
%
%   maxComplexity: maximum tree decomposition complexity of the returned
%     elimination order.  The complexity of a tree decomposition node is
%     log_2 of the product of the node's variables' domain sizes.  The
%     complexity of the whole tree decomposition is the maximum complexity
%     of any node.  Note that if all variables are binary (ie. domain size
%     = 2), then the complexity is the treewidth plus one.
%
%   clampRanks: variable clamping hint.  Clamped variables are those not
%     included in the elimination order.  This parameter must be a vector
%     indexed by variable.  Values indicate relative preferences for
%     clamping variables.  A negative value causes the corresponding
%     variable to be clamped immediately.  If, during the course of the
%     algorithm, a variable needs to be clamped, the choice is confined to
%     those variables with minimum clamp rank.  After clamping a variable,
%     the clamp ranks of all variables with higher clamp rank a decreased
%     by one.
%
%   heuristic: variable elimination heuristic to use.  Choices are:
%     'mindeg'   - minimum degree, 'wmindeg'  - weighted (by domain size)
%     degree, 'minfill'  - minimum number of new edges introduced,
%     'wminfill' - weighted min fill.
%
%   selectionScale: scaling factor for range size when randomly
%     eliminating or clamping variables.  Normally, the variable to be
%     eliminated is chosen randomly from all variables with lowest cost
%     whose elimination don't violate the maximum complexity bound. The
%     size of that set of variables is multiplied by this value, possibly
%     shrinking it (<1) or enlarging it (>1).  The set is enlarged by
%     adding variables with low elimination cost.  The situation is similar
%     when clamping a variable, though the variable clamp ranks are always
%     respected.
%     
%
% Output:
%
%   varOrder: the variable eliminiation order.  Missing variables are
%       clamped.
