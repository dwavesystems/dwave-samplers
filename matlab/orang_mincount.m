%orang_mincount count minimum-energy solutions
%
% Usage:
%
%   [c e] = orang_mincount(tables, varOrder, maxComplexity, relEps, x0)
%
%
%
% Inputs:
%   tables: problem data formatted as a structure array with fields 'vars',
%     'domSizes', and 'values'.  The 'vars' field must contain positive
%     integral values listed in increasing order.  Variable domain sizes
%     must be positive integers and consistent across all tables.  The
%     dimensions of each field are irrelevant; only the number of elements
%     matters.
%
%   varOrder: variable elimination order.  Missing variables are clamped.
%
%   maxComplexity: maximum tree decomposition complexity.  The complexity
%     of a tree decomposition node is log_2 of the product of the node's
%     variables' domain sizes.  The complexity of the whole tree
%     decomposition is the maximum complexity of any node.  Note that if
%     all variables are binary (ie. domain size = 2), then the complexity
%     is the treewidth plus one.
%
%   relEps: relative error tolerance when comparing objective values.
%     Specifically, "a < b" is interpreted as a < b - relEps * abs(b).
%     Use this parameter with high-BOP input.  The value given should be
%     much smaller than the "step size" of the input values.  If missing or
%     empty, it defaults to zero.
%
%   x0: initial state.  Values are domain indices (ie. in 1:N for a domain
%     of size N).  Only the values of clamped variables are used.  This
%     parameter is optional.  If missing or empty, it defaults to all ones.
%
%
% Outputs:
%
%   c: number of states with energy e.
%
%   e: minimum energy.
