%orang_minsum solve a min-sum problem
%
% Usage:
%
%   [e x t] = orang_minsum( ...
%       tables, varOrder, maxComplexity, maxSolutions, x0, minVars)
%
%   [e t] = orang_minsum(tables, varOrder, maxComplexity, 0, x0, minVars)
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
%   maxSolutions: maximum number of solutions and objective values to
%     return. If missing or empty, default value is 1.  If zero, no
%     solutions are returned, but the bucket tree node information is still
%     available.
%
%   x0: initial state.  Values are domain indices (ie. in 1:N for a domain
%     of size N).  Only the values of clamped variables are used.  This
%     parameter is optional.  If missing or empty, it defaults to all ones.
%
%   minVars: by default, the number of variables is the largest variable
%     index appearing in an input table.  This parameter can be used to
%     increase the number of variables (affecting things such as output
%     sample sizes, required x0 size, and variables allowed to appear in
%     varOrder).
%                  
%
% Outputs:
%
%   e: row vector of up to maxSolutions best objective values.  If
%     maxSolutions is zero, the optimal objective value is still returned.
%
%   x: matrix of best solutions.  Each column is one solution and
%     size(x,2) == numel(e).  Solutions with equal objective values are
%     sorted lexicographically.  Each solution entry is in the range 1:N,
%     where N is the size of the corresponding variable domain.
%
%   t: bucket tree node information.  This is a structure array with fields
%     'nodeVar', 'sepVars', and 'tables'.  The 'nodeVar' field is the
%     variable eliminated at the corresponding bucket tree node.  Each
%     non-clamped variable appears as a nodeVar exactly once.  The
%     'sepVars' field is the set of variables appearing in both the bucket
%     tree node and its parent. The 'tables' field contains all the tables
%     (input, lambda, and pi) for the node.  The format of the tables is
%     the same as the input parameter.
