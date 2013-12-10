%orang_sample sample states from a Gibbs measure
%
% This function returns states sampled from a Gibbs measure.  That is, a
% state s is sampled with probability proportional to exp(E(s)), where E(s)
% is the objective value of s (the sum of the table values indexed by s).
%
% Note that this probabilty is normally written as exp(-b*E(s)) for some
% free parameter b.  This function requires the -b scaling to be done and
% the table level.  For example, given a QUBO matrix Q, construct the
% tables t with t = quboTables(-b*Q).
%
% Usage:
%
%   [pf s m] = orang_sample( ...
%       tables, varOrder, maxComplexity, samples, x0, minVars, rngSeed)
%
%   [pf m] = orang_sample( ...
%       tables, varOrder, maxComplexity, 0, x0, minVars, rngSeed)
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
%   samples: number of samples to return.  If missing or empty, default
%     value is 1.  If zero, no samples are returned, but the bucket tree
%     node information is still available.
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
%   rngSeed: random nunmber generator seed.  If missing or empty, the seed
%     is based on the current date and time.
%
%
% Outputs:
%
%   pf: log of the partition function.  The partition function Z is
%     exp(E(s)) summed over all states s.  Then the probability of a state
%     s is exp(E(s))/Z.
%
%   s: samples.  Each column is one samples.  Each solution entry is in the
%     range 1:N, where N is the size of the corresponding variable domain.
%
%   m: marginal probabilities of individual variables and pairs of variables.
%     This is a structure array with fields 'vars' and 'values'.  The 'vars'
%     field contains one or two variable indices and the 'values' field
%     contains a vector or matrix of marginal probabilities.  For each index
%     k, if m(k) contains individual variable marginals, then m(k).values(di)
%     is the probability that variable m(k).vars (a scalar in this case)
%     takes on domain index di.  If m(k) contains variable pair marginals,
%     then m(k).values(di, dj) is the probability that variables m(k).vars(1)
%     and m(k).vars(2) take domain indices di and dj, respectively.
