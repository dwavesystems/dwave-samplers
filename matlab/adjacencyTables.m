function tables = adjacencyTables(A, domSizes)
% adjacencyTables converts an adjacency matrix into table format
%
% tables = adjacencyTables(A, domSizes)
%
% The return value is a structure array with 'vars' and 'domSizes' fields.
% Note that there is no 'values' field; the purpose of this function is to
% ease the creation of a 'tables' parameter of the orang_greedyvarorder
% function, which does not use the 'values' field required by other orang_*
% functions.
%
% Both upper- and lower-triangular regions of A are used.  The domSizes
% parameter is optional and defaults to all 2's.  

[r c] = find(tril(A|A', -1));
numTables = numel(r);
tables = struct('vars', cell(1,numTables), 'domSizes', []);

if nargin < 2
  domSizes = 2 * ones(1, numTables);
end

for ii=1:numTables
  tables(ii).vars = [c(ii) r(ii)];
  tables(ii).domSizes = domSizes(tables(ii).vars);
end

end
