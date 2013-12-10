function tables = quboTables(Q)
%quboTables converts a QUBO problem into table format
%
% tables = quboTables(Q)
%
% Both upper- and lower-triangular entries of Q are used.

Q = triu(Q) + triu(Q', 1);
[r c] = find(Q);

numTables = numel(r);
tables = struct('vars', cell(1,numTables), 'domSizes', [], 'values', []);

for ti=1:numTables
  if r(ti) == c(ti)
    tables(ti).vars = r(ti);
    tables(ti).domSizes = 2;
    tables(ti).values = full([0 Q(r(ti),r(ti))]);
  else
    tables(ti).vars = [ r(ti) c(ti) ];
    tables(ti).domSizes = [2 2];
    tables(ti).values = full([0 0 0 Q(r(ti),c(ti))]);
  end
end

end
