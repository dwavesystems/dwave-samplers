function tables = isingTables(h, J)
%isingTables converts an Ising problem into table format
%
% tables = isingTables(h, J)
%
% Both upper- and lower-triangular entires of J are used.  The diagonal is
% ignored.


J = triu(J + J', 1);
[jr jc] = find(J);
numj = numel(jr);

hi = find(h);
numh = numel(hi);

numTables = numj + numh;
tables = struct('vars', cell(1,numTables), 'domSizes', [], 'values', []);

for ti=1:numh
  ii = hi(ti);
  tables(ti).vars = ii;
  tables(ti).domSizes = 2;
  tables(ti).values = full(h(ii) * [-1 1]);
end

for ii=1:numj
  ti = ii + numh;
  tables(ti).vars = [jr(ii) jc(ii)];
  tables(ti).domSizes = [2 2];
  tables(ti).values = full(J(jr(ii),jc(ii)) * [1 -1; -1 1]);
end  

end
