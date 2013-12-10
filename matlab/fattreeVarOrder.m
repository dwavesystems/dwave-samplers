function varOrder = fattreeVarOrder(chiDims, horiz, trunkI, p1, p2)
%fattreeVarOrder construct a Chimera fattree variable elimination order
%
% varOrder = fattreeVarOrder(chiDims, horiz, trunkI, p1, p2)
%
% chiDims: chimera dimensions [m n l]
% horiz: true for horizontal tree, false for vertical
% trunkI: row/column index of trunk
% p1: parity of left/top fat branch indices (true = even, false = odd)
% p2: same as p2 but for right/bottom fat branches

p1 = 1 + logical(p1);
p2 = 1 + logical(p2);

chiM = chiDims(1);
chiN = chiDims(2);
chiL = chiDims(3);
chiIndices = reshape(1:2*prod(chiDims), [ chiL 2 chiN chiM ]);

if horiz
  chiM = chiDims(2);
  chiN = chiDims(1);
  chiIndices = permute(flipdim(chiIndices, 2), [1 2 4 3]);
end

varOrder = [ ...
  ... vertical qubits in left fat branches
  reshape(chiIndices(:,1,1:trunkI-1,p1:2:chiM), 1, []) ... 
  ... vertical qubits in right fat branches
  reshape(chiIndices(:,1,chiN:-1:trunkI+1,p2:2:chiM), 1, []) ...
  ... horizontal qubits in all left branches
  reshape(chiIndices(:,2,1:trunkI-1,:), 1, []) ...
  ... horizontal qubits in all right branches
  reshape(chiIndices(:,2,chiN:-1:trunkI+1,:), 1, []) ...
  ... horizontal trunk qubits
  reshape(chiIndices(:,2,trunkI,:), 1, []) ...
  ... vertical trunk qubits
  reshape(chiIndices(:,1,trunkI,:), 1, []) ...
  ];

end
