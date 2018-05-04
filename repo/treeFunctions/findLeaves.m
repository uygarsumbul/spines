function leaves = findLeaves(tree,nodeSubset)
if nargin < 2
  nodeSubset = 1:numel(tree);
end

leaves = [];
for node = 1:numel(nodeSubset)
  if isempty(tree{nodeSubset(node)}{2})
    leaves = [leaves nodeSubset(node)];
  end
end
