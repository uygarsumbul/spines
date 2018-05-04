function totalSubtreeLength = treeLength(tree,node,lengthOnly)
if nargin < 3
  lengthOnly = false;
  if nargin < 2
    node = rootFinder(tree);
  end
end
totalSubtreeLength = 0;
descendents = allDescendents(tree,node);
for kk = 1:numel(descendents)
  if lengthOnly
    totalSubtreeLength = totalSubtreeLength + sum(tree{descendents(kk)}{4}{2});
  else
    totalSubtreeLength = totalSubtreeLength + sum(tree{descendents(kk)}{4}{2}.*tree{descendents(kk)}{4}{4});
  end
end
