function descendents = allDescendents(tree, node)
descendents = tree{node}{2};
for child = 1:numel(tree{node}{2})
  descendents = [descendents allDescendents(tree, tree{node}{2}(child))];
end
