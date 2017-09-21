function ancestors = allAncestors(tree, node)
if isempty(tree{node}{1})
  ancestors = [];
  return;
end
ancestors = tree{node}{1};
if ~isempty(tree{ancestors}{1})
  ancestors = [ancestors allAncestors(tree, ancestors)];
end
