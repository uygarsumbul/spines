function root = rootFinder(tree);

root = 1;
while numel(tree{root}{1}) > 0
  root = tree{root}{1};
end
