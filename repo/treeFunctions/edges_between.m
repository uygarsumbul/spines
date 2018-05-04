function [edgeLengths, edgeOrientations, edgeAreas, pillars] = edges_between(descendent,ancestor,tree)

edgeLengths = []; edgeOrientations = []; edgeAreas = []; pillars = [];
while descendent ~= ancestor
  if ~isempty(edgeLengths)
    pillars = [pillars; uint32(tree{descendent}{4}{5}(2:end) + numel(edgeLengths))];
    edgeOrientations = [edgeOrientations tree{descendent}{4}{3}(:,2:end)];
  else
    pillars = uint32(tree{descendent}{4}{5});
    edgeOrientations = [edgeOrientations tree{descendent}{4}{3}];
  end
  edgeLengths = [edgeLengths; tree{descendent}{4}{2}];
  edgeAreas = [edgeAreas; tree{descendent}{4}{4}];
  descendent = tree{descendent}{1};
end
