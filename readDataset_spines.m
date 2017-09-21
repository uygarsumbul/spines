function allTrees = readDataset_spines(code,directory,options)

if ~isfield(options,'normalize') || isempty(options.normalize); normalize = false; else; normalize = options.normalize; end;
if ~isfield(options,'absoluteLengths') || isempty(options.absoluteLengths); absoluteLengths = true; else; absoluteLengths = options.absoluteLengths; end;
if ~isfield(options,'neuronParts') || isempty(options.neuronParts); neuronParts = [-10:10]; else; neuronParts = options.neuronParts; end;
if ~isfield(options,'anisotropyDivisors') || isempty(options.anisotropyDivisors); anisotropyDivisors = [1 1 1]; else; anisotropyDivisors = options.anisotropyDivisors; end;
if ~isfield(options,'SCALECONSTANT') || isempty(options.SCALECONSTANT); SCALECONSTANT = 1; else; SCALECONSTANT = options.SCALECONSTANT; end;
if ~isfield(options,'removeShortLeaves') || isempty(options.removeShortLeaves); removeShortLeaves = false; else; removeShortLeaves = options.removeShortLeaves; end;
if ~isfield(options,'pruneRatio') || isempty(options.pruneRatio); pruneRatio = 1; else; pruneRatio = options.pruneRatio; end;
if ~isfield(options,'removeBranches') || isempty(options.removeBranches); removeBranches = false; else; removeBranches = options.removeBranches; end;
if ~isfield(options,'resolution') || isempty(options.resolution); resolution = [1 1 1]; else; resolution = options.resolution; end;

allTrees = cell(0); allNodesAndEdges = cell(0);

for kk = 1:numel(code)
  thisFileName = strcat(directory,code{kk}{1}); %,code{kk}{2});
  [nodes,edges,radii,~,features,~] = readSWCfile(thisFileName,neuronParts,resolution);
  nodes(:,1) = nodes(:,1)/anisotropyDivisors(1); nodes(:,2) = nodes(:,2)/anisotropyDivisors(2); nodes(:,3) = nodes(:,3)/anisotropyDivisors(3); 
  [tree,rawLength,nodes,edges]=generateIrreducibleDoubleLinkedTree(nodes,edges,radii,features,false); % 1:soma, 2:axon, 3:basal dendrite, 4:apical dendrite
  if removeShortLeaves; tree = pruneShortLeaves(tree,treeLength(tree)*pruneRatio); end;
  if normalize; tree = normalizeNeuron(tree); end;
  allTrees{end+1} = tree;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges,radii,nodeTypes,features,abort] = readSWCfile(fileName,validNodeTypes,resolution)
abort = false; nodes = []; edges = []; nodeTypes = [];
validNodeTypes = setdiff(validNodeTypes,1);
[nodeID, nodeType, xPos, yPos, zPos, radii, parentNodeID,d1,d2,f1,f2,f3,f4,f5,f6,f7,f8,f9] = textread(fileName, '%u%d%f%f%f%f%d%d%d%f%f%f%f%f%f%f%f%f','commentstyle', 'shell');
xPos = xPos * resolution(1); yPos = yPos * resolution(2); zPos = zPos * resolution(3);
features = [d1,d2,f1,f2,f3,f4,f5,f6,f7,f8,f9];
if ~any(parentNodeID==-1)
  disp(strcat('root not found in ',fileName));
  nodes = []; edges = []; radii = []; nodeTypes = []; features = []; abort = true;
  return;
end
nodeType(find(parentNodeID==-1))=1; % Every tree should start from a node of type 1 (soma)
firstSomaNode = find(nodeType == 1 & parentNodeID == -1, 1);
somaNodes = find(nodeType == 1); somaX = mean(xPos(somaNodes)); somaY = mean(yPos(somaNodes)); somaZ = mean(zPos(somaNodes));
% set soma radius to 0, to avoid multiple counting of soma volume through different subsets. ignore soma volume
somaRadius = 0; %mean(radii(somaNodes));
xPos(firstSomaNode) = somaX; yPos(firstSomaNode) = somaY; zPos(firstSomaNode) = somaZ; radii(firstSomaNode) = somaRadius; features(firstSomaNode,:) = zeros(1,size(features,2));
parentNodeID(ismember(parentNodeID,somaNodes)) = firstSomaNode; % assign a single soma parent
nodesToDelete = setdiff(somaNodes,firstSomaNode); % delete all the soma nodes except for the firstSomaNode
nodeID(nodesToDelete)=[]; nodeType(nodesToDelete)=[]; xPos(nodesToDelete)=[]; yPos(nodesToDelete)=[]; zPos(nodesToDelete)=[]; radii(nodesToDelete)=[]; parentNodeID(nodesToDelete)=[]; features(nodesToDelete,:)=[];
for kk = 1:numel(nodeID)
  while ~any(nodeID==kk)
    nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
    parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
  end
end
validNodes = nodeID(ismember(nodeType,validNodeTypes));
additionalValidNodes = [];
for kk = 1:numel(validNodes)
  thisParentNodeID = parentNodeID(validNodes(kk)); thisParentNodeType = nodeType(thisParentNodeID);
  while ~ismember(thisParentNodeType,validNodeTypes)
    if thisParentNodeType == 1
      break;
    end
    additionalValidNodes = union(additionalValidNodes, thisParentNodeID); nodeType(thisParentNodeID) = validNodeTypes(1);
    thisParentNodeID = parentNodeID(thisParentNodeID); thisParentNodeType = nodeType(thisParentNodeID);
  end
end
validNodes = [firstSomaNode; validNodes; additionalValidNodes']; validNodes = unique(validNodes);
nodeID = nodeID(validNodes); nodeType = nodeType(validNodes); parentNodeID = parentNodeID(validNodes);
xPos = xPos(validNodes); yPos = yPos(validNodes); zPos = zPos(validNodes); radii = radii(validNodes); features = features(validNodes, :);
for kk = 1:numel(nodeID)
  while ~any(nodeID==kk)
    nodeID(nodeID>kk) = nodeID(nodeID>kk)-1;
    parentNodeID(parentNodeID>kk) = parentNodeID(parentNodeID>kk)-1;
  end
end
nodes = [xPos yPos zPos];
edges = [nodeID parentNodeID];
edges(any(edges==-1,2),:) = [];
nodeTypes = unique(nodeType)';
%[nodes,edges,radii,features] = removeImproperRootNodes(nodes,edges,radii,features); % remove the root nodes until the first node
%[nodes,edges,radii,features] = removeZeroLengthEdges(nodes,edges,radii,features); % remove zero branches
nodes = nodes + 1; % added on Mar 7, 2012 - band tracing starts at 1 - relevant if the tree will be warped

somaNode = find(nodeType == 1, 1);
% swap on edges
edges(edges==1) = 0; edges(edges==somaNode) = 1; edges(edges==0) = somaNode;
% swap on nodes
tmp = nodes(1,:); nodes(1,:) = nodes(somaNode,:); nodes(somaNode,:) = tmp;
% swap on radii
tmpR = radii(1); radii(1) = radii(somaNode); radii(somaNode) = tmpR;
% swap on features
tmpF = features(1,:); features(1,:) = features(somaNode,:); features(somaNode,:) = tmpF;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges,radii,features] = removeImproperRootNodes(nodes,edges,radii,features)
if size(nodes,1) < 2
  return;
end
% node 1 is assigned as the root node
root = 1;
expandedEdges = [edges edges(:,1)];
child = expandedEdges([zeros(size(edges,1),1) edges==root]>0);
while numel(child) < 2
  nodes = [nodes(1:root-1,:); nodes(root+1:end,:)];
  radii = [radii(1:root-1); radii(root+1:end)];
  features = [features(1:root-1,:); features(root+1:end,:)];
  edges(any(edges==root,2),:) = [];
  edges(edges>root) = edges(edges>root) - 1;
  expandedEdges = [edges edges(:,1)];
  if child > root
    child = child-1;
  end
  root = child;
  child = expandedEdges([zeros(size(edges,1),1) edges==root]>0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges,radii,features] = removeZeroLengthEdges(nodes,edges,radii,features)
kk = 1;
while kk < size(edges,1)
  node1 = edges(kk,1); node2 = edges(kk,2);
  if norm(nodes(node1,:)-nodes(node2,:)) < 1e-10 % BE CAREFUL / ARBITRARY
    edges(edges==node1) = node2;
    edges(edges>node1) = edges(edges>node1)-1;
    edges(kk,:) = []; % remove edge
    nodes(node1,:) = []; % remove node
    radii(node1) = []; % remove area
    features(node1,:) = []; % remove features
  else
    kk = kk + 1;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges,radii] = removeOutOfBoundsBranches(nodes,edges,radii,bounds)
if size(nodes,1) < 2
  return;
end
thisNode=find(nodes(:,1)<bounds.minX | nodes(:,1)>bounds.maxX | nodes(:,2)<bounds.minY | nodes(:,2)>bounds.maxY | nodes(:,3)<bounds.minZ | nodes(:,3)>bounds.maxZ,'first',1);
while ~isempty(thisNode)
  selfOrDescendant = selfOrDescendantFromEdges(edges, thisNode, thisNode);
  edges(ismember(edges(:,1),selfOrDescendant) | ismember(edges(:,2),selfOrDescendant),:) = [];
  nodes(selfOrDescendant,:) = [];
  radii(selfOrDescendant) = [];
  selfOrDescendent = sort(selfOrDescendant,'descend');
  for kk=1:numel(selfOrDescendant)
    nn=selfOrDescendant(kk);
    edges(edges>nn) = edges(edges>nn)-1;
  end
  thisNode=find(nodes(:,1)<bounds.minX | nodes(:,1)>bounds.maxX | nodes(:,2)<bounds.minY | nodes(:,2)>bounds.maxY | nodes(:,3)<bounds.minZ | nodes(:,3)>bounds.maxZ,'first',1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sod = selfOrDescendantFromEdges(edges, sod, node)
children = edges(find(edges(:,1)==node),2);
sod = [sod children];
if ~isempty(children)
  for kk=1:numel(children)
    sod = selfOrDescendantFromEdges(edges, sod, children(kk));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tree,rawLength,nodes,edges] = generateIrreducibleDoubleLinkedTree(nodes,edges,radii,features,uniformCSAreas)
areas = pi*radii.^2;
% 1 is assumed to be the root node. edges is nx2: each row is [child parent]
h = hist(edges(:,2),[1:size(nodes,1)]); irreducibleNodes = union(find(h~=1),1); % 1 is the root node
% put all the irreducible nodes at the beginning
for kk = 1:numel(irreducibleNodes)
  % swap on edges
  edges(edges==kk) = 0; edges(edges==irreducibleNodes(kk)) = kk; edges(edges==0) = irreducibleNodes(kk);
  % swap on nodes
  tmp = nodes(kk,:); nodes(kk,:) = nodes(irreducibleNodes(kk),:); nodes(irreducibleNodes(kk),:) = tmp;
  % swap on areas
  tmpA = areas(kk); areas(kk) = areas(irreducibleNodes(kk)); areas(irreducibleNodes(kk)) = tmpA;
  % swap on features
  tmpF = features(kk,:); features(kk,:) = features(irreducibleNodes(kk),:); features(irreducibleNodes(kk),:) = tmpF;
end
numelNodes = numel(irreducibleNodes);
%initialize tree with root as 1
tree{1}{1} = []; tree{1}{3} = nodes(1,:); tree{1}{4}{1} = nodes(1,:); tree{1}{4}{2} = 0; tree{1}{4}{3} = [[0;0;0] [0;0;0]]; tree{1}{4}{4} = 0;
for kk = 1:numelNodes
  tree{kk}{2} = [];
end
rawLength = 0;
for kk = 2:numelNodes
  tmpParent = edges(find(edges(:,1)==kk),2); % assume that the edges are ordered pairs: (child, parent)
  path = nodes(kk,:); thisFeature = features(kk, :);
  if uniformCSAreas
    pathAreas = 1;
  else
    pathAreas = (areas(kk)+areas(tmpParent)+sqrt(areas(kk)*areas(tmpParent)))/3;
  end
  while tmpParent > numelNodes
    newTmpParent = edges(find(edges(:,1)==tmpParent),2);
    path = [path; nodes(tmpParent,:)]; thisFeature = [thisFeature; features(tmpParent,:)];
    if uniformCSAreas
      pathAreas = [pathAreas; 1];
    else
      pathAreas = [pathAreas; (areas(tmpParent)+areas(newTmpParent)+sqrt(areas(tmpParent)*areas(newTmpParent)))/3]; % now modeled as a cylinder
    end
    tmpParent = newTmpParent;
  end
  path = [path; nodes(tmpParent,:)]; thisFeature = [thisFeature; features(tmpParent,:)]; rawPathLengths = sqrt(sum(diff(path,1,1).^2,2)); rawLength = rawLength + sum(rawPathLengths);
  tree{kk}{1} = tmpParent; tree{tmpParent}{2} = [tree{tmpParent}{2} kk];
  tree{kk}{3} = nodes(kk,:);
  tree{kk}{4}{1} = path; tree{kk}{4}{2} = rawPathLengths; tree{kk}{4}{4} = pathAreas; tree{kk}{4}{5} = thisFeature;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tree = pruneShortLeaves(tree,lengthThreshold)
if numel(tree)<2
  return;
end
root = rootFinder(tree);
leaves = findLeaves(tree); leaves = setdiff(leaves, tree{root}{2})';
leaves = [leaves zeros(size(leaves))];
for kk = 1:size(leaves,1)
  leaves(kk,2) = sum(tree{leaves(kk)}{4}{2});
end
leaves = sortrows(leaves,2); leaves = leaves(:,1);

counter = 1;
while counter <= numel(leaves)
  if (sum(tree{leaves(counter)}{4}{2}) > lengthThreshold)
    leaves(counter) = [];
  else
    counter = counter + 1;
  end
end
for counter = 1:numel(leaves)
  leaf = leaves(counter);
  if (sum(tree{leaf}{4}{2}) <= lengthThreshold)
    parent = tree{leaf}{1}; siblings = setdiff(tree{parent}{2},leaf);
    % remove the leaf from the child relationship
    tree{parent}{2} = setdiff(tree{parent}{2}, leaf);
    % change node numbers of parents/children if larger than this leaf's number
    for node = 1:numel(tree)
      if tree{node}{1} > leaf
        tree{node}{1} = tree{node}{1}-1;
      end
      children = tree{node}{2};
      tree{node}{2}(children > leaf) = children(children > leaf)-1;
    end
    % update the parent, sibling, leaf numbers
    parent(parent>leaf) = parent(parent>leaf)-1;
    siblings(siblings>leaf) = siblings(siblings>leaf)-1;
    leaves(leaves>leaf) = leaves(leaves>leaf)-1;
    % remove leaf from tree
    tree(leaf) = [];
    % update the root
    root = rootFinder(tree);
    % if the removal made the parent node a reducible node...
    if (numel(siblings) == 1) && (parent ~= root)
      % the grandparent removes the reducible parent from the list of children, and inherits the single child
      tree{tree{parent}{1}}{2} = [setdiff(tree{tree{parent}{1}}{2},parent) siblings];
      % the single child's parent becomes the grandparent
      tree{siblings}{1} = tree{parent}{1};
      % the single child adds the reducible node's geometric information to itself
      tree{siblings}{4}{1} = [tree{siblings}{4}{1}(1:end-1,:); tree{parent}{4}{1}]; % critical nodes
      tree{siblings}{4}{2} = [tree{siblings}{4}{2}; tree{parent}{4}{2}]; % lengths
%      tree{siblings}{4}{3} = [tree{siblings}{4}{3} tree{parent}{4}{3}(:,2:end)]; % orientations
      tree{siblings}{4}{4} = [tree{siblings}{4}{4}; tree{parent}{4}{4}]; % areas
      tree{siblings}{4}{5} = [tree{siblings}{4}{5}; tree{parent}{4}{5}]; % features
%      if isempty(tree{siblings}{4}{5})
%        tree{siblings}{4}{5} = tree{parent}{4}{5};
%      else
%        tree{siblings}{4}{5} = [tree{siblings}{4}{5}; uint32(tree{parent}{4}{5}(2:end)+numel(tree{siblings}{4}{5})-1)];
%      end
      % change node numbers of parents/children if larger than the parent's number
      for node = 1:numel(tree)
        if tree{node}{1} > parent
          tree{node}{1} = tree{node}{1}-1;
        end
        children = tree{node}{2};
        tree{node}{2}(children > parent) = children(children > parent)-1;
      end
      % update leaf numbers
      leaves(leaves>parent) = leaves(leaves>parent)-1;
      % remove the reducible parent from tree
      tree(parent) = [];
    end
  end
  tmpLeaves = leaves(counter+1:end); tmpLeaves = [tmpLeaves zeros(size(tmpLeaves))];
  if ~isempty(tmpLeaves)
    for kk = 1:size(tmpLeaves,1)
      tmpLeaves(kk,2) = sum(tree{tmpLeaves(kk,1)}{4}{2});
    end
    tmpLeaves = sortrows(tmpLeaves,2);
    leaves(counter+1:end) = tmpLeaves(:,1);
 end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [localMass,newNodes] = switchToMassRepresentation(nodes,edges)
localMass = zeros(size(nodes,1),1); newNodes = zeros(size(nodes,1),3);
for kk=1:size(nodes,1);
  parent = edges(find(edges(:,1)==kk),2);
  if ~isempty(parent)
    localMass(kk) = norm(nodes(parent,:)-nodes(kk,:)); newNodes(kk,:) = (nodes(parent,:)+nodes(kk,:))/2;
  else
    localMass(kk) = 0; newNodes(kk,:) = nodes(kk,:);
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function principalAxes = calculatePrincipalAxes(xyz, weights)
if nargin < 2; weights = ones(size(xyz,1),1); end;
% find the inertia tensor
inertiaTensor = zeros(3);
inertiaTensor(1,1) = sum(weights .* (xyz(:,2).^2 + xyz(:,3).^2)); inertiaTensor(2,2) = sum(weights .* (xyz(:,1).^2 + xyz(:,3).^2));
inertiaTensor(3,3) = sum(weights .* (xyz(:,1).^2 + xyz(:,2).^2)); inertiaTensor(1,2) = -sum(weights .* xyz(:,1) .* xyz(:,2));
inertiaTensor(1,3) = -sum(weights .* xyz(:,1) .* xyz(:,3)); inertiaTensor(2,3) = -sum(weights .* xyz(:,2) .* xyz(:,3));
inertiaTensor(2,1) = inertiaTensor(1,2); inertiaTensor(3,1) = inertiaTensor(1,3); inertiaTensor(3,2) = inertiaTensor(2,3);
% find the principal axes of the inertia tensor
[principalAxes, evMatrix] = eig(inertiaTensor);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tr = normalizeNeuron(tr)
% normalize the total dendritic length
lengthOnly = true;
totalLength = treeLength(tr,rootFinder(tr),lengthOnly);
for node = 1:numel(tr)
  if ~isempty(tr{node}{4}{2})
    tr{node}{4}{1} = tr{node}{4}{1}/totalLength; tr{node}{4}{2} = tr{node}{4}{2}/totalLength; tr{node}{4}{3} = tr{node}{4}{3}/totalLength;
  end
  tr{node}{3} = tr{node}{3}/totalLength;
end
% normalize the total dendritic volume by changing the cross-sectional areas, and keeping the lengths the same
lengthOnly = false;
totalLength = treeLength(tr,rootFinder(tr),lengthOnly);
for node = 1:numel(tr)
  if ~isempty(tr{node}{4}{4})
    tr{node}{4}{4} = tr{node}{4}{4}/totalLength;
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tree nodes areas rawLength] = integrateTree(tree,node,nodes,areas,rawLength,startingPos,interval)
tree{node}{3} = startingPos;
children = tree{node}{2};
for child = 1:numel(children)
  thisBranch = tree{children(child)}{4}{1}; thisBranch = thisBranch(end:-1:1,:);
  branchAreas = tree{children(child)}{4}{4}; branchAreas = branchAreas(end:-1:1);
  [thisBranch pathLengths branchAreas] = integrateBranch(thisBranch,branchAreas,startingPos,interval);
  tree{children(child)}{4}{1} = thisBranch(end:-1:1,:);
  nodes = [nodes; tree{children(child)}{4}{1}];
  tree{children(child)}{4}{2} = pathLengths(end:-1:1);
  tree{children(child)}{4}{4} = branchAreas;
  areas = [areas; tree{children(child)}{4}{4}];
  rawLength = rawLength+sum(pathLengths);
  [tree nodes areas rawLength] = integrateTree(tree,children(child),nodes,areas,rawLength,thisBranch(end,:),interval);
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [path pathLengths newAreas] = integrateBranch(branchCoordinates,areas,startingPos,interval)
path=startingPos; curPos=0; curEdgePos=0; nextPillar=2; curCoord = branchCoordinates(1,:);
edgeVectors = diff(branchCoordinates,1,1); edgeLengths=sqrt(sum(edgeVectors.^2,2)); totalLength=sum(edgeLengths);
curVec=edgeVectors(1,:); curVec = curVec/norm(curVec); nextPillarDistance=edgeLengths(1); newAreas = [];
while curPos<totalLength
  xyzSteps = abs((-curCoord+sqrt(curCoord.^2+2*curVec*interval))./(curVec+1e-6));
  step = min(max(xyzSteps),nextPillarDistance);
  path = [path; path(end,:)+curCoord*step+curVec*step^2/2];
  curEdgePos = curEdgePos+step;
  curCoord = curCoord+step*curVec;
  curPos = curPos+step;
  nextPillarDistance = nextPillarDistance-step;
  newAreas = [areas(min(nextPillar,size(edgeVectors,1))); newAreas]; % inverted orientation for areas to avoid inversion of array in integrateTree!
  if nextPillarDistance<totalLength/1e6 % 1e-6 precision
    if nextPillar>size(edgeVectors,1)
      break;
    end
    curCoord=branchCoordinates(nextPillar,:);
    curVec=edgeVectors(nextPillar,:); curVec=curVec/norm(curVec);
    nextPillarDistance = edgeLengths(nextPillar,:);
    nextPillar=nextPillar+1;
    curEdgePos=0;
  end
end
pathLengths = sqrt(sum(diff(path,1,1).^2,2));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tree = rotateTree(tree,rotMatrix)
for kk = 1:numel(tree)
  tree{kk}{3} = tree{kk}{3}*rotMatrix;
  if ~isempty(tree{kk}{4}{1})
    tree{kk}{4}{1} = tree{kk}{4}{1}*rotMatrix;
    tree{kk}{4}{3} = rotMatrix' * tree{kk}{4}{3};
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rotMatrix = findRotationMatrix(rotAxis,rotAngle)
ux=rotAxis(1); uy=rotAxis(2); uz=rotAxis(3);
rotMatrix = (cos(rotAngle)*eye(3) + sin(rotAngle)*[0 -uz uy;uz 0 -ux;-uy ux 0] + (1-cos(rotAngle))*kron(rotAxis,rotAxis'))';
if any(any(isnan(rotMatrix)))
  rotMatrix = eye(3);
end

