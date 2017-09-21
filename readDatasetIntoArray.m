function [allTrees, allPAs] = readDatasetIntoArray(code,directory,options)

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
  [nodes,edges,radii,nodeTypes,features,~] = readSWCfile(thisFileName,neuronParts,resolution);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nodes,edges,radii,nodeType,features,abort] = readSWCfile(fileName,validNodeTypes,resolution)
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
%nodeTypes = unique(nodeType)';
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
% swap on node types
tmpT = nodeType(1); nodeType(1) = nodeType(somaNode,:); nodeType(somaNode,:) = tmpT;
