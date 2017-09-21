function [branchLengths eventCounts validBranches] = branchLengthsAndEventCounts(trees, spineOrIS, compartment, primaryOrIntermediateOrTerminal, qq, ISonShaftSet, minLength)

switch spineOrIS
  case 'spine'
    feature = 12 - 7;
  case 'IS'
    feature = 17 - 7;
  otherwise
    warning('Unexpected spine-or-IS input')
end
switch compartment
  case 'basal'
    dendriteType = 7;
  case 'apical'
    dendriteType = [3 4 5];
  case 'apicalTuft'
    dendriteType = 3;
  case 'apicalTrunk'
    dendriteType = 4;
  case 'apicalOblique'
    dendriteType = 5;
  case 'basalOrApical'
    dendriteType = [3 4 5 7];
  otherwise
    warning('Unexpected dendrite type')
end
maxBranchCount                = max(cellfun(@numel, trees));
branchLengths                 = -1 * ones(numel(trees), maxBranchCount);
eventCounts                   = -1 * ones(numel(trees), maxBranchCount);
validBranches                 = false(numel(trees), maxBranchCount);
for tr = 1:numel(trees)
  nodesOfInterestVolumes      = [];
  for kk=2:numel(trees{tr})
    nodesOfInterest           = find(trees{tr}{kk}{4}{5}(1:end-1,feature)>0 & ismember(trees{tr}{kk}{4}{5}(1:end-1,19-7), dendriteType) & ismember(trees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
    nodesOfInterestVolumes    = [nodesOfInterestVolumes; trees{tr}{kk}{4}{5}(nodesOfInterest,feature)];
  end
  th                          = quantile(nodesOfInterestVolumes, qq);
  for br=2:numel(trees{tr})
    denTyp                    = trees{tr}{br}{4}{5}(1:end-1,19-7);
    if trees{tr}{br}{1}==1
      parentDenTyp            = 0;
    else
      parentDenTyp            = trees{tr}{trees{tr}{br}{1}}{4}{5}(1:end-1,19-7);
    end
    hasChildren               = ~isempty(trees{tr}{br}{2});
    [flag, lastNode]          = primaryOrIntermediateOrTerminal_function(denTyp, parentDenTyp, hasChildren, dendriteType, primaryOrIntermediateOrTerminal);    
    branchLength              = abs(diff(trees{tr}{br}{4}{5}([1 lastNode], 11-7)));
    if flag & (branchLength>=minLength)
      largeNodesOfInterest    = find(trees{tr}{br}{4}{5}(1:end-1,feature)>th & ismember(denTyp, dendriteType) & ismember(trees{tr}{br}{4}{5}(1:end-1,18-7), ISonShaftSet));
      branchLengths(tr, br)   = branchLength;
      eventCounts(tr, br)     = numel(largeNodesOfInterest);
      validBranches(tr, br)   = true;
    end
  end
end
