function [allPairwiseDistances, allXYZofNodesOfInterest, allSomaDistancesOfNodesOfInterest, pairAndSomaDistances] = pairwiseDistanceDistribution(allTrees, spineOrIS, ISonShaftSet, qq, basalOrApical)

switch spineOrIS
  case 'spine'
    feature = 12 - 7;
  case 'IS'
    feature = 17 - 7;
  otherwise
    warning('Unexpected spine-or-IS input')
end
switch basalOrApical
  case 'basal'
    branchType = 7;
  case 'apical'
    branchType = 4;
  case 'basalAndApical'
    branchType = [4 7];
  otherwise
    warning('Unexpected branch type')
end
allPairwiseDistances                             = cell(0);
for tr = 1:numel(allTrees)
  nodesOfInterestVolumes                         = [];
  for kk=2:numel(allTrees{tr})
    nodesOfInterest                              = find(allTrees{tr}{kk}{4}{5}(1:end-1,feature)>0 & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
    nodesOfInterestVolumes                       = [nodesOfInterestVolumes; allTrees{tr}{kk}{4}{5}(nodesOfInterest,feature)];
  end
  threshold                                      = quantile(nodesOfInterestVolumes, qq);
  for kk=2:numel(allTrees{tr})
    largeNodesOfInterest{kk}                     = find(allTrees{tr}{kk}{4}{5}(1:end-1,feature)>threshold & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
    lastNodes{kk}                                = allTrees{tr}{kk}{4}{5}(end,feature)>threshold & ismember(allTrees{tr}{kk}{4}{5}(end,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(end,18-7), ISonShaftSet);
  end
  for kk=2:numel(allTrees{tr})
    if lastNodes{kk}
      largeNodesOfInterest{kk}                   = [largeNodesOfInterest; size(allTrees{tr}{kk}{4}{5}, 1)];
      localDistances                             = cumsum([0; allTrees{tr}{kk}{4}{2}]);
      allPairwiseDistances{tr}{kk}               = diff(localDistances(largeNodesOfInterest));
      allXYZofNodesOfInterest{tr}{kk}            = allTrees{tr}{kk}{4}{1}(largeNodesOfInterest, :);
      allSomaDistancesOfNodesOfInterest{tr}{kk}  = allTrees{tr}{kk}{4}{5}(largeNodesOfInterest, 11-7);
    else
      ancestor                                   = allTrees{tr}{kk}{1};
      firstAncestralLargeNodeOfInterest          = [];
      firstAncestorOfInterest                    = [];
      while (ancestor ~= 1) & isempty(firstAncestralLargeNodeOfInterest)
        firstAncestralLargeNodeOfInterest        = find(allTrees{tr}{ancestor}{4}{5}(:,feature)>threshold & ismember(allTrees{tr}{ancestor}{4}{5}(:,19-7), branchType) & ismember(allTrees{tr}{ancestor}{4}{5}(:,18-7), ISonShaftSet), 1);
        if ~isempty(firstAncestralLargeNodeOfInterest)
          firstAncestorOfInterest                = ancestor;
        end
        ancestor                                 = allTrees{tr}{kk}{1};
      end
      if isempty(firstAncestralLargeNodeOfInterest)
	largeNodesOfInterest{kk}                   = [largeNodesOfInterest; size(allTrees{tr}{kk}{4}{5}, 1)];
        localDistances                             = cumsum([0; allTrees{tr}{kk}{4}{2}]);
        allPairwiseDistances{tr}{kk}               = diff(localDistances(largeNodesOfInterest));
        allXYZofNodesOfInterest{tr}{kk}            = allTrees{tr}{kk}{4}{1}(largeNodesOfInterest, :);
        allSomaDistancesOfNodesOfInterest{tr}{kk}  = allTrees{tr}{kk}{4}{5}(largeNodesOfInterest, 11-7);
      else
        localDistances                             = cumsum([0; allTrees{tr}{kk}{4}{2}]);
        offset                                     = allTrees{tr}{kk}{4}{5}(1, 11-7) - allTrees{tr}{firstAncestorOfInterest}{4}{5}(1, 11-7);
        localDistancesAncestor                     = cumsum([0; allTrees{tr}{firstAncestorOfInterest}{4}{2}]) + offset;
        allPairwiseDistances{tr}{kk}               = diff([localDistances(largeNodesOfInterest); localDistancesAncestor(firstAncestralLargeNodeOfInterest)]);
        allXYZofNodesOfInterest{tr}{kk}            = [allTrees{tr}{kk}{4}{1}(largeNodesOfInterest, :); allTrees{tr}{firstAncestorOfInterest}{4}{1}(firstAncestralLargeNodeOfInterest, :)];
        allSomaDistancesOfNodesOfInterest{tr}{kk}  = [allTrees{tr}{kk}{4}{5}(largeNodesOfInterest, 11-7); allTrees{tr}{firstAncestorOfInterest}{4}{5}(firstAncestralLargeNodeOfInterest,, 11-7)];
      end
    end
  end
end

pairAndSomaDistances = cell(0);
for tr = 1:numel(allTrees)
  pairAndSomaDistances{tr} = [];
  for kk = 2:numel(allTrees{tr})
    if numel(allPairwiseDistances{tr}{kk})>0
      pairAndSomaDistances{tr} = [pairAndSomaDistances{tr}; [allPairwiseDistances{tr}{kk} (allSomaDistancesOfNodesOfInterest{tr}{kk}(1:end-1)+allSomaDistancesOfNodesOfInterest{tr}{kk}(2:end))/2 kk*ones(size(allPairwiseDistances{tr}{kk}))]];
    end
  end
end
%[counts,centers]                       = hist(cell2mat(allPairwiseDistances'), binSize/2:binSize:maxmax-binSize/2);
