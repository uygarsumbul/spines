function [allPairwiseDistances, allXYZofNodesOfInterest, allSomaDistancesOfNodesOfInterest, pairAndSomaDistances] = pairwiseDistanceDistribution(allTrees, spineOrIS, ISonShaftSet, LARGER, qq, basalOrApical)

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
    parent                                       = allTrees{tr}{kk}{1};
    parentNotRoot                                = (parent~=1);
    if LARGER
      if parentNotRoot
        largeNodesOfInterest                     = find(allTrees{tr}{kk}{4}{5}(1:end-1,feature)>threshold & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
        firstParentalLargeNodeOfInterest         = find(allTrees{tr}{parent}{4}{5}(1:end-1,feature)>threshold & ismember(allTrees{tr}{parent}{4}{5}(1:end-1,19-7), branchType) & ismember(allTrees{tr}{parent}{4}{5}(1:end-1,18-7), ISonShaftSet), 1);
      else
        largeNodesOfInterest                     = find(allTrees{tr}{kk}{4}{5}(:,feature)>threshold & ismember(allTrees{tr}{kk}{4}{5}(:,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(:,18-7), ISonShaftSet));
      end
    else
      if parentNotRoot
        tempSizes                                = allTrees{tr}{kk}{4}{5}(1:end-1,feature);
        largeNodesOfInterest                     = find(tempSizes<=threshold & tempSizes>0 & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
        tempSizes                                = allTrees{tr}{parent}{4}{5}(1:end-1,feature);
        firstParentalLargeNodeOfInterest         = find(tempSizes<=threshold & tempSizes>0 & ismember(allTrees{tr}{parent}{4}{5}(1:end-1,19-7), branchType) & ismember(allTrees{tr}{parent}{4}{5}(1:end-1,18-7), ISonShaftSet), 1);
      else
        tempSizes                                = allTrees{tr}{kk}{4}{5}(:,feature);
        largeNodesOfInterest                     = find(tempSizes<=threshold & tempSizes>0 & ismember(allTrees{tr}{kk}{4}{5}(:,19-7), branchType) & ismember(allTrees{tr}{kk}{4}{5}(:,18-7), ISonShaftSet));
      end
    end
    localDistances                               = cumsum([0; allTrees{tr}{kk}{4}{2}]);
    if parentNotRoot
      localDistancesParent                       = cumsum([0; allTrees{tr}{parent}{4}{2}]) + localDistances(end);
      allPairwiseDistances{tr}{kk}               = diff([localDistances(largeNodesOfInterest); localDistancesParent(firstParentalLargeNodeOfInterest)]);
      allXYZofNodesOfInterest{tr}{kk}            = [allTrees{tr}{kk}{4}{1}(largeNodesOfInterest, :); allTrees{tr}{parent}{4}{1}(firstParentalLargeNodeOfInterest, :)];
      allSomaDistancesOfNodesOfInterest{tr}{kk}  = [allTrees{tr}{kk}{4}{5}(largeNodesOfInterest, 11-7); allTrees{tr}{parent}{4}{5}(firstParentalLargeNodeOfInterest, 11-7)];
    else
      allPairwiseDistances{tr}{kk}               = diff(localDistances(largeNodesOfInterest));
      allXYZofNodesOfInterest{tr}{kk}            = allTrees{tr}{kk}{4}{1}(largeNodesOfInterest, :);
      allSomaDistancesOfNodesOfInterest{tr}{kk}  = allTrees{tr}{kk}{4}{5}(largeNodesOfInterest, 11-7);
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
