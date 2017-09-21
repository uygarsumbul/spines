function [allPairwiseDistances, allXYZofNodesOfInterest, allSomaDistancesOfNodesOfInterest, pairAndSomaDistances] = pairwiseDistanceDistribution(trees, spineOrIS, ISonShaftSet, qq, basalOrApical, primalOrIntermediateOrTerminal)

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
    dendriteType = 7;
  case 'apical'
    dendriteType = 4;
  case 'basalAndApical'
    dendriteType = [4 7];
  otherwise
    warning('Unexpected dendrite type')
end

allPairwiseDistances                               = cell(0);
for tr = 1:numel(trees)
  nodesOfInterestVolumes                           = [];
  for kk=2:numel(trees{tr})
    nodesOfInterest                                = find(trees{tr}{kk}{4}{5}(1:end-1,feature)>0 & ismember(trees{tr}{kk}{4}{5}(1:end-1,19-7), dendriteType) & ismember(trees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
    nodesOfInterestVolumes                         = [nodesOfInterestVolumes; trees{tr}{kk}{4}{5}(nodesOfInterest,feature)];
  end
  th                                               = quantile(nodesOfInterestVolumes, qq);
  for kk=2:numel(trees{tr})
    flag = false;
    switch primalOrIntermediateOrTerminal
      case 'primal'
        flag = isempty(trees{tr}{kk}{1});
      case 'intermediate'
        flag = ~isempty(trees{tr}{kk}{1}) & ~isempty(trees{tr}{kk}{2});
      case 'terminal'
        flag = isempty(trees{tr}{kk}{2});
      case 'intermediateOrTerminal'
        flag = ~isempty(trees{tr}{kk}{1}); 
      case 'primalOrIntermediateOrTerminal'
        flag = true;
      otherwise
        warning('Unexpected branch type');
    end
    if flag
      largeNodesOfInterest{kk}                       = find(trees{tr}{kk}{4}{5}(1:end-1,feature)>th & ismember(trees{tr}{kk}{4}{5}(1:end-1,19-7), dendriteType) & ismember(trees{tr}{kk}{4}{5}(1:end-1,18-7), ISonShaftSet));
      lastNodes{kk}                                  = trees{tr}{kk}{4}{5}(end,feature)>th & ismember(trees{tr}{kk}{4}{5}(end,19-7), dendriteType) & ismember(trees{tr}{kk}{4}{5}(end,18-7), ISonShaftSet);
    else
      largeNodesOfInterest{kk}                       = [];
      lastNodes{kk}                                  = false;
    end
  end
  for kk=2:numel(trees{tr})
    if lastNodes{kk}
      largeNodesOfInterest{kk}                     = [largeNodesOfInterest{kk}; size(trees{tr}{kk}{4}{5}, 1)];
      allPairwiseDistances{tr}{kk}                 = diff(-trees{tr}{kk}{4}{5}(largeNodesOfInterest{kk}, 11-7));
      allXYZofNodesOfInterest{tr}{kk}              = trees{tr}{kk}{4}{1}(largeNodesOfInterest{kk}, :);
      allSomaDistancesOfNodesOfInterest{tr}{kk}    = trees{tr}{kk}{4}{5}(largeNodesOfInterest{kk}, 11-7);
    else
      ancestor                                     = trees{tr}{kk}{1};
      firstAncestralLargeNodeOfInterest            = [];
      firstAncestorOfInterest                      = [];
      while (ancestor ~= 1) & isempty(firstAncestralLargeNodeOfInterest)
        if ~isempty(largeNodesOfInterest{ancestor})
          firstAncestralLargeNodeOfInterest        = largeNodesOfInterest{ancestor}(1);
          firstAncestorOfInterest                  = ancestor;
        end
        ancestor                                   = trees{tr}{ancestor}{1};
      end
      if isempty(firstAncestralLargeNodeOfInterest)
        allPairwiseDistances{tr}{kk}               = diff(-trees{tr}{kk}{4}{5}(largeNodesOfInterest{kk}, 11-7));
        allXYZofNodesOfInterest{tr}{kk}            = trees{tr}{kk}{4}{1}(largeNodesOfInterest{kk}, :);
        allSomaDistancesOfNodesOfInterest{tr}{kk}  = trees{tr}{kk}{4}{5}(largeNodesOfInterest{kk}, 11-7);
      else
        negDist                                    = -[trees{tr}{kk}{4}{5}(largeNodesOfInterest{kk}, 11-7); trees{tr}{firstAncestorOfInterest}{4}{5}(firstAncestralLargeNodeOfInterest, 11-7)];
        allPairwiseDistances{tr}{kk}               = diff(negDist);

        allXYZofNodesOfInterest{tr}{kk}            = [trees{tr}{kk}{4}{1}(largeNodesOfInterest{kk}, :); trees{tr}{firstAncestorOfInterest}{4}{1}(firstAncestralLargeNodeOfInterest, :)];
        allSomaDistancesOfNodesOfInterest{tr}{kk}  = [trees{tr}{kk}{4}{5}(largeNodesOfInterest{kk}, 11-7); trees{tr}{firstAncestorOfInterest}{4}{5}(firstAncestralLargeNodeOfInterest, 11-7)];
      end
    end
  end
end

pairAndSomaDistances = cell(0);
for tr = 1:numel(trees)
  pairAndSomaDistances{tr} = [];
  for kk = 2:numel(trees{tr})
    if numel(allPairwiseDistances{tr}{kk})>0

      bo=1; br=kk; while true; br=trees{tr}{br}{1}; if br~=1; bo=bo+1; else; break; end; end;
      pairAndSomaDistances{tr} = [pairAndSomaDistances{tr}; [allPairwiseDistances{tr}{kk} (allSomaDistancesOfNodesOfInterest{tr}{kk}(1:end-1)+allSomaDistancesOfNodesOfInterest{tr}{kk}(2:end))/2 kk*ones(size(allPairwiseDistances{tr}{kk})) bo*ones(size(allPairwiseDistances{tr}{kk}))]];
    end
  end
end
%[counts,centers]                       = hist(cell2mat(allPairwiseDistances'), binSize/2:binSize:maxmax-binSize/2);
