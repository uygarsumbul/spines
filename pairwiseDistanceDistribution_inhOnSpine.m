function [counts,centers, allCounts, allPairwiseDistances] = pairwiseDistanceDistribution(allTrees, thisFeature, maxmax, binSize)

%code{1}{1}                               = '9_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{2}{1}                               = '18_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{3}{1}                               = '19_unsampled_reconstruction_stitched_connected_labelled.eswc';
%code{4}{1}                               = '20_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{5}{1}                               = '21_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{6}{1}                               = '22_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{7}{1}                               = '23_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{8}{1}                               = '24_upsampled_reconstruction_stitched_connected_labelled.eswc';
%code{9}{1}                               = '25_upsampled_reconstruction_stitched_connected_labelled.eswc';
%directory                                = '/home/uygar/spines/data/upsampledFinalReconstructions/';
%options.resolution                       = [0.12 0.12 0.1];
%allTrees                                 = readDataset_spines(code,directory,options);

allPairwiseDistances                     = cell(0);
for tr = 1:numel(allTrees)
  pairwiseDistances                      = [];
  for kk=2:numel(allTrees{tr})
    parent                               = allTrees{tr}{kk}{1};
    parentNotRoot                        = (parent~=1);
      if parentNotRoot
        largeNodesOfInterest             = find(allTrees{tr}{kk}{4}{5}(1:end-1,thisFeature-7)==0);
        firstParentalLargeNodeOfInterest = find(allTrees{tr}{parent}{4}{5}(1:end-1,thisFeature-7)==0, 1);
      else
        largeNodesOfInterest             = find(allTrees{tr}{kk}{4}{5}(:,thisFeature-7)==0);
      end
    localDistances                       = cumsum([0; allTrees{tr}{kk}{4}{2}]);
    if parentNotRoot
      localDistancesParent               = cumsum([0; allTrees{tr}{parent}{4}{2}]) + localDistances(end);
      pairwiseDistances                  = [pairwiseDistances; diff([localDistances(largeNodesOfInterest); localDistancesParent(firstParentalLargeNodeOfInterest)])];
    else
      pairwiseDistances                  = [pairwiseDistances; diff(localDistances(largeNodesOfInterest))];
    end
  end
  allPairwiseDistances{tr}             = pairwiseDistances;
end
[counts,centers]                       = hist(cell2mat(allPairwiseDistances'), binSize/2:binSize:maxmax-binSize/2);
allCounts                              = cell(1, tr);
for tr = 1:numel(allTrees)
  [tmp, ~]                             = hist(allPairwiseDistances{tr}, binSize/2:binSize:maxmax-binSize/2);
  allCounts{tr}                        = tmp;
end
