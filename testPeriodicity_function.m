function [counts,centers, allLSPD] = testPeriodicity_function(thisFeature, LARGER, qq, maxmax, binSize)
cd ~/spines/repo

code{1}{1}  = '9_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{2}{1}  = '18_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{3}{1}  = '19_unsampled_reconstruction_stitched_connected_labelled.eswc';
code{4}{1}  = '20_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{5}{1}  = '21_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{6}{1}  = '22_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{7}{1}  = '23_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{8}{1}  = '24_upsampled_reconstruction_stitched_connected_labelled.eswc';
code{9}{1}  = '25_upsampled_reconstruction_stitched_connected_labelled.eswc';
directory          = '/home/uygar/spines/data/upsampledFinalReconstructions/';

%code{1}{1}  = 'AIBS_19_AllSynapses_stitched_connected_pruned_labelled.eswc';
%code{2}{1}  = 'AIBS_21_AllSynapses_stitched_connected_pruned_labelled.eswc';
%code{3}{1}  = 'AIBS_25_AllSynapses_stitched_connected_pruned_labelled.eswc';
%code{4}{1}  = 'AIBS_20_AllSynapses_stitched_connected_pruned_labelled.eswc';
%code{5}{1}  = 'AIBS_22_AllSynapses_stitched_connected_pruned_labelled.eswc';
%directory   = 'data/totalReconstructions/';

options.resolution = 1*[0.12 0.12 0.1];
allTrees           = readDataset_spines(code,directory,options);

myxx = 0:binSize:maxmax; myxx=myxx-round(maxmax/2); myxx=myxx/(round(maxmax/2)/2.5);

allLSPD                              = cell(0);
for tr = 1:numel(allTrees)
  spineVolumes                       = []; % SHOULD WE AVERAGE OVER ALL THE NEURONS?
  for kk=2:numel(allTrees{tr})
    theseSpines                      = find(allTrees{tr}{kk}{4}{5}(1:end-1,thisFeature-7)>0);
    spineVolumes                     = [spineVolumes; allTrees{tr}{kk}{4}{5}(theseSpines,thisFeature-7)];
  end
%  largeSpineThreshold                = th*mean(spineVolumes);
  largeSpineThreshold                = quantile(spineVolumes, qq); disp(['fraction of mean size: ' num2str(largeSpineThreshold/mean(spineVolumes))])
  largeSpinePairDistances            = [];
  for kk=2:numel(allTrees{tr})
    if LARGER
      largeSpines                    = find(allTrees{tr}{kk}{4}{5}(1:end-1,thisFeature-7)>largeSpineThreshold);
    else
      largeSpines                    = find(allTrees{tr}{kk}{4}{5}(1:end-1,thisFeature-7)<=largeSpineThreshold & allTrees{tr}{kk}{4}{5}(1:end-1,thisFeature-7)>0);
    end

    temp = cumsum([0; allTrees{tr}{kk}{4}{2}]);
    largeSpinePairDistances          = [largeSpinePairDistances; diff(temp(largeSpines))];
  end
  allLSPD{tr} = largeSpinePairDistances;
end

[counts,centers]=hist(cell2mat(allLSPD'), 0:binSize:maxmax);
