function [allIonECounts allIonlyCounts totalIonECount totalIonlyCount] = returnAllRatios_inhOnSpines(allTrees, minFeatureCount)

totalIonECount             = 0;
totalIonlyCount            = 0;
for tr = 1:numel(allTrees)
  allIonECounts{tr}        = [];
  allIonlyCounts{tr}       = [];
  for kk=2:numel(allTrees{tr})
    n1                     = nnz(allTrees{tr}{kk}{4}{5}(1:end-1,18-7)==0);
    n2                     = nnz(allTrees{tr}{kk}{4}{5}(1:end-1,18-7)==1);
    totalIonECount         = totalIonECount + n1;
    totalIonlyCount        = totalIonlyCount + n2;
    if n1+n2>=minFeatureCount
      allIonECounts{tr}    = [allIonECounts{tr} n1];
      allIonlyCounts{tr}   = [allIonlyCounts{tr} n2];
    end
  end
end
