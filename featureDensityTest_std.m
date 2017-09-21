function featureDensityTest_std(spineOrIS, basalOrApical, primaryOrIntermediateOrTerminal, qq, ISonShaftSet, minLength, shuffleCount)
%spineOrIS                       = 'spine';
%basalOrApical                   = 'apicalTuft';
%primaryOrIntermediateOrTerminal = 'intermediate'; % 'intermediate'; % 'terminal'; % 'intermediateOrTerminal'; % 'primaryOrIntermediateOrTerminal';
%qq                              = 0.0000001;
%ISonShaftSet                    = [-1 0 1];
%minLength                       = 0.1;
%shuffleCount                    = 1000;
load /rigel/stats/users/us2157/spines/data/shuffledTrees/trueTrees12.mat
[branchLengths eventCounts validBranches] = branchLengthsAndEventCounts(trees, spineOrIS, basalOrApical, primaryOrIntermediateOrTerminal, qq, ISonShaftSet, minLength);
allTrueRatios = eventCounts ./ branchLengths;
allTrueRatios(~validBranches) = -1;
for tr = 1:numel(trees)
  tmp                    = eventCounts(tr, :);
  totalEventCount(tr)    = sum(tmp(validBranches(tr, :)));
  tmp                    = branchLengths(tr, :);
  validBranchLengths{tr} = tmp(validBranches(tr, :));
  cumValidBrLengths{tr}  = cumsum(validBranchLengths{tr});
  totalValidBrLength(tr) = cumValidBrLengths{tr}(end);
  validBranchCounts(tr)  = nnz(validBranches(tr, :));
  sampleSTD(tr)          = std(allTrueRatios(tr, validBranches(tr, :)));
end
for tt = 1:shuffleCount
  for tr = 1:size(eventCounts, 1)
    shuffledPos         = totalValidBrLength(tr)*rand(1, totalEventCount(tr));
    prevCount           = 0;
    for br = 1:validBranchCounts(tr)
      newCount                  = nnz(shuffledPos<=cumValidBrLengths{tr}(br));
      shuffEventCounts(tr, br)  = newCount - prevCount;
      shuffEventRatios(tr, br)  = shuffEventCounts(tr, br) ./ validBranchLengths{tr}(br);
      prevCount                 = newCount;
    end
    shuffledSTD(tr, tt)         = std(shuffEventRatios(tr, 1:validBranchCounts(tr)));
  end
  allShuffles{tt}               = shuffEventRatios;
end
allOneSidedP=[]; all5quantile=[]; all95quantile=[]; allMean=[];
for tr=1:numel(trees)
  allOneSidedP(tr)  = nnz(sampleSTD(tr)<shuffledSTD(tr, :))/numel(shuffledSTD(tr, :));
  all5quantile(tr)  = quantile(shuffledSTD(tr,:), 0.05);
  all95quantile(tr) = quantile(shuffledSTD(tr,:), 0.95);
  allMean(tr)       = mean(shuffledSTD(tr,:));        
end
pooledTrueRatios = [];
pooledShufRatios = cell(1, shuffleCount);
for tr = 1:numel(trees)
  pooledTrueRatios = [pooledTrueRatios allTrueRatios(tr, validBranches(tr, :))];
end
pooledSampleSTD = std(pooledTrueRatios);
for tt = 1:shuffleCount
  for tr = 1:numel(trees)
    pooledShufRatios{tt} = [pooledShufRatios{tt} allShuffles{tt}(tr, 1:validBranchCounts(tr))];
  end
  pooledShuffSTD(tt) = std(pooledShufRatios{tt});
end
pooledOneSidedP  = nnz(pooledSampleSTD<pooledShuffSTD)/numel(pooledShuffSTD);
pooled5quantile  = quantile(pooledShuffSTD, 0.05);
pooled95quantile = quantile(pooledShuffSTD, 0.95);
pooledMean       = mean(pooledShuffSTD);
disp([spineOrIS ', ' basalOrApical ', ' primaryOrIntermediateOrTerminal ', min br. length:' num2str(minLength) ', neuron count with p<0.05: ' num2str(nnz(allOneSidedP<0.05)) ', pooled p-value: ' num2str(pooledOneSidedP) 'pool size: ' num2str(numel(pooledTrueRatios)) 'branches'])
disp(['0.05 value of the shuffled st. dev. for individual neurons: ' num2str(all5quantile)]);
disp(['0.95 value of the shuffled st. dev. for individual neurons: ' num2str(all95quantile)]);
disp(['mean value of the shuffled st. dev. for individual neurons: ' num2str(allMean)]);
disp(['true st. dev. values for individual neurons: ' num2str(sampleSTD)]);
%[counts,centers]=hist(allTrueRatios(1,validBranches(1,:)),10); figure;plot(centers,counts,'LineWidth',2); hold;
%for tt=1:shuffleCount; [counts,centers]=hist(allShuffles{tt}(1,1:validBranchCounts(1)),10); plot(centers,counts,'Color',[0.7 0.7 0.7]); end
