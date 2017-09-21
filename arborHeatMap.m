function [allMu, mini,maxi] = arborHeatMap(tree, pairAndSomaDistances, windowSize, drawPlot)

allMu = cell(1,numel(tree));
cd /rigel/stats/users/us2157/spines/repo/treeFunctions/
for kk = 2:numel(tree)
  des = allDescendents(tree, kk);
  anc = allAncestors(tree, kk);
  plusBranches = [des kk anc];
  allMu{kk} = zeros(1, size(tree{kk}{4}{5},1));
  for mm = 1:size(tree{kk}{4}{5},1)
    somaDist = tree{kk}{4}{5}(mm, 11-7);
    relevant = find(ismember(pairAndSomaDistances(:,3), plusBranches) & pairAndSomaDistances(:,2)>=somaDist-windowSize/2 & pairAndSomaDistances(:,2)<somaDist+windowSize/2);
    relevantPairDistances = pairAndSomaDistances(relevant, 1);
    if ~isempty(relevantPairDistances)
      allMu{kk}(mm) = 1/(mean(relevantPairDistances)+eps); %numel(relevantPairDistances) / (sum(relevantPairDistances) + eps);
    else
      allMu{kk}(mm) = 0;
    end
  end
end
cd /rigel/stats/users/us2157/spines/repo/

allMeans = cell2mat(allMu); maxi = quantile(allMeans, 0.97); mini = quantile(allMeans, 0.03);
if drawPlot
  figure;hold;
  for kk = 2:numel(tree)
    for mm = 1:size(tree{kk}{4}{1},1)
      tt = tree{kk}{4}{1}(mm, :);
      ratio = min(1, max(0, (allMu{kk}(mm)-mini)/(maxi-mini)));
      plot3(tt(:,1),tt(:,2),tt(:,3),'Color',[0 1 0]*ratio+[1 0 1]*(1-ratio),'MarkerSize',5,'Marker','.');
    end
  end
end
