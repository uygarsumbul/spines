function [allMu, mini,maxi] = arborHeatMap_proportionOfTwoFeatures(tree, pairAndSomaDistances1, pairAndSomaDistances2, windowSize, drawPlot)

allMu = cell(1,numel(tree));
for kk = 2:numel(tree)
  des = allDescendents(tree, kk);
  anc = allAncestors(tree, kk);
  plusBranches = [des kk anc];
  allMu{kk} = zeros(1, size(tree{kk}{4}{5},1));
  for mm = 1:size(tree{kk}{4}{5},1)
    somaDist = tree{kk}{4}{5}(mm, 11-7);
    relevant1 = find(ismember(pairAndSomaDistances1(:,3), plusBranches) & pairAndSomaDistances1(:,2)>=somaDist-windowSize/2 & pairAndSomaDistances1(:,2)<somaDist+windowSize/2);
    relevantPairDistances1 = pairAndSomaDistances1(relevant1, 1);
    relevant2 = find(ismember(pairAndSomaDistances2(:,3), plusBranches) & pairAndSomaDistances2(:,2)>=somaDist-windowSize/2 & pairAndSomaDistances2(:,2)<somaDist+windowSize/2);
    relevantPairDistances2 = pairAndSomaDistances2(relevant2, 1);
    allMu{kk}(mm) = (numel(relevantPairDistances1) / (sum(relevantPairDistances1) + eps)) / (numel(relevantPairDistances2) / (sum(relevantPairDistances2) + eps) + eps);

if allMu{kk}(mm)>1; disp([kk mm]); end;

  end
end

allMeans = cell2mat(allMu); maxi = quantile(allMeans, 0.98); mini = quantile(allMeans, 0.02);
if drawPlot
  figure;hold;
  for kk = 2:numel(tree)
    for mm = 1:size(tree{kk}{4}{1},1)
      tt = tree{kk}{4}{1}(mm, :);
      ratio = min(1, max(0, (allMu{kk}(mm)-mini)/(maxi-mini)));
      plot3(tt(:,1),tt(:,2),tt(:,3),'Color',[1 0 0]*ratio+[0 0 1]*(1-ratio),'MarkerSize',5,'Marker','.');
    end
  end
end
