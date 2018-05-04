function eiBalanceHeatMap_6march2018

saturate = 0.01; disp(saturate)
hiCol = [0 1 1];
loCol = [1 0 0];
maSi  = 10;

load /home/uygar/misc/columbiaOnHabanero/spines/data/shuffledTrees/trueTrees12.mat
filePreamble                             = '/home/uygar/misc/columbiaOnHabanero/spines/results/heatMaps/april27_newColorMap_eiBalanceHeatMap_fullRange_';
denType                                  = 'basalAndApical';
brType                                   = 'primaryOrIntermediateOrTerminal';
cd /home/uygar/misc/columbiaOnHabanero/spines/repo/treeFunctions/
for tr = 1:numel(trees)
  trLen(tr) = treeLength(trees{tr},1,true);
end
cd /home/uygar/misc/columbiaOnHabanero/spines/repo/
synType1{1}     = 'spine'; synType2{1}     = 'IS';    IStype1{1}     = [-1 0 1]; IStype2{1}     = [0 1];    qqType1{1}     = 0.00001; qqType2{1}     = 0.00001; % spine, IS
synType1{end+1} = 'spine'; synType2{end+1} = 'IS';    IStype1{end+1} = [-1 0 1]; IStype2{end+1} = [0];      qqType1{end+1} = 0.00001; qqType2{end+1} = 0.00001; % spine, IS on spine
synType1{end+1} = 'spine'; synType2{end+1} = 'IS';    IStype1{end+1} = [-1 0 1]; IStype2{end+1} = [1];      qqType1{end+1} = 0.00001; qqType2{end+1} = 0.00001; % spine, IS on shaft

for kk = 1:numel(synType1)
  [~, ~, ~, pairAndSomaDistances1] = extractPairwiseDistances(trees, synType1{kk}, IStype1{kk}, qqType1{kk}, denType, brType);
  ISstr1 = num2str(IStype1{kk}); ISstr1(ISstr1==' ') = '';
  [~, ~, ~, pairAndSomaDistances2] = extractPairwiseDistances(trees, synType2{kk}, IStype2{kk}, qqType2{kk}, denType, brType);
  ISstr2 = num2str(IStype2{kk}); ISstr2(ISstr2==' ') = '';
  for tr = 1:numel(trees)

    if ~exist([filePreamble synType1{kk} '_' synType2{kk} '__ISonShaft' ISstr1 '_' ISstr2 '_neuron' num2str(tr) '_positiveLambda_absDiff.pdf'])
	     
    lambda1 = size(pairAndSomaDistances1{tr},1)/trLen(tr);
    lambda2 = size(pairAndSomaDistances2{tr},1)/trLen(tr);
    miniX = 99999;
    miniY = 99999;
    if lambda1 > 0 & lambda2 > 0
      [allMu1, ~, ~] = arborHeatMap(trees{tr}, pairAndSomaDistances1{tr}, 16/lambda1, false);
      [allMu2, ~, ~] = arborHeatMap(trees{tr}, pairAndSomaDistances2{tr}, 16/lambda2, false);
      allMeans1 = cell2mat(allMu1); maxi1 = quantile(allMeans1, 1-saturate); mini1 = quantile(allMeans1, saturate);
      allMeans2 = cell2mat(allMu2); maxi2 = quantile(allMeans2, 1-saturate); mini2 = quantile(allMeans2, saturate);
      figure;hold;
      for bb = 2:numel(trees{tr})
	for mm = 1:size(trees{tr}{bb}{4}{1},1)
	  tt = trees{tr}{bb}{4}{1}(mm, :); miniX = min(miniX, min(tt(1))); miniY = min(miniY, min(tt(2)));
          ratio1 = min(1, max(0, (allMu1{bb}(mm)-mini1)/(maxi1-mini1)));
          ratio2 = min(1, max(0, (allMu2{bb}(mm)-mini2)/(maxi2-mini2)));
          absDiff = abs(ratio1 - ratio2);
          plot3(tt(:,1),tt(:,2),tt(:,3),'Color',hiCol*absDiff+loCol*(1-absDiff),'MarkerSize',maSi,'Marker','.');
        end
      end
      text(miniX-100, miniY-40,'0','FontSize',16); text(miniX-20, miniY-40,'1','FontSize',16);
      for tt=0:99; plot3([miniX-100+tt miniX-100+tt+1], [miniY-20 miniY-20], [0 0],'LineWidth',10,'Color', hiCol*tt/99+loCol*(1-tt/99)); end;
      view(0, 90); set(gcf,'Color','w','Position',[100 100 800 800]); axis off; axis square;
      saveas(gcf,[filePreamble synType1{kk} '_' synType2{kk} '__ISonShaft' ISstr1 '_' ISstr2 '_neuron' num2str(tr) '_positiveLambda_absDiff.pdf']); close;
    end
    end
  end
end

