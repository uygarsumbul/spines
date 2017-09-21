function allHeatMapsAndSpecialNodes_function
load /rigel/stats/users/us2157/spines/data/shuffledTrees/trueTrees12.mat
filePreamble                             = '/rigel/stats/users/us2157/spines/results/heatMaps/jul28_3Percent_heatMapsAndSpecialNodes_';
synType                                  = {'spine', 'IS'};
IStype                                   = {[-1 0 1], 0, 1, [0 1]};
qqType                                   = {0.00001, 0.8};
denType                                  = 'basalAndApical';
brType                                   = 'primaryOrIntermediateOrTerminal';
cd /rigel/stats/users/us2157/spines/repo/treeFunctions/
for tr = 1:numel(trees)
  trLen(tr) = treeLength(trees{tr},1,true);
end
cd /rigel/stats/users/us2157/spines/repo/
synType1{1}     = 'IS';    synType2{1}     = 'IS';    IStype1{1}    = [0 1];     IStype2{1}     = 0;        qqType1{1}     = 0.00001; qqType2{1}     = 0.00001; % IS, punkta on spine
synType1{end+1} = 'spine'; synType2{end+1} = 'spine'; IStype1{end+1} = [-1 0 1]; IStype2{end+1} = [-1 0 1]; qqType1{end+1} = 0.00001; qqType2{end+1} = 0.8;     % spine, punkta large spine
synType1{end+1} = 'IS';    synType2{end+1} = 'IS';    IStype1{end+1} = [0 1];    IStype2{end+1} = [0 1];    qqType1{end+1} = 0.00001; qqType2{end+1} = 0.8;     % IS, punkta large spine
synType1{end+1} = 'IS';    synType2{end+1} = [];      IStype1{end+1} = [0 1];    IStype2{end+1} = [];       qqType1{end+1} = 0.00001; qqType2{end+1} = [];      % IS
synType1{end+1} = 'IS';    synType2{end+1} = [];      IStype1{end+1} = [0 1];    IStype2{end+1} = [];       qqType1{end+1} = 0.8;     qqType2{end+1} = [];      % large IS
synType1{end+1} = 'IS';    synType2{end+1} = 'IS';    IStype1{end+1} = [0 1];    IStype2{end+1} = [0 1];    qqType1{end+1} = 0.8;     qqType2{end+1} = 0.8;     % large IS, punkta large IS
synType1{end+1} = 'spine'; synType2{end+1} = [];      IStype1{end+1} = [-1 0 1]; IStype2{end+1} = [];       qqType1{end+1} = 0.00001; qqType2{end+1} = [];      % spine
synType1{end+1} = 'spine'; synType2{end+1} = [];      IStype1{end+1} = [-1 0 1]; IStype2{end+1} = [];       qqType1{end+1} = 0.8;     qqType2{end+1} = [];      % large spine
synType1{end+1} = 'spine'; synType2{end+1} = 'spine'; IStype1{end+1} = [-1 0 1]; IStype2{end+1} = [-1 0 1]; qqType1{end+1} = 0.8;     qqType2{end+1} = 0.8;     % large spine, punkta large spine


for kk = 1:numel(synType1)
  [~, ~                    , ~, pairAndSomaDistances1] = extractPairwiseDistances(trees, synType1{kk}, IStype1{kk}, qqType1{kk}, denType, brType);
  ISstr1 = num2str(IStype1{kk}); ISstr1(ISstr1==' ') = '';
  if ~isempty(synType2{kk})
    [~, xyzOfNodesOfInterest2, ~, pairAndSomaDistances2] = extractPairwiseDistances(trees, synType2{kk}, IStype2{kk}, qqType2{kk}, denType, brType);
    ISstr2 = num2str(IStype2{kk}); ISstr2(ISstr2==' ') = '';
  end
  for tr = 1:numel(trees)
    lambda = size(pairAndSomaDistances1{tr},1)/trLen(tr);
    if lambda > 16/100 %70
      [allMu, mini, maxi] = arborHeatMap(trees{tr}, pairAndSomaDistances1{tr}, 16/lambda, false);
      if ~isempty(synType2{kk})
        plotHeatMapAndSpecialNodes3(trees{tr}, allMu, cell2mat(xyzOfNodesOfInterest2{tr}'));
        saveas(gcf,[filePreamble synType1{kk} '_' synType2{kk} '__ISonShaft' ISstr1 '_' ISstr2 '__quartile' num2str(qqType1{kk}) '_' num2str(qqType2{kk}) '_neuron' num2str(tr) '_maxWidth100_bigPunkta.jpg']); close;
      else
        plotHeatMapAndSpecialNodes3(trees{tr}, allMu, []);
	saveas(gcf,[filePreamble synType1{kk} '__ISonShaft' ISstr1 '__quartile' num2str(qqType1{kk}) '_neuron' num2str(tr) '_maxWidth100.jpg']); close;
      end
    end
  end
end

