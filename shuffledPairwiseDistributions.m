cd ~/spines/repo
synType                                  = {{'spine', [-1 0 1], 0.00001}, {'spine', [-1 0 1], 0.8}, {'IS', 0, 0.00001}, {'IS', 0, 0.8}, {'IS', 1, 0.00001}, {'IS', 1, 0.8}, {'IS', [0 1], 0.00001}, {'IS', [0 1], 0.8}};
denType                                  = {'basal', 'apical'};
brType                                   = {'primaryOrIntermediateOrTerminal'};
preamble                                 = '/home/uygar/spines/results/branchLevelAnalysis/pairwiseDistancePlots/pairwiseDistance_';
load /home/uygar/spines/data/shuffledTrees/trueTrees.mat
for kk1=1:numel(synType); for kk2=1:numel(denType); for kk3=1:numel(brType); [~, ~, ~, pASD{kk1}{kk2}{kk3}] = extractPairwiseDistances(trees, synType{kk1}{1}, synType{kk1}{2}, synType{kk1}{3}, denType{kk2}, brType{kk3}); end; end; end;
for shuffle=1:100
  load(['/home/uygar/spines/data/shuffledTrees/shuffledTrees_' num2str(shuffle) '.mat']);
  for kk1=1:numel(synType); for kk2=1:numel(denType); for kk3=1:numel(brType); [~, ~, ~, pASD_shuffle{shuffle}{kk1}{kk2}{kk3}] = extractPairwiseDistances(trees, synType{kk1}{1}, synType{kk1}{2}, synType{kk1}{3}, denType{kk2}, brType{kk3}); end; end; end;
end

%for tr=1:numel(trees); try; [f,x] = ecdf(pASD{1}{1}{1}{tr}(:,1)); subplot(3,4,tr);plot(x,f);hold;for kk=1:100; [f,x] = ecdf(pASD_shuffle{kk}{1}{1}{1}{tr}(:,1)); plot(x,f,'Color','r'); end;xlim([0 5]); end; end;


%synType                                  = {{'spine', [-1 0 1], 0.00001}, {'spine', [-1 0 1], 0.8}};
%denType                                  = {'basal', 'apical'};
%brType                                   = {'primaryOrIntermediateOrTerminal'};
%load /home/uygar/spines/data/shuffledTreesValentine/trueTreesValentine.mat
%for kk1=1:numel(synType); for kk2=1:numel(denType); for kk3=1:numel(brType); [~, ~, ~, pASD{kk1}{kk2}{kk3}] = extractPairwiseDistances(trees, synType{kk1}{1}, synType{kk1}{2}, synType{kk1}{3}, denType{kk2}, brType{kk3}); end; end; end;
%for shuffle=1:100
%  load(['/home/uygar/spines/data/shuffledTreesValentine/shuffledTreesValentine_' num2str(shuffle) '.mat']);
%  for kk1=1:numel(synType); for kk2=1:numel(denType); for kk3=1:numel(brType); [~, ~, ~, pASD_shuffle{shuffle}{kk1}{kk2}{kk3}] = extractPairwiseDistances(trees, synType{kk1}{1}, synType{kk1}{2}, synType{kk1}{3}, denType{kk2}, brType{kk3}); end;
%end
%end
%end

for tr=1:numel(trees);
try
[f,x] = ecdf(pASD{1}{2}{1}{tr}(:,1)); subplot(3,4,tr);plot(x,f,'Color','k');hold;for kk=1:100; [f,x] = ecdf(pASD_shuffle{kk}{1}{2}{1}{tr}(:,1)); plot(x,f,'Color',[0.6 0.6 0.6]); end;xlim([0 0.2]); ylim([0 0.5]);
set(gca,'FontSize',16);
end
end
subplot(345); ylabel('empirical CDF','FontSize',16);subplot(3,4,10); xlabel('nearest neighbor distance (\mu m)','FontSize',16);set(gcf, 'Position',[100 100 900 550]); %saveas(gcf, 'allApicalSpines_zoomed.jpg');
