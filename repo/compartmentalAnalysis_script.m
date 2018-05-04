load ~/spines/data/shuffledTrees/trueTrees12.mat
spineOrIS={'spine', 'IS'};
synType            = {{'spine', [-1 0 1], 0, '-101'}, {'spine', [-1 0 1], 0.8, '-101'}, {'IS', 0, 0, '0'}, {'IS', 0, 0.8, '0'}, {'IS', 1, 0, '1'}, {'IS', 1, 0.8, '1'}, {'IS', [0 1], 0, '01'}, {'IS', [0 1], 0.8, '01'}};
denType            = {{'basal', {'primary', 'intermediate', 'terminal'}}, {'apical', {'primary', 'intermediate', 'terminal'}}, {'apicalTrunk', {'primaryOrIntermediateOrTerminal'}}, {'apicalOblique', {'intermediate', 'terminal'}}, {'apicalTuft', {'intermediate', 'terminal'}}};
for kk1=1:numel(synType)
  for kk2=1:numel(denType)
    for kk3=1:numel(denType{kk2}{2})
      [compartmentalLength eventCount] = compartments_lengthAndEventCount(trees, synType{kk1}{1},denType{kk2}{1},denType{kk2}{2}{kk3},synType{kk1}{3},synType{kk1}{2});
      for tr=1:numel(trees)
	  fn=['/home/uygar/spines/results/branchLevelAnalysis/compartments_' synType{kk1}{1} '_ISonShaft' synType{kk1}{4} '_quartile' num2str(synType{kk1}{3}) '_' denType{kk2}{1} '_' denType{kk2}{2}{kk3} '_neuron' num2str(tr)];
	  fn=[fn '_length' num2str(compartmentalLength(tr)) '_count' num2str(eventCount(tr)) '.txt'];
	  fid=fopen(fn,'w'); fwrite(fid,0); fclose(fid);
      end
    end
  end
end
