cd ~/spines/repo
code{1}{1}                               = 'aibs_9_stitched_connected_labeled.eswc';
code{2}{1}                               = 'aibs_18_stitched_connected_labeled.eswc';
code{3}{1}                               = 'aibs_19_stitched_connected_labeled.eswc';
code{4}{1}                               = 'aibs_20_stitched_connected_labeled.eswc';
code{5}{1}                               = 'aibs_21_stitched_connected_labeled.eswc';
code{6}{1}                               = 'aibs_22_stitched_connected_labeled.eswc';
code{7}{1}                               = 'aibs_23_stitched_connected_labeled.eswc';
code{8}{1}                               = 'aibs_24_stitched_connected_labeled.eswc';
code{9}{1}                               = 'aibs_25_stitched_connected_labeled.eswc';
directory                                = '/home/uygar/spines/data/upsampledFinalReconstructions/update_3mar2017/';
options.resolution                       = [0.12 0.12 0.1];
options.resolution20x                    = [0.62 0.62 0.75];
filePreamble                             = '/home/uygar/spines/results/batchHeatMaps';
synType                                  = {'spine', 'IS'};
IStype                                   = {[-1 0 1], 0, 1, [0 1]};
qqType                                   = {0.00001, 0.8};
denType                                  = {'basalAndApical'}; % {'basal', 'apical', 'basalAndApical'};
brType                                   = {'primal', 'intermediate', 'terminal', 'intermediateOrTerminal', 'primalOrIntermediateOrTerminal'};


maxiMinusMini=zeros(2, 4, 2, 100, 9);

for iter=1:100
trees                                    = readDataset_spines_mod_shuffleNodesWithinBranchTypes(code,directory,options);

for tr = 1:numel(trees)
  trLen(tr) = treeLength(trees{tr},1,true);
end
for kk1=1:numel(synType)
  for kk2=1:numel(IStype)
    for kk3=1:numel(qqType)
      for kk4=1:numel(denType)
          [~, ~, ~, pairAndSomaDistances] = pairwiseDistanceDistribution_scratch2(trees, synType{kk1}, IStype{kk2}, qqType{kk3}, denType{kk4}, 'primalOrIntermediateOrTerminal'); %brType{kk5});
          for tr = 1:numel(trees)
if ~isempty(pairAndSomaDistances{tr})
              [allMu, mini, maxi] = arborHeatMap(trees{tr}, pairAndSomaDistances{tr}, 50, false);
              maxiMinusMini(kk1,kk2,kk3,iter,tr)=maxi-mini;
end
          end
      end
    end
  end
end

end
