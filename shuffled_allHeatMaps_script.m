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
trees                                    = readDataset_spines_mod_shuffleNodesWithinBranchTypes(code,directory,options);
filePreamble                             = '/home/uygar/spines/results/batchHeatMaps';
synType                                  = {'spine', 'IS'};
IStype                                   = {[-1 0 1], 0, 1, [0 1]};
qqType                                   = {0.00001, 0.8};
denType                                  = {'basalAndApical'}; % {'basal', 'apical', 'basalAndApical'};
brType                                   = {'primal', 'intermediate', 'terminal', 'intermediateOrTerminal', 'primalOrIntermediateOrTerminal'};
for tr = 1:numel(trees)
  trLen(tr) = treeLength(trees{tr},1,true);
end
for kk1=1:numel(synType)
  for kk2=1:numel(IStype)
    for kk3=1:numel(qqType)
      for kk4=1:numel(denType)
%        for kk5=1:numel(brType)
          [~, ~, ~, pairAndSomaDistances] = pairwiseDistanceDistribution_scratch2(trees, synType{kk1}, IStype{kk2}, qqType{kk3}, denType{kk4}, 'primalOrIntermediateOrTerminal'); %brType{kk5});
          for tr = 1:numel(trees)
            lambda = size(pairAndSomaDistances{tr},1)/trLen(tr);
            if lambda > 16/100 %70
              [allMu, mini, maxi] = arborHeatMap(trees{tr}, pairAndSomaDistances{tr}, 16/lambda, false);
              plotHeatMap(trees{tr}, allMu);
              switch kk2
                case 1
                  saveas(gcf,[filePreamble synType{kk1} '_ISonShaft-101_quartile' num2str(qqType{kk3}) num2str(denType{kk4}) '_neuron' num2str(tr) '_maxWidth100_shuffled.jpg']); %brType{kk5} '.jpg']);
                case 4
                  saveas(gcf,[filePreamble synType{kk1} '_ISonShaft01_quartile' num2str(qqType{kk3}) num2str(denType{kk4}) '_neuron' num2str(tr) '_maxWidth100_shuffled.jpg']); %brType{kk5} '.jpg']);
                otherwise
                  saveas(gcf,[filePreamble synType{kk1} '_ISonShaft' num2str(IStype{kk2}) '_quartile' num2str(qqType{kk3}) denType{kk4} '_neuron' num2str(tr) '_maxWidth100_shuffled.jpg']); %brType{kk5} '.jpg']);
              end
              close;
            end
          end
%        end
      end
    end
  end
end
