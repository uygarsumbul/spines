function save1000ShuffledTrees

code{1}{1}                               = 'aibs_8_spinerecon_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_9_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_13_mouse2_spinerecon_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_14_stitched_connected_corrected_labeled.eswc'; % added on May 19, 2017
code{end+1}{1}                           = 'aibs_18_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_19_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_20_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_21_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_22_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_23_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_24_stitched_connected_labeled.eswc';
code{end+1}{1}                           = 'aibs_25_stitched_connected_labeled.eswc';
directory                                = '/home/uygar/spines/data/upsampledFinalReconstructions/'; % _updated3-22-17/';
options.resolution                       = [0.12 0.12 0.1];
options.resolution20x                    = [0.62 0.62 0.75];

parfor kk = 1:1000
  trees                                  = readDataset_spines_mod_shuffleNodes(code,directory,options);
  mySaveFunction(kk, trees);  
end
trees                                    = readDataset_spines_mod(code,directory,options);
save /home/uygar/spines/data/shuffledTrees/trueTrees.mat trees

function mySaveFunction(kk, trees)
filename = ['/home/uygar/spines/data/shuffledTrees/shuffledTrees_' num2str(kk) '.mat'];
save(filename, 'trees');
