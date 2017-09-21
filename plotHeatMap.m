function plotHeatMap(tree, allMu)

allMeans = cell2mat(allMu); maxi = quantile(allMeans, 0.98); mini = quantile(allMeans, 0.02);
disp([mini maxi])

  figure;hold;
  for kk = 2:numel(tree)
    for mm = 1:size(tree{kk}{4}{1},1)
      tt = tree{kk}{4}{1}(mm, :);
      ratio = min(1, max(0, (allMu{kk}(mm)-mini)/(maxi-mini)));
      plot3(tt(:,1),tt(:,2),tt(:,3),'Color',[1 0 0]*ratio+[0 0 1]*(1-ratio),'MarkerSize',16,'Marker','.');
    end
  end
