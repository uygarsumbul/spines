function plotHeatMapAndSpecialNodes(tree, allMu, xyz)

  saturate = 0.01; disp(saturate)
miniX = 99999;
miniY = 99999;
hiCol = [1 0 0];
loCol = [0 1 1];
maSi  = 10;
punCo = [1 1 0];
allMeans = cell2mat(allMu); maxi = quantile(allMeans, 1-saturate); mini = quantile(allMeans, saturate);
disp([mini maxi])
figure;hold;
for kk = 2:numel(tree)
  for mm = 1:size(tree{kk}{4}{1},1)
    tt = tree{kk}{4}{1}(mm, :);
    miniX = min(miniX, min(tt(1)));
    miniY = min(miniY, min(tt(2)));
    ratio = min(1, max(0, (allMu{kk}(mm)-mini)/(maxi-mini)));
    plot3(tt(:,1),tt(:,2),tt(:,3),'Color',hiCol*ratio+loCol*(1-ratio),'MarkerSize',maSi,'Marker','.');
  end
end
for kk = 1:size(xyz, 1)
  plot3(xyz(kk,1), xyz(kk,2), xyz(kk,3), 'Color', punCo,'MarkerSize',18,'Marker','.');
end
text(miniX-100, miniY-40,num2str(round(100*mini)/100),'FontSize',16); text(miniX-20, miniY-40,num2str(round(100*maxi)/100),'FontSize',16);
for kk=0:99; plot3([miniX-100+kk miniX-100+kk+1], [miniY-20 miniY-20], [0 0],'LineWidth',10,'Color', hiCol*kk/99+loCol*(1-kk/99)); end;
view(0, 90); set(gcf,'Color','w','Position',[100 100 800 800]); axis off; axis square;
