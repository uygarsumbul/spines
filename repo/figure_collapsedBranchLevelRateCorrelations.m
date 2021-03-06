load ~/misc/columbiaOnHabanero/spines/data/shuffledTrees/trueTrees12.mat
cd ~/Desktop/spines_27april2018/repo
minBranchLength    = 5;
cc                 = {[0 0 0],[32 188 16]/255,[180 0 180]/255};
sy                 = {'o','o','o'};
preamble           = '~/misc/columbiaOnHabanero/spines/repo/figures/april27_may2_branchRatePlot_';
preamble2          = '~/misc/columbiaOnHabanero/spines/repo/figures/april27_may2_branchRatePlot_noZeros_';
fontsize           = 16;
fontname           = 'Arial';
linwid             = 1;
per                = 0.02; %0; %0.05;
trendlim           = 0.02; %0.15;
parSet{1}     = {{'primary', 'intermediate', 'terminal'}, 'spine', 'IS'   , [-1 0 1], [0 1]   , 'basal'        , 1e-6, 1e-6};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'spine', 'IS'   , [-1 0 1], [0 1]   , 'apical'       , 1e-6, 1e-6};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'spine', 'IS'   , [-1 0 1], [0 1]   , 'basal'        , 0.8 , 0.8};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'spine', 'IS'   , [-1 0 1], [0 1]   , 'apical'       , 0.8 , 0.8};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'spine', 'spine', [-1 0 1], [-1 0 1], 'basal'        , 1e-6, 0.8};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'spine', 'spine', [-1 0 1], [-1 0 1], 'apical'       , 1e-6, 0.8};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'IS'   , 'IS'   , [0 1]   , [0 1]   , 'basal'        , 1e-6, 0.8};
parSet{end+1} = {{'primary', 'intermediate', 'terminal'}, 'IS'   , 'IS'   , [0 1]   , [0 1]   , 'apical'       , 1e-6, 0.8};

for pS = 1:numel(parSet)
  brType = parSet{pS}{1}; denType = parSet{pS}{6}; synType1 = parSet{pS}{2}; ISset1 = parSet{pS}{4}; qq1 = parSet{pS}{7}; tmp1=num2str(ISset1); tmp1=tmp1(tmp1~=' ');
                                                   synType2 = parSet{pS}{3}; ISset2 = parSet{pS}{5}; qq2 = parSet{pS}{8}; tmp2=num2str(ISset2); tmp2=tmp2(tmp2~=' ');

  if ~exist([preamble2 '-' synType1 '-' synType2 '-onShaft' tmp1 '-onShaft' tmp2 '-q' num2str(qq1) '-q' num2str(qq2) '-' denType '-allNeurons.pdf'])


  for kk4 = 1:numel(brType)
    [~,XYZsp{kk4},~, ~,dTsp{kk4}] = extractNodesOfInterest(trees, synType1, ISset1, qq1, denType, brType{kk4});
    [~,XYZis{kk4},~, ~,dTis{kk4}] = extractNodesOfInterest(trees, synType2, ISset2, qq2, denType, brType{kk4});
  end
  h1 = figure(1);hold; h3 = figure(3);hold; allTT1 = []; allTT1int = []; allTT1ter = []; allTT3 = []; allTT3int = []; allTT3ter = [];
  for tr = 1:numel(trees)
    if ~(tr==1 | tr==3)
      h2 = figure(2);hold; allTT2 = [];
      for kk4=1:numel(brType)
        for mm=2:numel(trees{tr})
          if dTsp{kk4}{tr}{mm} & abs(diff(trees{tr}{mm}{4}{5}([1 end],4)))>minBranchLength
            tt=[size(XYZsp{kk4}{tr}{mm},1); size(XYZis{kk4}{tr}{mm},1)]/sum(trees{tr}{mm}{4}{2});
            allTT1 = [allTT1 tt]; allTT2 = [allTT2 tt]; if kk4==2; allTT1int = [allTT1int tt]; end; if kk4==3; allTT1ter = [allTT1ter tt]; end;
            figure(1); plot(tt(1)',tt(2)','Marker',sy{kk4},'LineStyle','none','MarkerSize',8,'MarkerFaceColor',cc{kk4},'MarkerEdgeColor','k', 'LineWidth',linwid);
            figure(2); plot(tt(1)',tt(2)','Marker',sy{kk4},'LineStyle','none','MarkerSize',8,'MarkerFaceColor',cc{kk4},'MarkerEdgeColor','k', 'LineWidth',linwid);
            if ~any(tt==0)
              allTT3 = [allTT3 tt]; if kk4==2; allTT3int = [allTT3int tt]; end; if kk4==3; allTT3ter = [allTT3ter tt]; end;
              figure(3); plot(tt(1)',tt(2)','Marker',sy{kk4},'LineStyle','none','MarkerSize',8,'MarkerFaceColor',cc{kk4},'MarkerEdgeColor','k', 'LineWidth',linwid);
            end
          end
        end
      end
      figure(2); set(gca,'FontSize',fontsize, 'FontName',fontname, 'LineWidth', 2, 'box', 'off'); xlabel([synType1 ' rate (\mu m^{-1})'],'FontSize',fontsize); ylabel([synType2 ' rate (\mu m^{-1})'],'FontSize',fontsize); set(gcf, 'Color', 'w');
      xlim([quantile(allTT2(1, :), per) quantile(allTT2(1, :), 1-per)+eps]); ylim([quantile(allTT2(2, :), per) quantile(allTT2(2, :), 1-per)+eps]);
      saveas(h2, [preamble '-' synType1 '-' synType2 '-onShaft' tmp1 '-onShaft' tmp2 '-q' num2str(qq1) '-q' num2str(qq2) '-' denType '-neuron' num2str(tr) '.pdf']); close(h2);
    end
  end

  figure(1); % ALL
  lims1 = [quantile(allTT1(1, :), per) quantile(allTT1(1, :), 1-per)+eps quantile(allTT1(2, :), per) quantile(allTT1(2, :), 1-per)+eps];
  try
    coeffs = polyfit(allTT1int(1,:), allTT1int(2,:), 1); fittedX = linspace(lims1(1), lims1(2), 200); fittedY = polyval(coeffs, fittedX); plot(fittedX, fittedY, 'Color', cc{2}, 'LineStyle', '--', 'LineWidth', 4);
    [rrr, ppp, ~, ~] = corrcoef(allTT1int(1,:), allTT1int(2,:)); disp(['r_int=' num2str(rrr(1,2)) ', p_int=' num2str(ppp(1,2))]);
  end
  try
    coeffs = polyfit(allTT1ter(1,:), allTT1ter(2,:), 1); fittedX = linspace(lims1(1), lims1(2), 200); fittedY = polyval(coeffs, fittedX); plot(fittedX, fittedY, 'Color', cc{3}, 'LineStyle', '--', 'LineWidth', 4);
    [rrr, ppp, ~, ~] = corrcoef(allTT1ter(1,:), allTT1ter(2,:)); disp(['r_ter=' num2str(rrr(1,2)) ', p_ter=' num2str(ppp(1,2))]);
  end
  figure(1); set(gca,'FontSize',fontsize, 'FontName',fontname, 'LineWidth', 2, 'box', 'off'); xlabel([synType1 ' rate (\mu m^{-1})'],'FontSize',fontsize); ylabel([synType2 ' rate (\mu m^{-1})'],'FontSize',fontsize); set(gcf, 'Color', 'w');
  axis(lims1);
  saveas(h1, [preamble '-' synType1 '-' synType2 '-onShaft' tmp1 '-onShaft' tmp2 '-q' num2str(qq1) '-q' num2str(qq2) '-' denType '-allNeurons.pdf']); close(h1);

  figure(3); % NO ZEROS
  lims3 = [quantile(allTT3(1, :), per) quantile(allTT3(1, :), 1-per)+eps quantile(allTT3(2, :), per) quantile(allTT3(2, :), 1-per)+eps];
  try
    coeffs = polyfit(allTT3int(1,:), allTT3int(2,:), 1); fittedX = linspace(lims3(1), lims3(2), 200); fittedY = polyval(coeffs, fittedX); plot(fittedX, fittedY, 'Color', cc{2}, 'LineStyle', '--', 'LineWidth', 4);
    [rrr, ppp, ~, ~] = corrcoef(allTT3int(1,:), allTT3int(2,:)); disp(['r_int_noZeros=' num2str(rrr(1,2)) ', p_int_noZeros=' num2str(ppp(1,2))]);
  end
  try
    coeffs = polyfit(allTT3ter(1,:), allTT3ter(2,:), 1); fittedX = linspace(lims3(1), lims3(2), 200); fittedY = polyval(coeffs, fittedX); plot(fittedX, fittedY, 'Color', cc{3}, 'LineStyle', '--', 'LineWidth', 4);
    [rrr, ppp, ~, ~] = corrcoef(allTT3ter(1,:), allTT3ter(2,:)); disp(['r_ter_noZeros=' num2str(rrr(1,2)) ', p_ter_noZeros=' num2str(ppp(1,2))]);
  end
  figure(3); set(gca,'FontSize',fontsize, 'FontName',fontname, 'LineWidth', 2, 'box', 'off'); xlabel([synType1 ' rate (\mu m^{-1})'],'FontSize',fontsize); ylabel([synType2 ' rate (\mu m^{-1})'],'FontSize',fontsize); set(gcf, 'Color', 'w');
  axis(lims3);
  saveas(h3, [preamble2 '-' synType1 '-' synType2 '-onShaft' tmp1 '-onShaft' tmp2 '-q' num2str(qq1) '-q' num2str(qq2) '-' denType '-allNeurons.pdf']); close(h3);

end

end
