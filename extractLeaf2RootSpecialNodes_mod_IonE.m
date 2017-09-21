function [allLSPD, largeSpineXYZ] = extractLeaf2RootSpecialNodes_mod_IonE(allTrees, dendriteID, distMin, distMax, myWindow, resolution)

allLSPD                                = cell(0);
largeSpineXYZ                          = cell(0);
thisFeature = 18;
for tr = 1:numel(allTrees)

  leaf2RootPaths                       = extractLeaf2RootPaths(allTrees{tr});
  largeSpinePairDistances              = [];
  largeSpineXYZ{tr}{kk}                = [];
  for kk = 1:numel(leaf2RootPaths)
    largeSpines                        = [];
    largeSpineXYZ{tr}{kk}              = cell(1,numel(leaf2RootPaths{kk}));
    for mm = 1:numel(leaf2RootPaths{kk})
      thisSeg                          = leaf2RootPaths{kk}(mm);
      withinWindow                     = allTrees{tr}{thisSeg}{4}{5}(:,11-7)>distMin & allTrees{tr}{thisSeg}{4}{5}(:,11-7)<distMax & ismember(allTrees{tr}{thisSeg}{4}{5}(:,19-7), dendriteID);
      largeSpines                      = find(allTrees{tr}{thisSeg}{4}{5}(:,thisFeature-7)==0 & withinWindow);
      largeSpinePairDistances          = [largeSpinePairDistances; abs(diff(allTrees{tr}{thisSeg}{4}{5}(largeSpines,11-7)))];

      withinWindow                     = allTrees{tr}{thisSeg}{4}{5}(:,11-7)>=(distMin+distMax)/2-myWindow/2-resolution/2 & allTrees{tr}{thisSeg}{4}{5}(:,11-7)<=(distMin+distMax)/2+myWindow/2+resolution/2 & ismember(allTrees{tr}{thisSeg}{4}{5}(:,19-7), dendriteID);

      largeSpineXYZ{tr}{kk}{mm}        = allTrees{tr}{thisSeg}{4}{1}(withinWindow, :);

    end
  end
  allLSPD{tr} = largeSpinePairDistances;
end
