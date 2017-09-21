function [flag, lastNode] = primaryOrIntermediateOrTerminal_function(denTyp, parentDenTyp, hasChildren, validDendriteTypes, primaryOrIntermediateOrTerminal)
  % denTyp       = trees{tr}{br}{4}{5}(1:end-1,19-7);
  % parentDenTyp = trees{tr}{thisParent}{4}{5}(1:end-1,19-7);
  % hasChildren  = ~isempty(trees{tr}{br}{2});
  flag                      = false;
  membership                = ismember(denTyp, validDendriteTypes);
  memberCount               = nnz(membership);
  somaticCount              = nnz(denTyp==0);
  remainingCount            = numel(denTyp) - somaticCount - memberCount;
  lastNode                  = numel(denTyp);
  switch primaryOrIntermediateOrTerminal
    case 'primary'
      if hasChildren & (any(denTyp==0) || all(parentDenTyp==0)) & any(membership)
	lastNode            = find(membership, 1);
	flag                = true;
      end
    case 'intermediate'
      if any(parentDenTyp>0)
	flag                = hasChildren & (memberCount>remainingCount);
      end
    case 'terminal'
      flag                  = ~hasChildren & (memberCount>remainingCount);
      if any(parentDenTyp==0)
	lastNode            = find(membership, 1, 'last');
      end
    case 'intermediateOrTerminal'
      if any(parentDenTyp>0)
	flag                = (memberCount>remainingCount);
      else
	flag                = ~hasChildren & (memberCount>remainingCount);
	if any(parentDenTyp==0)
	  lastNode          = find(membership, 1, 'last');
	end
      end
    case 'primaryOrIntermediateOrTerminal'
      flag                  = (memberCount>remainingCount);
      if any(denTyp==0)
	lastNode            = find(membership, 1, 'last');
      end
    otherwise
      warning('Unexpected branch type');
  end
  
  
