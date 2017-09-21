function newModel= addProteinDilution2ModelFromSequence(newModel,geneInfo)
gene  = geneInfo.name;

eq=sprintf('%s_folding [cytoplasm] => ',gene);
rxnID=sprintf('%s_dilution',gene);
rxnID=strrep(rxnID,'-','');
rxnName=sprintf('dilution of %s Peptide',gene);




