function newModel= addProteinDegradation2ModelFromSequence(newModel,geneInfo)

gene=geneInfo.name;  
[eq_product number_of_aa]= degrade(geneInfo.seq);       
 
 %termination
nATP =floor(1.3*number_of_aa);
eq=sprintf('%d H2O [cytoplasm] + %d ATP [cytoplasm] + %s_subunit [cytoplasm] => %s + %d ADP [cytoplasm] + %d phosphate [cytoplasm] + %d H+ [cytoplasm]',nATP, nATP,gene,eq_product,nATP,nATP,nATP);

rxnID=sprintf('%s_subunit_degradation',gene);
 rxnID=strrep(rxnID,'-','');
 rxnName=sprintf('degradation of %s Peptide',gene);

 newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'Proteasome'});

 

fprintf('finishing\n');

