function newModel= addProtein2ModelFromSequence(newModel,geneInfo)

 utr_length=geneInfo.UTR5;
 gene=geneInfo.name;
 %intial
 step_average = randomScanning(utr_length,.9,0.002);
 nATP_intial=1+step_average;
 nGTP_intial=2;
 %elongation
 elnogation_seq=geneInfo.seq(4:numel(geneInfo.seq)-3);
 [eq_substrates eq_product nATP_elnogation nGTP_elnogation] = elnogation(elnogation_seq);       
 
 %termination
 nATP_termination=1;
 nGTP_termination=2;
 total_atp=nATP_intial + nATP_elnogation + nATP_termination;
 total_gtp= nGTP_intial + nGTP_elnogation + nGTP_termination;
 total_h2o = total_atp + total_gtp;
 total_pi = total_atp + total_gtp;
 total_h = total_atp + total_gtp;
eq=sprintf('%d H2O [cytoplasm] + %d ATP [cytoplasm] + %d GTP [cytoplasm] + Met1-tRNA(AUG1) [cytoplasm] + %s => %s_peptide [cytoplasm] + %d ADP [cytoplasm] + %d GDP [cytoplasm] + %d phosphate [cytoplasm] + %d H+ [cytoplasm] + tRNA(AUG1) [cytoplasm] + %s', total_h2o, total_atp,total_gtp,eq_substrates,gene,total_atp,total_gtp,total_pi,total_h,eq_product);

 rxnID=sprintf('%s_translation',gene);
 rxnID=strrep(rxnID,'-','');
 rxnName=sprintf('biosynthesis of %s Peptide',gene);

 newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'Ribosome'});
 
 
fprintf('finishing\n');

