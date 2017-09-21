function newModel= addMitoProtein2Model(newModel)

% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start Extracting\n');

fptr=fopen('utr_seq.txt','w');
fptr1=fopen('reactions.txt','w');

[a b geneData]=xlsread('TableS1.xlsx','Annotation');
n=numel(geneData(:,7));

for i=2:1294
    gene  = cell2mat(geneData(i,2));
    Start = cell2mat(geneData(i,12));
    End   = cell2mat(geneData(i,13));
    direction = cell2mat(geneData(i,5));
    TSS=cell2mat(geneData(i,11));
    intron=cell2mat(geneData(i,10));
    if (strcmp(cell2mat(geneData(i,4)),'chrM'))
        %intial
        seq=cell2mat(geneData(i,8));
        
        
        nATP_intial=0;
        nGTP_intial=1;
        %elongation
        
        %seq=cell2mat(geneData(i,8));
        elnogation_seq=seq(4:numel(seq)-3);
        [eq_substrates eq_product nATP_elnogation nGTP_elnogation] = elnogationMito(elnogation_seq);
        
        %termination
        nATP_termination=0;
        nGTP_termination=1;
        total_atp=nATP_intial + nATP_elnogation + nATP_termination;
        total_gtp= nGTP_intial + nGTP_elnogation + nGTP_termination;
        total_h2o = total_atp + total_gtp;
        total_pi = total_atp + total_gtp;
        total_h = total_atp + total_gtp;
        eq=sprintf('%d H2O [mitochondrion] + %d ATP [mitochondrion] + %d GTP [mitochondrion] + Met1-tRNA(AUG1) [mitochondrion] + %s => %s_peptide [mitochondrion] + %d ADP [mitochondrion] + %d GDP [mitochondrion] + %d phosphate [mitochondrion] + %d H+ [mitochondrion] + tRNA(AUG1) [mitochondrion] + %s', total_h2o, total_atp,total_gtp,eq_substrates,gene,total_atp,total_gtp,total_pi,total_h,eq_product);
        %eq=sprintf('%d H2O [cytoplasm] + %d ATP [cytoplasm] + %d GTP [cytoplasm] + Met1-tRNA(AUG1) [cytoplasm] => %s_peptide [cytoplasm] + %d ADP [cytoplasm] + %d GDP [cytoplasm] + %d phosphate [cytoplasm] + %d H+ [cytoplasm] + %s', total_h2o, total_atp,total_gtp,gene,total_atp,total_gtp,total_pi,total_h,'tRNA(AUG1) [cytoplasm]');
%        eq=sprintf(' => %s_peptide [mitochondrion]',gene);
        rxnID=sprintf('%s_translation_mitochondrion',gene);
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('biosynthesis of %s Peptide into mitochondrion',gene,gene);
        
        newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'Ribosome'});
        
%         eq=sprintf('%s_peptide [mitochondrion] => %s_folding [mitochondrion]',gene,gene);
%         rxnID=sprintf('%s_folding_mitochondrion',gene);
%         rxnID=strrep(rxnID,'-','');
%         rxnName=sprintf('folding of %s into mitochondrion',gene);
%         
%         newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'Ribosome'});
%         
%         %msfolding
%         eq=sprintf('%s_peptide [mitochondrion] => %s_folding [mitochondrion]',gene,gene);
%         rxnID=sprintf('%s_folding_mitochondrion',gene);
%         rxnID=strrep(rxnID,'-','');
%         rxnName=sprintf('folding of %s into mitochondrion',gene);
%         
%         newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'Ribosome'});
        
        
        fprintf(fptr1,'%s\n',eq);
        fprintf(fptr1,'\tfolding of %s Peptide\t%s_peptide [cytoplasm] => %s_folding [cytoplasm]\n',gene,gene,gene);
 
        
    end
end
