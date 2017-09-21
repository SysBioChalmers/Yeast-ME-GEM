function newModel= addProteinDilution2Model(newModel,n,m)

% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start dilution\n');



[a b geneData]=xlsread('TableS1.xlsx','Annotation');


for i=n:m
   gene  = cell2mat(geneData(i,2));
   Start = cell2mat(geneData(i,12));
   End   = cell2mat(geneData(i,13));
   direction = cell2mat(geneData(i,5));
   TSS=cell2mat(geneData(i,11));
   intron=cell2mat(geneData(i,10));
   if ~(strcmp(cell2mat(geneData(i,4)),'chrM'))
      chr= getChromosomeNumber(cell2mat(geneData(i,4)));

   
   eq=sprintf('%s_peptide [cytoplasm] => ',gene);

  rxnID=sprintf('%s_dilution_peptide',gene);
  rxnID=strrep(rxnID,'-','');
  rxnName=sprintf('dilution of %s Peptide',gene);

  newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{''});

 end % end if chrM
     
end

fprintf('End dilution\n');

