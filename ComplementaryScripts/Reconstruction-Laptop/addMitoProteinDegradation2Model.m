function newModel= addMitoProteinDegradation2Model(newModel)

% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start Extracting\n');

fptr=fopen('utr_seq.txt','w');
fptr1=fopen('reactions.txt','w');

[a b geneData]=xlsread('TableS1.xlsx','Annotation');


for i=1:numel(geneData(:,1))
    
   gene  = cell2mat(geneData(i,2));
   if (strcmp(cell2mat(geneData(i,4)),'chrM')==1)
       
       deg_seq=cell2mat(geneData(i,8));
       
       
       [eq_product number_of_aa]= degradeMito(deg_seq);
       
       %termination
       eq=sprintf('%s_subunit [mitochondrion] => %s',gene,eq_product);
       
       rxnID=sprintf('%s_subunit_degradation',gene);
       rxnID=strrep(rxnID,'-','');
       rxnName=sprintf('degradation of %s Peptide',gene);
       
       newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'MitoDegradome'});
   end
 end 
     

fclose(fptr1);
fprintf('finishing\n');

