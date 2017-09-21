% fptr=fopen('results.txt','w');
% cost='2 H2O [cytoplasm] + 1 ATP [cytoplasm] + 1 GTP [cytoplasm]'; 
% cost_relase='2 H+ [cytoplasm] + 2 phosphate [cytoplasm] + 1 ADP [cytoplasm] + 1 GDP [cytoplasm] + protein [cytoplasm]';
% 
% c=3;
[n a ctRNA]=xlsread('tRNA_modification_position_1.xlsx','test_codon');
% met1='protein [cytoplasm]';
for i=1:numel(ctRNA(:,1))
%     newModel1=newModel;
%     newModel1.c=zeros(numel(newModel1.c),1);
    met=cell2mat(ctRNA(i,2));
    met0=cell2mat(ctRNA(i,1));
%     eq=sprintf('%s => protein [cytoplasm] + %s',met0,met);
%     
%     rxnID='test_protein2';
%     rxnNames=sprintf('threonylcarbamoyladenylate synthase');
%     newModel1=addYeastReaction(newModel1,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf('%s => ',met);
    
    rxnID=sprintf('relase_tRNA%d',i);
    rxnNames=sprintf('threonylcarbamoyladenylate synthase');
    newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,1000,0,{''});

%     sol=testMetabolite(newModel1,met);
%    fprintf(fptr,'%s %f\n',met,sol.f);
% %     
%     sol1=testMetabolite(newModel1,met0);
%     fprintf(fptr,'%s %f\n',met0,sol1.f);
%     
%    sol2=testMetabolite(newModel1,'protein [cytoplasm]');
%    sol2=testMetabolite(newModel1,'YAL012W_peptide [cytoplasm]');
%    fprintf(fptr,'%s %f\n','protein [cytoplasm]',sol2.f);
% %     
  
end
% fclose(fptr);
%saveToExcel(newModel1,'model_test_protein_full_3.xlsx');
