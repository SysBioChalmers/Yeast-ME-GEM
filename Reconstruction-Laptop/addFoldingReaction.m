function model=addFoldingReaction(model,n,m)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
p=proteins(n:m,2);
for i=1:numel(p)
 eq= sprintf('%s_peptide [cytoplasm] => %s_folding [cytoplasm]',cell2mat(p(i)),cell2mat(p(i)));
 rxnID=sprintf('%s_folding',cell2mat(p(i)));
 rxnID=strrep(rxnID,'-','');
 rxnName=sprintf('folding of the protein %s',cell2mat(p(i)));
 model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
%  
%  eq= sprintf('%s_folding [cytoplasm] => %s_subunit [cytoplasm]',cell2mat(p(i)),cell2mat(p(i)));
%  rxnID=sprintf('%s_subunit',cell2mat(p(i)));
%  rxnID=strrep(rxnID,'-','');
%  rxnName=sprintf('folding of the protein %s',cell2mat(p(i)));
%  model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
 
 
end
