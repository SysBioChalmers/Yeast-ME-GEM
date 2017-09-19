function sol=testMetabolite(model,met)
model.c=zeros(numel(model.c),1);
eq=sprintf('%s => ',met);
rxnID='test_biomass';
rxnNames=sprintf('test biomass');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,1,{''});
sol=solveLP(model);
%saveToExcel(model,'model_test_protein_full_4.xlsx');