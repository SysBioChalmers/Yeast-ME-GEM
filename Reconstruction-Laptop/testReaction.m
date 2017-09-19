function sol=testReaction(testModel,rxn)
testModel.c=zeros(numel(testModel.c),1);
testModel.c(find(ismember(testModel.rxns,rxn)))=1;
sol=solveLP(testModel);