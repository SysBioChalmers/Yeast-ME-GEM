function model=addTranslationFactorsDilution(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
I1=find(ismember(proteins(:,1),'eEF3'));
I2=find(ismember(proteins(:,1),'RLI1'));
n=min(I1);
m=max(I2);
I=n:m;
p=proteins(I,1:4);

for i=1:numel(I)
    eq = sprintf('%s_folding [cytoplasm] => ',cell2mat(p(i,2)));
    rxnID=sprintf('%s_dilution',cell2mat(p(i,2)));
    
    rxnName=sprintf('Dilution of %s',cell2mat(p(i,2)));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
end
