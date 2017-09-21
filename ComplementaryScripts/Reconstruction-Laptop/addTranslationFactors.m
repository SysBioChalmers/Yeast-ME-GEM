function model=addTranslationFactors(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
n=find(ismember(proteins(:,1),'eIF1'));
m=find(ismember(proteins(:,1),'HBS1'));

I=n:m;
model= addProtein2Model_test(model,n,m);
model= addProteinDegradation2Model(model,n,m);

p=proteins(I,1:4);
k=1;

for i=1:numel(p(:,2))
           
    %transport the subunit to cytoplasm for degradation
    r=sprintf('%s_peptide [cytoplasm] => %s_folding [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_folding_cytoplasm',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s folding',cell2mat(p(i,2)));
    
    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
%misfolding
    r=sprintf('%s_folding [cytoplasm] => %s_misfolding [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_misfolding_cytoplasm',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
%refolding
    r=sprintf('%s_misfolding [cytoplasm] => %s_folding [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_refolding_cytoplasm',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s refolding',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});

    %misfolding degradation
    r=sprintf('%s_misfolding [cytoplasm] => %s_subunit [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_degradation_misfolding_cytoplasm',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding degradation',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    %misfolding dilution
    r=sprintf('%s_misfolding [cytoplasm] => ',cell2mat(p(i,2)));
    rxnID=sprintf('%s_dilution_misfolding_cytoplasm',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s misfolding dilution',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});

    
    %complex formation
    
    r=sprintf('%s_folding [cytoplasm] => %s_complex_formation [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_complex_formation',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s complex formation',cell2mat(p(i,3)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    % add digradtion
    eq = sprintf('%s_complex_formation [cytoplasm] => %s_subunit [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_degradation',cell2mat(p(i,2)));
    rxnName=sprintf('Degradation of %s',cell2mat(p(i,2)));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
  
    
    
    % add dilution
    eq = sprintf('%s_complex_formation [cytoplasm] => ',cell2mat(p(i,2)));
    rxnID=sprintf('%s_dilution',cell2mat(p(i,2)));
    rxnName=sprintf('dilution of %s',cell2mat(p(i,2)));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
          
end

% add complex
[a b rxn]=xlsread('TableS1.xlsx','Translation_factor_complex');

for i=1:numel(rxn(:,1))
    eq = cell2mat(rxn(i,3));
    rxnID= cell2mat(rxn(i,1));
    rxnName= cell2mat(rxn(i,2));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
end
