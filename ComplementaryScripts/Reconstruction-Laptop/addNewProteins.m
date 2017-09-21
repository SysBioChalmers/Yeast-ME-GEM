function model=addNewProteins(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
I=find(ismember(proteins(:,1),'New_protein'));
n=min(I);
m=max(I);
model= addProtein2Model_test(model,n,m);
model= addProteinDegradation2Model(model,n,m);

p=proteins(I,1:4);
k=1;

for i=1:numel(p(:,2))
    %transporting
    r=sprintf('%s_peptide [cytoplasm] => %s_peptide [mitochondrial membrane]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_importing_cytoplasm_IMS',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('Importing %s into IMS cytoplasm',cell2mat(p(i,2)));
    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

    r=sprintf('%s_peptide [mitochondrial membrane] => %s_peptide [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_importing_cytoplasm_Matrix',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('Importing %s into Matrix cytoplasm',cell2mat(p(i,2)));
    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    %folding
    r=sprintf('%s_peptide [mitochondrion] => %s_folding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_folding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s folding',cell2mat(p(i,3)));
    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
    %misfolding
      r=sprintf('%s [mitochondrion] => %s_misfolding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
      rxnID=sprintf('%s_misfolding_mitochondrion',cell2mat(p(i,2)));
      rxnID=strrep(rxnID,'-','');
      rxnName=sprintf('%s Misfolding',cell2mat(p(i,2)));
      model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});

      %refolding
      r=sprintf('%s_misfolding [mitochondrion] => %s [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
      rxnID=sprintf('%s_refolding_mitochondrion',cell2mat(p(i,2)));
      rxnID=strrep(rxnID,'-','');
      rxnName=sprintf('%s Misfolding',cell2mat(p(i,2)));
      model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});

      
    %misfolding degradation
    r=sprintf('%s_misfolding [mitochondrion] => %s_subunit [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_degradation_misfolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding degradation',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    %misfolding dilution
    r=sprintf('%s_misfolding [mitochondrion] => ',cell2mat(p(i,2)));
    rxnID=sprintf('%s_dilution_misfolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding dilution',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    
    
    %complex formation
    r=sprintf('%s_folding [mitochondrion] => %s [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_complex_formation',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s import to mitochondrion',cell2mat(p(i,3)));
    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    
    % add digradtion
    eq = sprintf('%s [mitochondrion] => %s_subunit [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_degradation',cell2mat(p(i,2)));
    rxnName=sprintf('Degradation of %s',cell2mat(p(i,2)));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
    eq = sprintf('%s_subunit [mitochondrion] => %s_subunit [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_subunit_transport',cell2mat(p(i,2)));
    rxnName=sprintf('subunit export %s',cell2mat(p(i,2)));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
    
     % add dilution
    eq = sprintf('%s [mitochondrion] => ',cell2mat(p(i,2)));
    rxnID=sprintf('%s_dilution',cell2mat(p(i,2)));
    rxnName=sprintf('dilution of %s',cell2mat(p(i,2)));
    model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
    
          
end
