function model=addHsp60(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
I=find(ismember(proteins(:,1),'HSP60'));
I1=find(ismember(proteins(:,1),'HSP60'));
n=min(I);
m=max(I);
model= addProtein2Model_test(model,n,m); 
model= addProteinDegradation2Model(model,n,m);
I=n:m;
p=proteins(I,1:4);
k=1;

for i=1:numel(p(:,2))
           
    %transport the subunit to mitochondrion for degradation
    r=sprintf('%s_peptide [cytoplasm] => %s_peptide [mitochondrial membrane]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_importing_%s_IMS_mitochondrion',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s transporting',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    %
    %transport the subunit to mitochondrion for degradation
    r=sprintf('%s_peptide [mitochondrial membrane] => %s_peptide [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_importing_%s_Matrix',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s importing into Matrix',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    
    %
    %transport the subunit to mitochondrion for degradation
    r=sprintf('%s_peptide [mitochondrion] => %s_folding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_folding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s folding',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    

    %misfolding
    r=sprintf('%s_folding [mitochondrion] => %s_misfolding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_misfolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    
    %refolding
    r=sprintf('%s_misfolding [mitochondrion] => %s_folding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_refolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s refolding',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    
    %misfolding degradation
    r=sprintf('%s_misfolding [mitochondrion] => %s_subunit [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_degradation_misfolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding degradation',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    r=sprintf('%s_subunit [mitochondrion] => %s_subunit [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_subunit_export_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s subunit exporting',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    
    
    %misfolding dilution
    r=sprintf('%s_misfolding [mitochondrion] => ',cell2mat(p(i,2)));
    rxnID=sprintf('%s_dilution_misfolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s misfolding dilution',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});

    
          
end

subunit=proteins(I,2);

eq=sprintf('%s_folding [mitochondrion]', cell2mat(subunit(1)));
eq_deg=sprintf('%s_subunit [mitochondrion]', cell2mat(subunit(1)));


eq = sprintf('14 %s + 14 YOR020C_folding [mitochondrion] => mtHSP10_mtHSP60 [mitochondrion]',eq);
rxnID=sprintf('mtHSP10_mtHSP60_biosynthesis');

rxnName=sprintf('biosynthesis of mtHSP10_mtHSP60');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

% add digradtion
eq = sprintf('mtHSP10_mtHSP60 [mitochondrion] => 14 %s + 14 YOR020C_subunit [mitochondrion]',eq_deg);
rxnID=sprintf('mtHSP10_mtHSP60_degradation');
rxnName=sprintf('Degradation of mtHSP10_mtHSP60');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

% add dilution
eq = sprintf('mtHSP10_mtHSP60 [mitochondrion] => ');
rxnID=sprintf('mtHSP10_mtHSP60_dilution');
rxnName=sprintf('Dilution of mtHSP10_mtHSP60');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
