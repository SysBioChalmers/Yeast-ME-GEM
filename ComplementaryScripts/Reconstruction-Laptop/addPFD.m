function model=addPFD(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
I=find(ismember(proteins(:,1),'PFD'));
I1=find(ismember(proteins(:,1),'PFD'));
n=min(I);
m=max(I);
model= addProtein2Model_test(model,n,m);
model= addProteinDegradation2Model(model,n,m);
I=n:m;
p=proteins(I,1:4);
k=1;

for i=1:numel(p(:,2))
           
    %transport the subunit to cytoplasm for degradation
    r=sprintf('%s_peptide [cytoplasm] => %s_folding [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_folding_cytoplasm',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s folding',cell2mat(p(i,2)));
    
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
     
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

    
          
end

subunit=proteins(I,2);

eq=sprintf('%s_folding [cytoplasm]', cell2mat(subunit(1)));
eq_deg=sprintf('%s_subunit [cytoplasm]', cell2mat(subunit(1)));

for i=2:numel(subunit)
   eq=sprintf('%s + %s_folding [cytoplasm]',eq,cell2mat(subunit(i)));
   eq_deg=sprintf('%s + %s_subunit [cytoplasm]',eq_deg,cell2mat(subunit(i)));
end

eq = sprintf('%s => PFD [cytoplasm]',eq);
rxnID=sprintf('PFD_biosynthesis');

rxnName=sprintf('biosynthesis of PFD');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

% add digradtion
eq = sprintf('PFD [cytoplasm] => %s',eq_deg);
rxnID=sprintf('PFD_degradation');

rxnName=sprintf('Degradation of PFD');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

% add dilution
eq = sprintf('PFD [cytoplasm] => ');
rxnID=sprintf('PFD_dilution');

rxnName=sprintf('Dilution of PFD');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
