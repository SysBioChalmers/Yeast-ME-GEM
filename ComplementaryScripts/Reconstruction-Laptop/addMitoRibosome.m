function model=addMitoRibosome(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
I1=find(ismember(proteins(:,1),'Mitchochondrial Ribosomal Large Subunit'));
I2=find(ismember(proteins(:,1),'Mitchochondrial Ribosomal Small Subunit'));
I=union(I1,I2);
p=proteins(I,1:4);
k=1;
for i=1:numel(p(:,2))
    %transport the subunit to cytoplasm for degradation
    chr=cell2mat(p(i,4));
    if ~strcmp(chr,'chrM')
      
        
    r=sprintf('%s_peptide [cytoplasm] => %s_peptide [mitochondrial membrane]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_importing_mitochondrion_IMS',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('importing %s into mitochondrion IMS',cell2mat(p(i,2)));
    
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{'TOM-TIM'});
    
    
    r=sprintf('%s_peptide [mitochondrial membrane] => %s_peptide [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_importing_mitochondrion_Matrix',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('importing %s into mitochondrion',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{'TOM-TIM'});
    
    
    
    %complex
    r=sprintf('%s_peptide [mitochondrion] => %s [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_folding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('folding %s into mitochondrion',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{'TOM-TIM'});
    
    %misfolding
    r=sprintf('%s_peptide [mitochondrion] => %s_misfolding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_misfolding_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('%s Misfolding',cell2mat(p(i,2)));
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    %refolding
    r=sprintf('%s_misfolding [mitochondrion] => %s_peptide [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
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
    
    %misfolding dilution
     r=sprintf('%s_misfolding [mitochondrion] => ',cell2mat(p(i,2)));
     rxnID=sprintf('%s_dilution_misfolding_mitochondrion',cell2mat(p(i,2)));
     rxnID=strrep(rxnID,'-','');
     rxnName=sprintf('%s Misfolding dilution',cell2mat(p(i,2)));
     model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    

   %export
    r=sprintf('%s_subunit [mitochondrion] => %s_subunit [cytoplasm]',cell2mat(p(i,2)),cell2mat(p(i,2)));
    rxnID=sprintf('%s_exporting_from_mitochondrion',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('exporting %s from mitochondrion',cell2mat(p(i,2)));
    
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});
    else
        eq=sprintf('%s_peptide [mitochondrion] => %s_folding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
        rxnID=sprintf('%s_folding_mitochondrion',cell2mat(p(i,2)));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('folding of %s into mitochondrion',cell2mat(p(i,2)));
        
        model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
        
        %misfolding
        eq=sprintf('%s_folding [mitochondrion] => %s_misfolding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
        rxnID=sprintf('%s_misfolding_mitochondrion',cell2mat(p(i,2)));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('misfolding of %s into mitochondrion',cell2mat(p(i,2)));
        model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
        
        %refolding
        eq=sprintf('%s_misfolding [mitochondrion] => %s_folding [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
        rxnID=sprintf('%s_refolding_mitochondrion',cell2mat(p(i,2)));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('refolding of %s into mitochondrion',cell2mat(p(i,2)));
        model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
  
        %misfolding degradation
        eq=sprintf('%s_misfolding [mitochondrion] => %s_subunit [mitochondrion]',cell2mat(p(i,2)),cell2mat(p(i,2)));
        rxnID=sprintf('%s_degradation_misfolding_mitochondrion',cell2mat(p(i,2)));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('degradation misfolding of %s into mitochondrion',cell2mat(p(i,2)));
        model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
        
        %misfolding dilution
        eq=sprintf('%s_misfolding [mitochondrion] => ',cell2mat(p(i,2)));
        rxnID=sprintf('%s_dilution_misfolding_mitochondrion',cell2mat(p(i,2)));
        rxnID=strrep(rxnID,'-','');
        rxnName=sprintf('dilution misfolding of %s into mitochondrion',cell2mat(p(i,2)));
        model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
        
        
     end
    
end

subunit=proteins(I,2);

eq=sprintf('%s [mitochondrion]', cell2mat(subunit(1)));
eq_deg=sprintf('%s_subunit [mitochondrion]', cell2mat(subunit(1)));
for i=2:numel(subunit)
     chr=cell2mat(p(i,4));
    if strcmp(chr,'chrM')
         eq=sprintf('%s + %s_folding [mitochondrion]',eq,cell2mat(subunit(i)));
          eq_deg=sprintf('%s + %s_subunit [mitochondrion]',eq_deg,cell2mat(subunit(i)));
    else
        eq=sprintf('%s + %s [mitochondrion]',eq,cell2mat(subunit(i)));
    eq_deg=sprintf('%s + %s_subunit [mitochondrion]',eq_deg,cell2mat(subunit(i)));
    end
end

eq = sprintf('%s => Ribosome [mitochondrion]',eq);
rxnID=sprintf('MitoRibosome_biosynthesis');

rxnName=sprintf('biosynthesis of mitochondrial ribosome');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

% add digradtion
eq = sprintf('Ribosome [mitochondrion] => %s',eq_deg);
rxnID=sprintf('MitoRibosome_degradation');

rxnName=sprintf('Degradation of mitochondrial ribosome');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

% add digradtion
eq = sprintf('Ribosome [mitochondrion] => ');
rxnID=sprintf('MitoRibosome_dilution');

rxnName=sprintf('Dilution of mitochondrial ribosome');
model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
