function model=addRibosome(model)
[a b proteins]=xlsread('TableS1.xlsx','Annotation');
I=find(ismember(proteins(:,1),'Ribosome'));
p=proteins(I,1:3);
k=1;
for i=1:numel(p(:,2))
    aa=regexp(cell2mat(p(i,3)),'S\d*','match');
    if isempty(aa)
         aa=regexp(cell2mat(p(i,3)),'L\d*','match');
         if isempty(aa)
             aa=regexp(cell2mat(p(i,3)),'P\d*','match');
             
             if isempty(aa)
                 aa = p(i,3);
             end
         end
    end
    rbosomeSubunit(i) =aa;
         
    r=sprintf('%s_folding [cytoplasm] => %s [cytoplasm]',cell2mat(p(i,2)),cell2mat(rbosomeSubunit(i)));
    
    rxnID=sprintf('%s_biosynthesis',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('biosynthesis %s from %s',cell2mat(rbosomeSubunit(i)),cell2mat(p(i,2)));
    
    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
    
    %add degradation
    r=sprintf('%s_subunit [cytoplasm] => %s_subunit [cytoplasm]',cell2mat(rbosomeSubunit(i)),cell2mat(p(i,2)));
    
    rxnID=sprintf('%s_degradation',cell2mat(p(i,2)));
    rxnID=strrep(rxnID,'-','');
    rxnName=sprintf('biosynthesis %s from %s',cell2mat(rbosomeSubunit(i)),cell2mat(p(i,3)));
    
    model=addYeastReaction(model,r,{rxnID},{rxnName},0,1000,0,{''});

    %folding
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

subunit=unique(rbosomeSubunit);

eq=sprintf('%s [cytoplasm]', cell2mat(subunit(1)));
eq_deg=sprintf('%s_subunit [cytoplasm]', cell2mat(subunit(1)));

for i=2:numel(subunit)
    eq=sprintf('%s + %s [cytoplasm]',eq,cell2mat(subunit(i)));
    eq_deg=sprintf('%s + %s_subunit [cytoplasm]',eq_deg,cell2mat(subunit(i)));
end

eq = sprintf('%s => Ribosome [cytoplasm]',eq);
rxnID=sprintf('Ribosome_biosynthesis');

rxnName=sprintf('biosynthesis of ribosome');
 model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});

%degradation
eq = sprintf('Ribosome [cytoplasm] => %s',eq_deg);
rxnID=sprintf('Ribosome_degradation');

rxnName=sprintf('degradation of ribosome');
 model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
 
%dilution
eq = sprintf('Ribosome [cytoplasm] => ');
rxnID=sprintf('Ribosome_dilution');

rxnName=sprintf('dilution of ribosome');
 model=addYeastReaction(model,eq,{rxnID},{rxnName},0,1000,0,{''});
