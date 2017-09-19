function model=addhemoglobin(model)

%sequence extraction
p=fastaread('HemoglobinSequence.txt');
HBA.seq=p(1).Sequence;
HBA.name= 'HBA';
HBA.UTR5=5;
HBA.MW= 15258;

model= addProtein2ModelFromSequence(model,HBA);
model= addProteinDegradation2ModelFromSequence(model,HBA);

HBB.seq = p(2).Sequence;
HBB.name= 'HBB';
HBB.UTR5=5;
HBB.MW=15998;

%total MW
MW=2*HBA.MW + 2*HBB.MW; %= 62512

%protein translation and degradation
model= addProtein2ModelFromSequence(model,HBB);
model= addProteinDegradation2ModelFromSequence(model,HBB);

%folding

r=sprintf('HBA_peptide [cytoplasm] => HBA_folding [cytoplasm]');
rxnID=sprintf('HBA_folding_cytoplasm');
rxnName=sprintf('HBA folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
%misfolding
r=sprintf('HBA_folding [cytoplasm] => HBA_subunit [cytoplasm]');
rxnID=sprintf('HBA_degradation_misfolding_cytoplasm');
rxnName=sprintf('HBA misfolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});


r=sprintf('HBB_peptide [cytoplasm] => HBB_folding [cytoplasm]');
rxnID=sprintf('HBB_degradation_misfolding_cytoplasm');
rxnName=sprintf('HBB folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('HBB_peptide [cytoplasm] => HBB_subunit [cytoplasm]');
rxnID=sprintf('HBB_misfolding_cytoplasm');
rxnName=sprintf('HBB misfolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});


%Hemoglobin
r=sprintf('2 HBA_folding [cytoplasm] + 2 HBB_folding [cytoplasm] + 4 ferroheme b [cytoplasm] => Hemoglobin [cytoplasm]');
rxnID=sprintf('Hemoglobin_complex_formation');
rxnName=sprintf('Hemoglobin complex formation');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('Hemoglobin [cytoplasm] => 2 HBA_subunit [cytoplasm] + 2 HBB_subunit [cytoplasm] + 4 ferroheme b [cytoplasm]');
rxnID=sprintf('Hemoglobin_complex_degradation');
rxnName=sprintf('Hemoglobin complex degradation');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('Hemoglobin [cytoplasm] => ');
rxnID=sprintf('Hemoglobin_complex_dilution');
rxnName=sprintf('Hemoglobin complex dilution');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('ferroheme b [mitochondrion] => ferroheme b [cytoplasm]');
rxnID=sprintf('ferroheme_b_transport');
rxnName=sprintf('ferroheme_b_transport');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
