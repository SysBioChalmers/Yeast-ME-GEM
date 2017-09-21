function model=addFluorescentProtein(model)

%sequence extraction
p=fastaread('YFP.txt');
YFP.seq=p(1).Sequence;
YFP.name= 'YFP';
YFP.UTR5=5;
YFP.MW= 15258;

model= addProtein2ModelFromSequence(model,YFP);
model= addProteinDegradation2ModelFromSequence(model,YFP);

%GFP
p=fastaread('GFP.txt');
GFP.seq = p.Sequence;
GFP.name= 'GFP';
GFP.UTR5=5;
GFP.MW=15998;

model= addProtein2ModelFromSequence(model,GFP);
model= addProteinDegradation2ModelFromSequence(model,GFP);

%mcherry
p=fastaread('mCherry.txt');
mCherry.seq = p.Sequence;
mCherry.name= 'mCherry';
mCherry.UTR5=5;
mCherry.MW=15998;

model= addProtein2ModelFromSequence(model,mCherry);
model= addProteinDegradation2ModelFromSequence(model,mCherry);

%mcherry intial cost
p=fastaread('mCherry.txt');
mCherry.seq = p.Sequence;
mCherry.name= 'mCherry1';
mCherry.UTR5=20;
mCherry.MW=15998;

model= addProtein2ModelFromSequence(model,mCherry);
model= addProteinDegradation2ModelFromSequence(model,mCherry);
%YFP
%folding

r=sprintf('YFP_peptide [cytoplasm] => YFP_folding [cytoplasm]');
rxnID=sprintf('YFP_folding_cytoplasm');
rxnName=sprintf('YFP folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('YFP_folding [cytoplasm] => ');
rxnID=sprintf('YFP_complex_dilution');
rxnName=sprintf('YFP dilution');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});


%misfolding
r=sprintf('YFP_folding [cytoplasm] => YFP_misfolding [cytoplasm]');
rxnID=sprintf('YFP_misfolding_cytoplasm');
rxnName=sprintf('YFP misfolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%refolding
r=sprintf('YFP_misfolding [cytoplasm] => YFP_folding [cytoplasm]');
rxnID=sprintf('YFP_refolding_cytoplasm');
rxnName=sprintf('YFP refolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});


%misfolding dilution
r=sprintf('YFP_folding [cytoplasm] => YFP_subunit [cytoplasm]');
rxnID=sprintf('YFP_degradation_misfolding_cytoplasm');
rxnName=sprintf('YFP misfolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%misfolding dilution
r=sprintf('YFP_misfolding [cytoplasm] => ');
rxnID=sprintf('YFP_dilution_misfolding_cytoplasm');
rxnName=sprintf('YFP misfolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});


%GFP
r=sprintf('GFP_peptide [cytoplasm] => GFP_folding [cytoplasm]');
rxnID=sprintf('GFP_cytoplasm');
rxnName=sprintf('GFP folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('GFP_folding [cytoplasm] => ');
rxnID=sprintf('GFP_complex_dilution');
rxnName=sprintf('GFP dilution');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%misfolding
r=sprintf('GFP_folding [cytoplasm] => GFP_subunit [cytoplasm]');
rxnID=sprintf('GFP_complex_degradation');
rxnName=sprintf('GFP misfolding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%mCherry
r=sprintf('mCherry_peptide [cytoplasm] => mCherry_folding [cytoplasm]');
rxnID=sprintf('mCherry_folding_cytoplasm');
rxnName=sprintf('mCherry folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});


r=sprintf('mCherry_folding [cytoplasm] => ');
rxnID=sprintf('mCherry_complex_dilution');
rxnName=sprintf('mCherry dilution');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%misfolding
r=sprintf('mCherry_folding [cytoplasm] => mCherry_subunit [cytoplasm]');
rxnID=sprintf('mCherry_complex_degradation');
rxnName=sprintf('mCherry complex degradation');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%mCherry1
r=sprintf('mCherry1_peptide [cytoplasm] => mCherry1_folding [cytoplasm]');
rxnID=sprintf('mCherry1_folding_cytoplasm');
rxnName=sprintf('mCherry1 folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

r=sprintf('mCherry1_folding [cytoplasm] => ');
rxnID=sprintf('mCherry1_complex_dilution');
rxnName=sprintf('mCherry1 dilution');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

%misfolding
r=sprintf('mCherry1_folding [cytoplasm] => mCherry1_subunit [cytoplasm]');
rxnID=sprintf('mCherry1_degradation_misfolding_cytoplasm');
rxnName=sprintf('mCherry folding');
model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

