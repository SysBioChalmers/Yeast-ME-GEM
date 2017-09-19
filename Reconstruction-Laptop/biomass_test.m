model=MyModel;
model.c=zeros(numel(model.c),1);
eq='40 H2O [cytoplasm] + 40 ATP [cytoplasm] + lipid [cytoplasm] + 1.14 (1-_GT_3)-beta-D-glucan [cytoplasm] + 0.821 mannan [cytoplasm] + 0.519 glycogen [cytoplasm] + 0.0234 trehalose [cytoplasm] + 0.051 AMP [cytoplasm] + 0.00359 dAMP [cytoplasm] + 0.051 GMP [cytoplasm] + 0.05 CMP [cytoplasm] + 0.067 UMP [cytoplasm] + 0.00243 dCMP [cytoplasm] + 0.00243 dGMP [cytoplasm] + 0.00359 dTMP [cytoplasm] => 40 H+ [cytoplasm] + 40 phosphate [cytoplasm] + 40 ADP [cytoplasm] + biomass_test [cytoplasm]';
rxnID='biomass_test';
rxnNames=sprintf('biomass test');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,1,{''});
sol=testMetabolite(model,'biomass_test [cytoplasm]')