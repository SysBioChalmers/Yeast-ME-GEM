clear all
model=readCbModel('yeast_7.6_cobra.xml')
tRNA_model=initModel_yeast();
[tRNA_model,reaction_list complex_list compartment substrate_list product_list  lb ub type_list ] = getComplexfromReaction(model,tRNA_model,true);

tRNA_model.ub(find(ismember(tRNA_model.rxns,'r_1714_reverse')))=1;
tRNA_model.lb(find(ismember(tRNA_model.rxns,'r_1714_reverse')))=1;
solveLP(tRNA_model,true)

tRNA_model = addtRNA2Model(tRNA_model);

% Adding the transporter reaction to nucleus
eq='S-adenosyl-L-methionine [cytoplasm] => S-adenosyl-L-methionine [nucleus]';
rxnID='SAM_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='S-adenosyl-L-methionine [nucleus] => S-adenosyl-L-methionine [cytoplasm]';
rxnID='SAM_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='S-adenosyl-L-homocysteine [cytoplasm] => S-adenosyl-L-homocysteine [nucleus]';
rxnID='SAH_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='S-adenosyl-L-homocysteine [nucleus] => S-adenosyl-L-homocysteine [cytoplasm]';
rxnID='SAH_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});


eq='NADP(+) [cytoplasm] => NADP(+) [nucleus]';
rxnID='NADP_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='NADP(+) [nucleus] => NADP(+) [cytoplasm]';
rxnID='NADP_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='NADPH [cytoplasm] => NADPH [nucleus]';
rxnID='NADPH_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='NADPH [nucleus] => NADPH [cytoplasm]';
rxnID='NADPH_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='CTP [cytoplasm] => CTP [nucleus]';
rxnID='CTP_transporting_forward';
rxnNames=sprintf('CTP transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='CTP [nucleus] => CTP [cytoplasm]';
rxnID='CTP_transporting_reverse';
rxnNames=sprintf('CTP transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='UTP [cytoplasm] => UTP [nucleus]';
rxnID='UTP_transporting_forward';
rxnNames=sprintf('UTP transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='UTP [nucleus] => UTP [cytoplasm]';
rxnID='UTP_transporting_reverse';
rxnNames=sprintf('UTP transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [cytoplasm] => GTP [nucleus]';
rxnID='GTP_transporting_forward';
rxnNames=sprintf('GTP transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [nucleus] => GTP [cytoplasm]';
rxnID='GTP_transporting_reverse';
rxnNames=sprintf('GTP transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [cytoplasm] => GTP [nucleus]';
rxnID='GTP_transporting_forward';
rxnNames=sprintf('GTP transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [nucleus] => GTP [cytoplasm]';
rxnID='GTP_transporting_reverse';
rxnNames=sprintf('GTP transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='diphosphate [cytoplasm] => diphosphate [nucleus]';
rxnID='diphosphate_transporting_forward';
rxnNames=sprintf('diphosphate transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='diphosphate [nucleus] => diphosphate [cytoplasm]';
rxnID='diphosphate_transporting_reverse';
rxnNames=sprintf('diphosphate transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='phosphate [cytoplasm] => phosphate [nucleus]';
rxnID='phosphate_transporting_forward';
rxnNames=sprintf('phosphate transporting into nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});


eq='phosphate [nucleus] => phosphate [cytoplasm]';
rxnID='phosphate_transporting_reverse';
rxnNames=sprintf('phosphate transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='phosphate [nucleus] => phosphate [cytoplasm]';
rxnID='phosphate_transporting_reverse';
rxnNames=sprintf('phosphate transporting from nucleus');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='5''-deoxyadenosine [cytoplasm] => ';
rxnID='five_prime_deoxyadenosine_sink';
rxnNames=sprintf('5''-deoxyadenosine leaves the model as dead metabolite');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

%We need to recylcle Arr>p that is generaed from tRNA splicing
eq='Arr_GT_p [cytoplasm] + H2O [cytoplasm] => Arr1p [cytoplasm] + H+ [cytoplasm]';
rxnID='RHEA_14492';
rxnNames=sprintf('2,3''-cyclic-nucleotide 3''-phosphodiesterase');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{'YGR247W'});

eq='Arr1p [cytoplasm] + H2O [cytoplasm] => ADP-ribose [cytoplasm] + diphosphate [cytoplasm]';
rxnID='RHEA_14492';
rxnNames=sprintf('2,3''-cyclic-nucleotide 3''-phosphodiesterase');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{'YBR022W'});

eq='ADP-ribose [cytoplasm] + H2O [cytoplasm] => ribose-5-phosphate [cytoplasm] + AMP [cytoplasm] + 2 H+ [cytoplasm]';
rxnID='RHEA_10415';
rxnNames=sprintf('ADP-ribose diphosphatase');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{'YBR111C'});

%TC-AMP for t6A modification
eq='L-threonine [cytoplasm] + ATP [cytoplasm] + bicarbonate [cytoplasm] => TC-AMP [cytoplasm] + diphosphate [cytoplasm] + H2O [cytoplasm]';
rxnID='TC_AMP';
rxnNames=sprintf('threonylcarbamoyladenylate synthase');
tRNA_model=addYeastReaction(tRNA_model,eq,{rxnID},{rxnNames},0,100,0,{''});

saveToExcel(tRNA_model,'model_test_protein_tRNA_only.xlsx');