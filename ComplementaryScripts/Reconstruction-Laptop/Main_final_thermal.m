cc=0;
if cc==0
clear all
% Reading the model
model=readCbModel('yeast_7.6_cobra.xml');

newModel1=initModel_yeast();
[newModel1,reaction_list complex_list compartment substrate_list product_list  lb ub type_list ] = getComplexfromReaction(model,newModel1,true);

newModel1.ub(find(ismember(newModel1.rxns,'r_1714_reverse')))=4.224;
newModel1.lb(find(ismember(newModel1.rxns,'r_1714_reverse')))=4.224;
sol=solveLP(newModel1,true)

 newModel=newModel1;



model_protein=newModel;
newModel_one=newModel;
newModele=newModel_one;
newModel=model_protein;


 


%Adding the transporter reaction to nucleus
eq='S-adenosyl-L-methionine [cytoplasm] => S-adenosyl-L-methionine [nucleus]';
rxnID='SAM_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='S-adenosyl-L-methionine [nucleus] => S-adenosyl-L-methionine [cytoplasm]';
rxnID='SAM_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='S-adenosyl-L-homocysteine [cytoplasm] => S-adenosyl-L-homocysteine [nucleus]';
rxnID='SAH_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='S-adenosyl-L-homocysteine [nucleus] => S-adenosyl-L-homocysteine [cytoplasm]';
rxnID='SAH_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});


eq='NADP(+) [cytoplasm] => NADP(+) [nucleus]';
rxnID='NADP_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='NADP(+) [nucleus] => NADP(+) [cytoplasm]';
rxnID='NADP_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='NADPH [cytoplasm] => NADPH [nucleus]';
rxnID='NADPH_transporting_forward';
rxnNames=sprintf('SAM transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='NADPH [nucleus] => NADPH [cytoplasm]';
rxnID='NADPH_transporting_reverse';
rxnNames=sprintf('SAM transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='CTP [cytoplasm] => CTP [nucleus]';
rxnID='CTP_transporting_forward';
rxnNames=sprintf('CTP transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='CTP [nucleus] => CTP [cytoplasm]';
rxnID='CTP_transporting_reverse';
rxnNames=sprintf('CTP transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='UTP [cytoplasm] => UTP [nucleus]';
rxnID='UTP_transporting_forward';
rxnNames=sprintf('UTP transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='UTP [nucleus] => UTP [cytoplasm]';
rxnID='UTP_transporting_reverse';
rxnNames=sprintf('UTP transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [cytoplasm] => GTP [nucleus]';
rxnID='GTP_transporting_forward';
rxnNames=sprintf('GTP transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [nucleus] => GTP [cytoplasm]';
rxnID='GTP_transporting_reverse';
rxnNames=sprintf('GTP transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [cytoplasm] => GTP [nucleus]';
rxnID='GTP_transporting_forward';
rxnNames=sprintf('GTP transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='GTP [nucleus] => GTP [cytoplasm]';
rxnID='GTP_transporting_reverse';
rxnNames=sprintf('GTP transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='diphosphate [cytoplasm] => diphosphate [nucleus]';
rxnID='diphosphate_transporting_forward';
rxnNames=sprintf('diphosphate transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='diphosphate [nucleus] => diphosphate [cytoplasm]';
rxnID='diphosphate_transporting_reverse';
rxnNames=sprintf('diphosphate transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='phosphate [cytoplasm] => phosphate [nucleus]';
rxnID='phosphate_transporting_forward';
rxnNames=sprintf('phosphate transporting into nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});


eq='phosphate [nucleus] => phosphate [cytoplasm]';
rxnID='phosphate_transporting_reverse';
rxnNames=sprintf('phosphate transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='phosphate [nucleus] => phosphate [cytoplasm]';
rxnID='phosphate_transporting_reverse';
rxnNames=sprintf('phosphate transporting from nucleus');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='5''-deoxyadenosine [cytoplasm] => ';
rxnID='five_prime_deoxyadenosine_sink';
rxnNames=sprintf('5''-deoxyadenosine leaves the model as dead metabolite');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

%We need to recylcle Arr>p that is generaed from tRNA splicing
eq='Arr_GT_p [cytoplasm] + H2O [cytoplasm] => Arr1p [cytoplasm] + H+ [cytoplasm]';
rxnID='RHEA_14492';
rxnNames=sprintf('2,3''-cyclic-nucleotide 3''-phosphodiesterase');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{'YGR247W'});

eq='Arr1p [cytoplasm] + H2O [cytoplasm] => ADP-ribose [cytoplasm] + diphosphate [cytoplasm]';
rxnID='RHEA_14492';
rxnNames=sprintf('2,3''-cyclic-nucleotide 3''-phosphodiesterase');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{'YBR022W'});

eq='ADP-ribose [cytoplasm] + H2O [cytoplasm] => ribose-5-phosphate [cytoplasm] + AMP [cytoplasm] + 2 H+ [cytoplasm]';
rxnID='RHEA_10415';
rxnNames=sprintf('ADP-ribose diphosphatase');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{'YBR111C'});

%TC-AMP for t6A modification
eq='L-threonine [cytoplasm] + ATP [cytoplasm] + bicarbonate [cytoplasm] => TC-AMP [cytoplasm] + diphosphate [cytoplasm] + H2O [cytoplasm]';
rxnID='TC_AMP';
rxnNames=sprintf('threonylcarbamoyladenylate synthase');
newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,100,0,{''});

[n a ctRNA]=xlsread('tRNA_modification_position_1.xlsx','test_codon');
for i=1:numel(ctRNA(:,1))

    met=cell2mat(ctRNA(i,2));
    met0=cell2mat(ctRNA(i,1));

    eq=sprintf('%s => ',met);
    
    rxnID=sprintf('relase_tRNA%d',i);
    rxnNames=sprintf('threonylcarbamoyladenylate synthase');
    newModel=addYeastReaction(newModel,eq,{rxnID},{rxnNames},0,1000,0,{''});
 
  
end




end


newModel2= addProtein2Model_test(newModel,2,1294);%
%newModel3=addFoldingReaction(newModel2,2,1294);
newModel5=addRibosome(newModel2);
newModel6= addProteinDegradation2Model(newModel5,9,1294);
newModel7= addProteinDilution2Model(newModel6,2,1294);

newModel10=addMitoRibosome(newModel7);
newModel11=addAssemblyFactors(newModel10);

newModel13=addNewProteins(newModel11);
newModel14 =addTranslationFactors(newModel13);
newModel15=addProteasome(newModel14);
newModel16=addChaperone(newModel15);
newModel17=addMitoProtein2Model(newModel16);
newModel19= addMetabolicComplex(newModel17,complex_list,compartment,type_list);

newModel20 = addTOM(newModel19);
newModel20 =  addTIM22(newModel20);
newModel20 =  addTIM23(newModel20);

newModel21=addUnmodeledProtein(newModel20);
newModel22 = addtRNA2Model_test(newModel21);

newModel23= addMitoProteinDegradation2Model(newModel22)
newModel24=addMitotRNA2Model_test(newModel23)
MyModel=newModel24;

eq='GDP [cytoplasm] => GDP [mitochondrion]';
rxnID='GDP_transporting_forward';
rxnNames=sprintf('GDP transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='CMP [cytoplasm] => CMP [mitochondrion]';
rxnID='CMP_transporting_forward';
rxnNames=sprintf('CMP transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});
% 
% 
eq='UMP [cytoplasm] => UMP [mitochondrion]';
rxnID='UMP_transporting_forward';
rxnNames=sprintf('UMP transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='ADP [cytoplasm] => ADP [mitochondrion]';
rxnID='ADP_transporting_forward';
rxnNames=sprintf('ADP transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

%add dilution tRNA reactions
%TC-AMP for t6A modification
eq='L-threonine [mitochondrion] + ATP [mitochondrion] + bicarbonate [mitochondrion] => TC-AMP [mitochondrion] + diphosphate [mitochondrion] + H2O [mitochondrion]';
rxnID='TC_AMP_mitochondrion';
rxnNames=sprintf('threonylcarbamoyladenylate synthase in mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='prenyl diphosphate [cytoplasm] => prenyl diphosphate [mitochondrion]';
rxnID='prenyl_diphosphate_transporting_forward';
rxnNames=sprintf('prenyl diphosphate transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='L-cysteine [cytoplasm] => L-cysteine [mitochondrion]';
rxnID='L_cysteine_transporting_forward';
rxnNames=sprintf('L-cysteine transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='L-glutamine [cytoplasm] => L-glutamine [mitochondrion]';
rxnID='L_glutamine_transporting_forward';
rxnNames=sprintf('L-glutamine transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

eq='L-methionine [cytoplasm] => L-methionine [mitochondrion]';
rxnID='L_methionine_transporting_forward';
rxnNames=sprintf('L-methionine transporting into mitochondrion');
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

%add AARS
MyModel.grRules(15652)={'YDR341C'};
MyModel.grRules(15671)={'YDR341C'};
MyModel.grRules(15690)={'YDR341C'};
MyModel.grRules(15708)={'YDR341C'};
MyModel.grRules(15723)={'YDR341C'};
MyModel.grRules(15743)={'YDR341C'};
MyModel.grRules(15763)={'YDR341C'};

MyModel= addMetabolicComplex(MyModel,{'YOR335C','YBR121C','YHR011W'},{'mitochondrion','mitochondrion','mitochondrion'},[1 1 1]);



save YeastNewThermal

saveToExcel(MyModel,'YeastNewGenerationThermal.xlsx');

MyModel=addhemoglobin(MyModel);
MyModel=addFluorescentProtein(MyModel);
MyModel=addHsp104(MyModel);
MyModel=addHsp78(MyModel);
MyModel=addHsp60(MyModel);
save YeastnextGeneration

% MyModel= addProteinDilution2Model(MyModel,2,1756)
% MyModel=addProteinDegradation2ModelRatio(MyModel,2,1756)
saveToExcel(MyModel,'NewYeastnextGeneration.xlsx');

save NewYeastnextGeneration


