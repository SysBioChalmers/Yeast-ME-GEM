
newModel1.ub(find(ismember(newModel1.rxns,'r_1714_reverse')))=4.224;
newModel1.lb(find(ismember(newModel1.rxns,'r_1714_reverse')))=0;

model = newModel1;

model.ub(find(ismember(model.rxns,'r_1992_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1992_reverse')))=0;

model.S(strcmp(model.mets,'m2220'),strcmp(model.rxns,'r_4041_forward')) = 0;

eq=sprintf(' => ergosterol [extracellular]');
    rxnID='ergosterol_import';
    rxnNames=sprintf('ergosterol_import');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf(' => lanosterol [extracellular]');
    rxnID='lanosterol_import';
    rxnNames=sprintf('lanosterol import');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf(' => palmitoleate [extracellular]');
    rxnID='palmitoleate_import';
    rxnNames=sprintf('palmitoleate_import');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf(' => zymosterol [extracellular]');
    rxnID='zymosterol_import';
    rxnNames=sprintf('zymosterol_import');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf(' => 14-demethyllanosterol [extracellular]');
    rxnID='14demethyllanosterol_import';
    rxnNames=sprintf('14-demethyllanosterol_import');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf(' => ergosta-5,7,22,24(28)-tetraen-3beta-ol [extracellular]');
    rxnID='rgosta57_import';
    rxnNames=sprintf('ergosta-5,7,22,24(28)-tetraen-3beta-olimport');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    
    eq=sprintf(' => oleate [extracellular]');
    rxnID='oleate_import';
    rxnNames=sprintf('oleate_import');
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
    

sol=solveLP(model,true)