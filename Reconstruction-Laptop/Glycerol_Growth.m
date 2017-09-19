function model = Glycerol_Growth(model,max_uptake_rate,c)
% the aim of this function is to fix ethanol growth media for the model
% 

%% close glucose
model.ub(find(ismember(model.rxns,'r_1714_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1714_reverse')))=0;
model.c(find(ismember(model.rxns,'r_1714_reverse')))=0;

%% ADD Glycerol Importer
eq=sprintf(' => glycerol [extracellular]');
rxnID='Glycerol_importer';
rxnNames=sprintf('Glycerol importer');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,max_uptake_rate,c,{''});

%% open glyoxylate shunt
model.ub(find(ismember(model.rxns,'r_0662_forward')))=100;
model.lb(find(ismember(model.rxns,'r_0662_forward')))=0;