function model = Glucose_Growth(model,max_uptake_rate,c)
% the aim of this function is to fix glucose growth media for the model

%% close glucose
model.ub(find(ismember(model.rxns,'r_1714_reverse')))=max_uptake_rate;
model.lb(find(ismember(model.rxns,'r_1714_reverse')))=0;
model.c(find(ismember(model.rxns,'r_1714_reverse')))=c;

%fixed pyruvate transporter with lactate
model.ub(find(ismember(model.rxns,'r_1138_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1138_reverse')))=0;



