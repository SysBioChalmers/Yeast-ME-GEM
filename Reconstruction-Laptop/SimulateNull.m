load nullmodel.mat;
model.ub(find(ismember(model.rxns,'r_1714_reverse')))=max_uptake_rate;
model.lb(find(ismember(model.rxns,'r_1714_reverse')))=0;

%max glucose

