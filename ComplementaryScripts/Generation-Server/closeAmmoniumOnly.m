function model = closeAmmoniumOnly(model,uptake_rate,c)

%% close Ammonium
model.ub(find(ismember(model.rxns,'r_1654_reverse')))=uptake_rate;
model.lb(find(ismember(model.rxns,'r_1654_reverse')))=0;
model.c(find(ismember(model.rxns,'r_1654_reverse')))=c;

