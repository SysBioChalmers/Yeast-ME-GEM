function model= fixReaction(model)
%prevent H+ export without cost
model.ub(find(ismember(model.rxns,'r_1824_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1824_reverse')))=0;
  
%
model.ub(find(ismember(model.rxns,'r_1250_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1250_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1259_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1259_forward')))=0;

 model.ub(find(ismember(model.rxns,'r_2033_forward')))=0;
 model.lb(find(ismember(model.rxns,'r_2033_forward')))=0;

%r_2083_forward	tryptophol [extracellular] => 
%model.ub(find(ismember(model.rxns,'r_2083_forward')))=0;
%model.lb(find(ismember(model.rxns,'r_2083_forward')))=0;