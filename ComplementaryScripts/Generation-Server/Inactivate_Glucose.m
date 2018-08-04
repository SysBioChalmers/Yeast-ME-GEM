function model = Inactivate_Glucose(model)
% %close ADP => ATP

%H+ 
model.ub(find(ismember(model.rxns,'r_1824_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1824_reverse')))=0;

%Pyruvate is transported without lactate
model.ub(find(ismember(model.rxns,'r_1138_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1138_reverse')))=0;


% % bicarbonate for only CO2 exporting
 model.ub(find(ismember(model.rxns,'r_1663_forward')))=0;
 model.lb(find(ismember(model.rxns,'r_1663_forward')))=0;
 
 
% %% Rui
model.ub(find(ismember(model.rxns,'r_2116_forward')))=0;
model.lb(find(ismember(model.rxns,'r_2116_forward')))=0;
% 
model.ub(find(ismember(model.rxns,'r_0659_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0659_forward')))=0;
% 
% %% fixed AA exporter
 model.ub(find(ismember(model.rxns,'r_1810_forward')))=0;
 model.lb(find(ismember(model.rxns,'r_1810_forward')))=0;
 
%prevent (R)-pantothenate exporting
model.ub(find(ismember(model.rxns,'r_1548_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1548_forward')))=0;

% BAT
model.ub(find(ismember(model.rxns,'r_1097_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1097_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1572_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1572_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1865_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1865_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1598_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1598_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1899_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1899_forward')))=0;


model.ub(find(ismember(model.rxns,'r_1914_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1914_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1631_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1631_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1862_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1862_forward')))=0;

%Block L-alanine transporter
model.ub(find(ismember(model.rxns,'r_1873_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1873_forward')))=0;

%block glycerol
%Block GCY1 
model.ub(find(ismember(model.rxns,'r_0487_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0487_forward')))=0;

model.ub(find(ismember(model.rxns,'r_0472_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0472_forward')))=0;


model.ub(find(ismember(model.rxns,'r_0713_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_0713_reverse')))=0;


model.ub(find(ismember(model.rxns,'r_0714_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_0714_reverse')))=0;





