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

%acetate
%model.ub(find(ismember(model.rxns,'r_0216_reverse')))=0;
%model.lb(find(ismember(model.rxns,'r_0216_reverse')))=0;





% %NADH regenerate
model.ub(find(ismember(model.rxns,'r_0182_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0182_forward')))=0;

model.ub(find(ismember(model.rxns,'r_0183_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0183_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1870_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1870_forward')))=0;
% 
 model.ub(find(ismember(model.rxns,'r_1087_forward')))=0;
 model.lb(find(ismember(model.rxns,'r_1087_forward')))=0;
%

model.ub(find(ismember(model.rxns,'r_1989_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1989_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1815_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1815_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1686_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_1686_reverse')))=0;




% 
% %adjust TCA in high growth rate
% 
% model.ub(find(ismember(model.rxns,'r_0002_forward')))=0;
% model.lb(find(ismember(model.rxns,'r_0002_forward')))=0;
% 
% model.ub(find(ismember(model.rxns,'r_0001_forward')))=0;
% model.lb(find(ismember(model.rxns,'r_0001_forward')))=0;
% 
% model.ub(find(ismember(model.rxns,'r_0004_forward')))=0;
% model.lb(find(ismember(model.rxns,'r_0004_forward')))=0;
% 
% % FMN-dependent NAD(P)H
% %suc
% model.ub(find(ismember(model.rxns,'r_2056_forward')))=0;
% model.lb(find(ismember(model.rxns,'r_2056_forward')))=0;
% 
% %glu import
% 
