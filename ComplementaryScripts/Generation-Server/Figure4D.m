%This file generates glucose-limit growth
clc
clear all
load NewYeastnextGeneration.mat

mu1= [0.1:0.01:0.16] ;%
%mu1= [0.02:0.02:0.38] ;
atpCost.gam=33;
 atpCost.ngam=1;
 T1=[33 36 38];
 overExpressed.syn='mCherry_complex';
 overExpressed.dilultion='mCherry_complex_dilution';
 overExpressed.degradation='mCherry_degradation';
 mu=mu1;
 T=33;%T=T1(k);
PredictedProteinExpresstion=zeros(1520,numel(mu1));
ribosomeRatio=zeros(1,numel(mu1));

fileAbundance.fileName='growth_rate_0.2.xlsx';
fileAbundance.sheet='growth_rate_0.2';
sol={};
ss={};
copies=0;
i=1;

ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',0.1);
AT=0.5:0.05:1;%1:-0.05:0.5;

for active_translation=AT
for kk=1:1
   fileName=sprintf('%g_Translation_activation.lp',active_translation);

%    [sol{i} glc  ]=gro  wthSearchDynamic(MyModel, 'minimum',atpCost,0,T, copies, overExpressed, complex_list,compartment,0, 0.6,fileName,'Soplex-1.6',active_translation)
% [sol model complexFileName total_protein_molecules_model ]=modelGenerationGECKO	(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp',mu,T,0,kk),complex_list,compartment,'Soplex-1.6',kk,active_translation);
[sol{i} glc  ]=growthSearchDynamic(MyModel, 'minimum',atpCost,30,T, copies, overExpressed, complex_list,compartment,0, 0.6,sprintf('Translation_active_%f_GECKO_Glucose_T%d_Copies_%d_k_%d.lp',active_translation,T,0,kk),'Soplex-1.6',active_translation)
[ss{i} model1 complexFileName ]=modelGenerationGECKO(MyModel,'minimum', sol{i}.mu,atpCost,'min_carbon',30,T,0, overExpressed, fileName,complex_list,compartment,'Soplex-1.6',1,active_translation);

i=i+1;
end

end %translation

for i=1:length(AT)
    fprintf('\t%f\t%.15g\n',sol{i}.mu,sol{i}.X(4321));
end

for i=1:length(AT)
    fprintf('\t%f\t%.15g\n',sol{i}.mu,sol{i}.X(8079));
end

save Figure4D