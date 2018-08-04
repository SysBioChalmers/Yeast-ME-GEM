%This file generates glucose-limit growth
clc
clear all
load NewYeastnextGeneration.mat

mu1= 0.1;%[0.1:0.02:.2] ;%
%mu1= [0.02:0.02:0.38] ;
atpCost.gam=40;
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
copies=0;
active_translation=1;

ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
fileName=sprintf('%g_Translation_activation.lp',active_translation);
[ a b c]=xlsread('mutant.xlsx');
gene=c(:,2);
sol={''};
for k=1:1%numel(gene)

for kk=1:1
   mu=mu1(k);
   model_mutant=MyModel;
   gene_reaction_name = sprintf('%s_translation', cell2mat(gene(k)));
   index=find(ismember( model_mutant.rxns,gene_reaction_name));
%    model_mutant.lb(index)=0;
%    model_mutant.ub(index)=0;

[sol1 glc  ]=growthSearchDynamic_mutant(model_mutant, 'minimum',atpCost,21,T, copies, overExpressed, complex_list,compartment,0.1, 0.4,fileName,'Soplex-1.6',active_translation) 
[sol{k} model complexFileName total_protein_molecules_model ]=modelGenerationGeko_mutant(model_mutant,'minimum', sol1.mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Std_Max_Mu_GEKO_Glucose_NH4_%f_T%d_Copies_%d_k_%d.lp',mu,T,0,kk),complex_list,compartment,'Soplex-3.0',kk,active_translation);

end
end

sol
save mutant