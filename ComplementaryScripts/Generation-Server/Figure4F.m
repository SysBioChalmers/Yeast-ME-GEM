%This file generates glucose-limit growth
clc
clear all
load NewYeastnextGeneration.mat



mu1= [0.2:0.01:0.36] ;%
atpCost.gam=33;
 atpCost.ngam=1;
 T1=[33 36 38];
  
 overExpressed.syn='YOR100C_complex_formation_mitochondrion';
 overExpressed.dilultion='YOR100C_complex_dilution_mitochondrion';
 overExpressed.degradation='YOR100C_complex_degradation_mitochondrion';
 
 mu=mu1;
 T=33;%T=T1(k);
PredictedProteinExpresstion=zeros(1520,numel(mu1));
ribosomeRatio=zeros(1,numel(mu1));

fileAbundance.fileName='growth_rate_0.2.xlsx';
fileAbundance.sheet='growth_rate_0.2';
sol={};
for c=3:10
 copies=c*1e5;
active_translation=1;
for k=1:length(mu1)

for kk=1:1
   mu=mu1(k);
   ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
   fileName=sprintf('%g_Translation_activation.lp',active_translation);

%[sol glc  ]=growthSearchDynamicNH4(MyModel, 'minimum',atpCost,21,T, copies, overExpressed, complex_list,compartment,0, 0.15,fileName,'Soplex-1.6',active_translation)
[sol model complexFileName total_protein_molecules_model ]=modelGenerationGECKO(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,copies, overExpressed, sprintf('CRC1_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp',mu,T,copies,kk),complex_list,compartment,'Soplex-1.6',kk,active_translation);
  % Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk), ExpressionFileName,total_protein_molecules_model,'Soplex-1.6');
  %[geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk),mu,0.045,'Soplex-1.6');
end
end
system(sprintf('./myResults1 %s >CRC1%g','CRC1_Mu_GEKO',copies),'-echo');

end
