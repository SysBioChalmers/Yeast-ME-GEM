%This file generates glucose-limit growth
%clc
clear all
load NewYeastnextGeneration.mat

mu1= 0.2;%0.1:0.01:0.2;%0.22:0.02:0.3;% [0.2 0.25];%[0.1 0.2];%[0.1:0.02:.28] ;%
%mu1= [0.02:0.02:0.38] ;
atpCost.gam=70; %return to 65
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
   
   
for k=1:length(mu1)

for kk=1:1
   mu=mu1(k);
   
 eq=sprintf('H+ [mitochondrion] + NADH [mitochondrion] + 0.5 oxygen [mitochondrion] => NAD [mitochondrion] + H2O [mitochondrion]');
 rxnID='AOX';
 rxnNames=sprintf('AOX reaction');
 MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,1000,0,{''});

%[sol glc  ]=growthSearchDynamicNH4(MyModel, 'minimum',atpCost,21,T, copies, overExpressed, complex_list,compartment,0.1, 0.4,fileName,'Soplex-1.6',active_translation)
[sol model complexFileName total_protein_molecules_model ]=modelGenerationGeko_NH4(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Std_Max_Mu_GEKO_Glucose_NH4_%f_T%d_Copies_%d_k_%d.lp',mu,T,0,kk),complex_list,compartment,'Soplex-1.6',kk,active_translation);
  % Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk), ExpressionFileName,total_protein_molecules_model,'Soplex-1.6');
  %[geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk),mu,0.045,'Soplex-1.6');
end
end

sol
% %model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
% system(sprintf('./myResults1 %s >translation_actvitiy%g','Std_Max_Mu_GEKO_Glucose_NH4',active_translation),'-echo');
% system('echo "Ribosome"; grep "X8079" Std_*.out > std_ribosome_result.txt','-echo');
% system('echo "MitoRibosome"; grep "X11289" *.out > std_mito_ribosome_result.txt','-echo');
% system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');
 system(sprintf('./myResults1 %s','Std_Max_Mu_GEKO_Glucose_NH4'),'-echo');
% 
% 
% %modelProteinAllocation(MyModel,'Max_Mu_Glucose_0.350000_T33_Copies_0.lp.out',mu,0.045)
% plot(mu1,ribosomeRatio)
save nitrogen-limit1