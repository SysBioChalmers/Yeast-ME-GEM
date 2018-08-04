%This file generates glucose-limit growth
%clc
clear all
load NewYeastnextGeneration.mat

%mu1=[ 0.025 0.05:0.05:0.35 0.38 0.28];%[0.025 0.05:0.05:0.35 0.38  0.28];%[0.25:0.01:0.4]% 0.32 0.34 0.35 0.36 0.38 0.4] ;%
mu1= [0.02:0. 02:0.38] ;
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
kk=1;


for k=1:length(mu1)

   mu=mu1(k);
   model=MyModel;
   model.grRules=strrep(model.grRules,'Proteasome','');
   model.grRules=strrep(model.grRules,'Ribosome','');
   model.grRules=strrep(model.grRules,'TOM-TIM','');
   model.grRules=strrep(model.grRules,'MitoDegradome','');
   

   
% eq=sprintf('H+ [mitochondrion] + NADH [mitochondrion] + 0.5 oxygen [mitochondrion] => NAD [mitochondrion] + H2O [mitochondrion]');
% rxnID='AOX';
% rxnNames=sprintf('AOX reaction');
% model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});



 ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
 [sol model1 complexFileName total_protein_molecules_model ]=modelGenerationGECKO(model,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('AOX_Ribosome_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp',mu,T,0,kk),complex_list,compartment,'Soplex-1.6',kk,1);
%    Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk), ExpressionFileName,total_protein_molecules_model,'Soplex-3.0',kk);
 [geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('AOX_Ribosome_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk),mu,0.045,'Soplex-1.6');
end

%model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Ribosome_Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
system(sprintf('./myResults1 %s','AOX_Ribosome_Max_Mu_GEKO_Glucose'),'-echo');
system('echo "Ribosome"; grep "X8079" *.out > ribosome_result.txt','-echo');
system('echo "MitoRibosome"; grep "X11289" *.out > mito_ribosome_result.txt','-echo');
system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');


%modelProteinAllocation(MyModel,'Max_Mu_Glucose_0.350000_T33_Copies_0.lp.out',mu,0.045)
plot(mu1,ribosomeRatio)