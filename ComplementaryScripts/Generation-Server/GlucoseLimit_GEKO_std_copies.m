%This file generates glucose-limit growth
clc
clear all
load NewYeastnextGeneration.mat

r=sprintf('YFP_folding [cytoplasm] => YFP_subunit [cytoplasm]');
rxnID=sprintf('YFP_complex_degradation');
rxnName=sprintf('YFP degradation');
MyModel=addYeastReaction(MyModel,r,{rxnID},{rxnName},0,1000,0,{''});


mu1= [0.1:0.02:0.3] ;%
%mu1= [0.02:0.02:0.38] ;
atpCost.gam=33;
 atpCost.ngam=1;
 T1=[33 36 38];
  overExpressed.syn='YFP_folding_cytoplasm';
 overExpressed.dilultion='YFP_complex_dilution';
 overExpressed.degradation='YFP_complex_degradation';
 
 mu=mu1;
 T=33;%T=T1(k);
PredictedProteinExpresstion=zeros(1520,numel(mu1));
ribosomeRatio=zeros(1,numel(mu1));

fileAbundance.fileName='growth_rate_0.2.xlsx';
fileAbundance.sheet='growth_rate_0.2';
sol={};
for c=1:2
 copies=c*1e6;
active_translation=1;
for k=1:length(mu1)

for kk=1:1
   mu=mu1(k);
   ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
   fileName=sprintf('%g_Translation_activation.lp',active_translation);

%[sol glc  ]=growthSearchDynamicNH4(MyModel, 'minimum',atpCost,21,T, copies, overExpressed, complex_list,compartment,0, 0.15,fileName,'Soplex-1.6',active_translation)
[sol model complexFileName total_protein_molecules_model ]=modelGenerationGECKO(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,copies, overExpressed, sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp',mu,T,copies,kk),complex_list,compartment,'Soplex-1.6',kk,active_translation);
  % Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk), ExpressionFileName,total_protein_molecules_model,'Soplex-1.6');
  %[geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('Std_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk),mu,0.045,'Soplex-1.6');
end
end

end
%model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
system(sprintf('./myResults1 %s >translation_actvitiy%g','Std_Max_Mu_GEKO_Glucose',active_translation),'-echo');
system('echo "Ribosome"; grep "X8079" Std_*.out > std_ribosome_result.txt','-echo');
system('echo "MitoRibosome"; grep "X11289" *.out > std_mito_ribosome_result.txt','-echo');
system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');
system(sprintf('./myResults1 %s','Std_Max_Mu_GEKO_Glucose'),'-echo');


%modelProteinAllocation(MyModel,'Max_Mu_Glucose_0.350000_T33_Copies_0.lp.out',mu,0.045)
plot(mu1,ribosomeRatio)