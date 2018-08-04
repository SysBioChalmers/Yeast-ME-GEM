%This file generates glucose-limit growth
clear all
load NewYeastnextGeneration.mat

atpCost.gam=20;
 atpCost.ngam=0.5;
 T1=[33 36 38];
 overExpressed.syn='mCherry_complex';
 overExpressed.dilultion='mCherry_complex_dilution';
 overExpressed.degradation='mCherry_degradation';

 T=33;%T=T1(k);

% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


copies =0;
%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol glc  ]=growthSearchDynamic(MyModel, 'minimum',atpCost,20,T, copies, overExpressed, complex_list,compartment,0, 0.4)
   % Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d.lp.out',mu,T,0), ExpressionFileName);
  %  [geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d.lp.out',mu,T,0),mu,0.045);

%model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
system(sprintf('./myResults1 %s','Max_Mu_GEKO_Glucose'),'-echo');
system('echo "Ribosome"; grep "X8079" *.out > ribosome_result.txt','-echo');
system('echo "MitoRibosome"; grep "X11289" *.out > mito_ribosome_result.txt','-echo');
system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');

