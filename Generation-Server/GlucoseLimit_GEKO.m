%This file generates glucose-limit growth
clear all
load NewYeastnextGeneration.mat
mu1=[0.025 0.05:0.05:0.4] ;
atpCost.gam=20;
 atpCost.ngam=.0;
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


for k=1:length(mu1)
     mu=mu1(k);
    ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
    [sol model complexFileName ]=modelGenerationGeko(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
   % Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d.lp.out',mu,T,0), ExpressionFileName);
   % [geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d.lp.out',mu,T,0),mu,0.045);
end
%model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
system(sprintf('./myResults1 %s','Max_Mu_GEKO_Glucose'),'-echo');
system('echo "Ribosome"; grep "X8079" *.out > ribosome_result.txt','-echo');
system('echo "MitoRibosome"; grep "X11289" *.out > mito_ribosome_result.txt','-echo');

%modelProteinAllocation(MyModel,'Max_Mu_Glucose_0.350000_T33_Copies_0.lp.out',mu,0.045)
