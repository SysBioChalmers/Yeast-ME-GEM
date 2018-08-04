%This file generates glucose-limit growth
%clc
clear all
load NewYeastnextGeneration.mat

mu1= 0.38 ;%[0.025 0.03:0.01:0.4];%[0.25:0.01:0.4]% 0.32 0.34 0.35 0.36 0.38 0.4] ;%
%mu1= [0.02:0. 02:0.38] ;
atpCost.gam=30;
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
gene={'YIL155C' 'YPL134C' 'YML054C' };
for k=1:length(mu1)

for kk=1:1
   mu=mu1(k);
   model=MyModel;
   for gg=1:numel(gene)
       gene_reaction_name = sprintf('%s_translation', cell2mat(gene(gg)));
       index=find(ismember( model.rxns,gene_reaction_name));
       model.lb(index)=0;
       model.ub(index)=0;
    end
% %    
%    % close oxaloacetate-malate shuttle
  model.ub(find(ismember(model.rxns,'r_0714_reverse')))=0.1;
  model.lb(find(ismember(model.rxns,'r_0714_reverse')))=0.1;
% 
%   model.ub(find(ismember(model.rxns,'r_0713_reverse')))=0;
%   model.lb(find(ismember(model.rxns,'r_0713_reverse')))=0;
%  
%glyceron
model.ub(find(ismember(model.rxns,'r_0487_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0487_forward')))=0;

%
model.ub(find(ismember(model.rxns,'r_0472_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0472_forward')))=0;

model.ub(find(ismember(model.rxns,'r_0472_reverse')))=0;
model.lb(find(ismember(model.rxns,'r_0472_reversed')))=0;

%FMN
model.ub(find(ismember(model.rxns,'r_1000_forward')))=0;
model.lb(find(ismember(model.rxns,'r_1000_forward')))=0;
% 
% %ODC
model.ub(find(ismember(model.rxns,'r_2132_forward')))=10;
model.lb(find(ismember(model.rxns,'r_2132_forward')))=0;

model.ub(find(ismember(model.rxns,'r_1126_reverse')))=10;
model.lb(find(ismember(model.rxns,'r_1126_reverse')))=0;
% 
% 
% %g3p [m]
% 
model.ub(find(ismember(model.rxns,'r_0492_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0492_forward')))=0;

model.ub(find(ismember(model.rxns,'r_0169_forward')))=0;
model.lb(find(ismember(model.rxns,'r_0169_forward')))=0;

%citrate transporter
% model.ub(find(ismember(model.rxns,'r_1126_reverse')))=0.2;
% model.lb(find(ismember(model.rxns,'r_1126_reverse')))=0;

% 
% %glycerol
% model.ub(find(ismember(model.rxns,'r_1808_forward')))=0.04;
% model.lb(find(ismember(model.rxns,'r_1808_forward')))=0.04;
% % 
% %  
% model.ub(find(ismember(model.rxns,'r_0770_forward')))=2;
% model.lb(find(ismember(model.rxns,'r_0770_forward')))=0;
% 
   
% eq=sprintf('H+ [mitochondrion] + NADH [mitochondrion] + 0.5 oxygen [mitochondrion] => NAD [mitochondrion] + H2O [mitochondrion]');
% rxnID='AOX';
% rxnNames=sprintf('AOX reaction');
% model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});



ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol model1 complexFileName total_protein_molecules_model ]=modelGenerationGeko_glycerol(model,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('GLYcerol_Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp',mu,T,0,kk),complex_list,compartment,'Soplex-1.6',kk,1);
%    Sfactor= Estimate_abundance(model,mu, complexFileName,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk), ExpressionFileName,total_protein_molecules_model,'Soplex-3.0',kk);
%    [geneList PredictedProteinExpresstion(:,k) ribosomeRatio(k)] =modelProteinAllocation(model,sprintf('Max_Mu_GEKO_Glucose_%f_T%d_Copies_%d_k_%d.lp.out',mu,T,0,kk),mu,0.045,'Soplex-3.0');
end
end
%model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
system(sprintf('./myResults1 %s','GLYcerol_Max_Mu_GEKO_Glucose'),'-echo');
system('echo "Ribosome"; grep "X8079" *.out > ribosome_result.txt','-echo');
system('echo "MitoRibosome"; grep "X11289" *.out > mito_ribosome_result.txt','-echo');
system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');


%modelProteinAllocation(MyModel,'Max_Mu_Glucose_0.350000_T33_Copies_0.lp.out',mu,0.045)
%plot(mu1,ribosomeRatio)