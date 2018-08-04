%This file generates glucose-limit growth
clear all
load NewYeastnextGeneration.mat

[ carbone b c]=xlsread('other_carbon.xlsx','Data_minimum');
sol={};

 atpCost.gam=33;
 atpCost.ngam=1;
  T=33;

 overExpressed.syn='mCherry_complex';
 overExpressed.dilultion='mCherry_complex_dilution';
 overExpressed.degradation='mCherry_degradation';

 
copies =0;



parfor i=1:1%numel(c(:,1))
 % add carbon transporter
 
eq=sprintf(' => %s',cell2mat(c(i,3)));
rxnID=sprintf('%s_transporter',cell2mat(c(i,1)));
rxnNames=sprintf('test biomass');

new_MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

fileName=sprintf('%s_batch.lp',cell2mat(c(i,1)));



% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol{i} glc  ]=growthSearchDynamic(new_MyModel, 'YPD',atpCost,0,T, copies, overExpressed, complex_list,compartment,0, 0.6,fileName)

 
end

%minmum uptake
parfor i=1:1%numel(c(:,1))
 % add carbon transporter
 
eq=sprintf(' => %s',cell2mat(c(i,3)));
rxnID=sprintf('%s_transporter',cell2mat(c(i,1)));
rxnNames=sprintf('test biomass');

new_MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,cell2mat(c(i,4)),-1,{''});

fileName=sprintf('%s_batch.lp',cell2mat(c(i,1)));


[ss model complexFileName ]=modelGenerationGeko(new_MyModel,'YPD', sol{i}.mu,atpCost,'min_carbon',0,T,0, overExpressed, fileName,complex_list,compartment);

 
end

% %model=modelGenerationNew(MyModel,'minimum', mu,atpCost,'min_glucose',22.6,T,0, overExpressed, sprintf('Max_Mu_Glucose_%f_T%d_Copies_%d.lp',mu,T,0),complex_list,compartment);
% system(sprintf('./myResults1 %s','Max_Mu_GEKO_Glucose'),'-echo');
% system('echo "Ribosome"; grep "X8079" *.out > ribosome_result.txt','-echo');
% system('echo "MitoRibosome"; grep "X11289" *.out > mito_ribosome_result.txt','-echo');
% system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');
for i=1:1%numel(c(:,1))
    fprintf('%s\t%f\n',cell2mat(c(i,1)),sol{i}.mu);
end
