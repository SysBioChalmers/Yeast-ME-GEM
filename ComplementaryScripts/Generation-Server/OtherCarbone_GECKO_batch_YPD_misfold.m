%This file generates glucose-limit growth
clear all
load NewYeastnextGeneration.mat

r=sprintf('YFP_folding [cytoplasm] => YFP_subunit [cytoplasm]');
rxnID=sprintf('YFP_complex_degradation');
rxnName=sprintf('YFP degradation');
MyModel=addYeastReaction(MyModel,r,{rxnID},{rxnName},0,1000,0,{''});


[ carbone b c]=xlsread('other_carbon.xlsx','Data_minimum');
sol={};

 atpCost.gam=30;
 atpCost.ngam=1;
  T=33;

 overExpressed.syn='YFP_folding_cytoplasm';
 overExpressed.dilultion='YFP_complex_dilution';
 overExpressed.degradation='YFP_complex_degradation';

 
parfor k=1:2

copies =2*(k-1)/10;    

eq=sprintf(' => %s',cell2mat(c(1,3)));
rxnID=sprintf('%s_transporter',cell2mat(c(1,1)));
rxnNames=sprintf('transporter');

new_MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

fileName=sprintf('%s_batch_copies_%d.lp',cell2mat(c(1,1)),copies);



% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol{k} glc  ]=growthSearchDynamic(new_MyModel, 'YPD',atpCost,0,T, copies, overExpressed, complex_list,compartment,0, 0.6,fileName)

 
end

