function [sol glc]=OtherCarbone_GECKO_batch_SC_proteome_cost_translation(k,translaction_activety)
%This file generates glucose-limit growth

load NewYeastnextGeneration.mat
SoplexVersion='Soplex-1.6';
r=sprintf('YFP_folding [cytoplasm] => YFP_subunit [cytoplasm]');
rxnID=sprintf('YFP_complex_degradation');
rxnName=sprintf('YFP degradation');
MyModel=addYeastReaction(MyModel,r,{rxnID},{rxnName},0,1000,0,{''});

[ carbone b c]=xlsread('other_carbon.xlsx','Data_minimum');


 atpCost.gam=30;
 atpCost.ngam=1;
  T=33;

 overExpressed.syn='YFP_folding_cytoplasm';
 overExpressed.dilultion='YFP_complex_dilution';
 overExpressed.degradation='YFP_complex_degradation';

 
copies =(k-1)*500000;    

eq=sprintf(' => %s',cell2mat(c(1,3)));
rxnID=sprintf('%s_transporter',cell2mat(c(1,1)));
rxnNames=sprintf('transporter');

new_MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});

fileName=sprintf('%s_batch_copies_%d.lp',cell2mat(c(1,1)),copies);



% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol glc  ]=growthSearchDynamic(new_MyModel, 'minimum',atpCost,0,T, copies, overExpressed, complex_list,compartment,0, 0.6,fileName,SoplexVersion,translaction_activety);

 
fileName=sprintf('Protein_Cost_SC%d.sol',k);
fptr=fopen(fileName,'w');
fprintf(fptr,'%d\t%f\n',copies,sol.mu);
fclose(fptr); 
