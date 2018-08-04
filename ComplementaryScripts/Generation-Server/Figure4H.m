%This file generates glucose-limit growth
%clear all
load NewYeastnextGeneration.mat

r=sprintf('YFP_folding [cytoplasm] => YFP_subunit [cytoplasm]');
rxnID=sprintf('YFP_complex_degradation');
rxnName=sprintf('YFP degradation');
MyModel=addYeastReaction(MyModel,r,{rxnID},{rxnName},0,1000,0,{''});


eq=sprintf(' => %s','D-galactose [extracellular]');
rxnID=sprintf('Galactose_transporter');
rxnNames=sprintf('test transporter');

new_MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,100,0,{''});



[ carbone b c]=xlsread('other_carbon.xlsx','Data_minimum');
sol={};

 atpCost.gam=33;
 atpCost.ngam=1;
  T=33;

 overExpressed.syn='YFP_folding_cytoplasm';
 overExpressed.dilultion='YFP_complex_dilution';
 overExpressed.degradation='YFP_complex_degradation';

n=20; 
m=1
start_mu=0;
for k=n:-1:m

copies =(k-1)*1e+06;    

fileName=sprintf('galactose_batch_copies_%d.lp',copies);




% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol{k} glc  ]=growthSearchDynamic(new_MyModel, 'minimum',atpCost,0,T, copies, overExpressed, complex_list,compartment,start_mu, 0.5,fileName,'Soplex-1.6',1)
[ss{k} model1 complexFileName ]=modelGenerationGECKO(MyModel,'minimum', sol{k}.mu,atpCost,'min_carbon',30,T,0, overExpressed, fileName,complex_list,compartment,'Soplex-1.6',1,1);
start_mu=sol{k}.mu;
 
end


for i=1:n
    fprintf('\t%f\t%.15g\n',sol{i}.mu,sol{i}.X(4321));
end

for i=1:n
    fprintf('\t%f\t%.15g\n',sol{i}.mu,sol{i}.X(8079));
end


save Protein_Cost_galactose
 
