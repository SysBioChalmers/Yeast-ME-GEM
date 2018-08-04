%This file generates glucose-limit growth
clear all
load NewYeastnextGeneration.mat
%%Fructose transporter
eq=sprintf('D-fructose [extracellular] => D-fructose [cytoplasm]');
rxnID=sprintf('fructose_diffusion');
rxnNames=rxnID;
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,1000,0,{'YDR342C'}); %cell2mat(metabolite(i,3))
    

[ carbone b c]=xlsread('other_carbon.xlsx');
sol={};
ss={};
NN=numel(c(:,1));
MM=1;
%NN=1;
 atpCost.gam=33;
 atpCost.ngam=1;
  T=33;

 overExpressed.syn='mCherry_complex';
 overExpressed.dilultion='mCherry_complex_dilution';
 overExpressed.degradation='mCherry_degradation';

 
copies =0;

%adjust media
[a b metabolite]=xlsread('other_carbon.xlsx','Medium');
model=MyModel;
for i=2:numel(metabolite(:,1))
    eq=sprintf(' => %s',cell2mat(metabolite(i,2)));
    rxnID=sprintf('%s_import',cell2mat(metabolite(i,1)));
    rxnNames=rxnID;
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,cell2mat(metabolite(i,3)),0,{''}); %
end

   model.grRules=strrep(model.grRules,'Proteasome','');
   model.grRules=strrep(model.grRules,'Ribosome','');
   model.grRules=strrep(model.grRules,'TOM-TIM','');
   model.grRules=strrep(model.grRules,'MitoDegradome','');

for i=MM:NN % add carbon transporter
 
eq=sprintf(' => %s',cell2mat(c(i,3)));
rxnID=sprintf('%s_transporter',cell2mat(c(i,1)));
rxnNames=sprintf('transporter');

new_MyModel=addYeastReaction(model,eq,{rxnID},{rxnNames},0,cell2mat(c(i,4)),-1,{''});

fileName=sprintf('%s_batch.lp',cell2mat(c(i,1)));



% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol{i} glc  PredictedProteinMass newModel]=growthSearchDynamic(new_MyModel, 'minimum',atpCost,0,T, copies, overExpressed, complex_list,compartment,0.0, 0.6,fileName,'Soplex-1.6',1,'')
[ss{i} model1 complexFileName ]=modelGenerationGECKO(new_MyModel,'minimum', sol{i}.mu,atpCost,'min_carbon',0,T,0, overExpressed, fileName,complex_list,compartment,'Soplex-1.6',1,1);

 
end


% system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');
for i=MM:NN
    fprintf('%s\t%f\t%.15g\n',cell2mat(c(i,1)),sol{i}.mu,ss{i}.X(8079));
end
save YNB
