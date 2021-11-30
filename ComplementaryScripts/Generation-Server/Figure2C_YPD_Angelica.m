%This file generates glucose-limit growth
%clear all
load NewYeastnextGeneration.mat
%%Fructose transporter
eq=sprintf('D-galactose [extracellular] => D-galactose [cytoplasm]');
rxnID=sprintf('fructose_diffusion');
rxnNames=rxnID;
MyModel=addYeastReaction(MyModel,eq,{rxnID},{rxnNames},0,3,0,{'YLR081W'}); % PMID: 7852299
    

[ carbone b c]=xlsread('other_carbon.xlsx','Data_minimum');
[ carbone b c]=xlsread('other_carbon.xlsx');

sol={};
%ss={};
NN=numel(c(:,1));
MM=1;%
%NN=4;
 atpCost.gam=30;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
 atpCost.ngam=1;
  T=33;

 overExpressed.syn='mCherry_complex';
 overExpressed.dilultion='mCherry_complex_dilution';
 overExpressed.degradation='mCherry_degradation';

 
copies =0;
model=MyModel;
% %adjust media
% [a b metabolite]=xlsread('other_carbon.xlsx','Medium');
% model=MyModel;
% for i=2:numel(metabolite(:,1))
%     eq=sprintf(' => %s',cell2mat(metabolite(i,2)));
%     rxnID=sprintf('%s_import',cell2mat(metabolite(i,1)));
%     rxnNames=rxnID;
%     model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,cell2mat(metabolite(i,3)),0,{''}); 
% end

   model.grRules=strrep(model.grRules,'Proteasome','');
   model.grRules=strrep(model.grRules,'Ribosome','');
   model.grRules=strrep(model.grRules,'TOM-TIM','');
   model.grRules=strrep(model.grRules,'MitoDegradome','');
uptakes=[2];
c=c([12],:);
PredictedProteinMass=zeros(836,1);

for i=1:1 % add carbon transporter
 
eq=sprintf(' => %s',cell2mat(c(i,3)));
rxnID=sprintf('%s_transporter',cell2mat(c(i,1)));
rxnNames=sprintf('transporter');

new_MyModel=addYeastReaction(model,eq,{rxnID},{rxnNames},0,cell2mat(c(i,4)),-1,{''});
new_MyModel=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,-1,{''});

fileName=sprintf('%s_batch_NDI_uptake.lp',cell2mat(c(i,1)));



% fileAbundance.fileName='growth_rate_0.2.xlsx';
% fileAbundance.sheet='growth_rate_0.2';


%ExpressionFileName=sprintf('growth_rate_%.2g.xlsx',mu);
[sol{i} glc  PredictedProteinMass(:,i) newModel]=growthSearchDynamic(new_MyModel, 'minimum',atpCost,0,T, copies, overExpressed, complex_list,compartment,0.0,0.5,fileName,'Soplex-1.6',1,'')
validateSolution(sol{1}.fileName,newModel)

end

fptr_metabolic_proteome = fopen('metabolic_proteome_Angelica_YPD.txt','w');
for j=1:836 %number of complex
  fprintf(fptr_metabolic_proteome,'%.15f\n',PredictedProteinMass(j,1));


end



% system('echo "Intialtion_Factors_dilution"; grep "X12782" *.out > Intialtion_Factors_dilution_result.txt','-echo');
for i=1
    fprintf('%s\t%f\t%f\t%f\t%f\t%.15g\t%.15g\t%.15g\n',cell2mat(c(i,1)),sol{i}.mu,sol{i}.mu,sol{i}.X(4321),sol{i}.X(4321),sol{i}.X(8079),sol{i}.X(18085),sol{i}.X(18085));
end
save Batch_min_angelica_uptake