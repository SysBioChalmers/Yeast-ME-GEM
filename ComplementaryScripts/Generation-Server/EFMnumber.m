function [nullModel struct model efms]  = EFMnumber(model,fileName,mu)
% how many are EFMs in a LP solution

% 1- read Soplex solution 


sol=readSoplex3_results(fileName,model);
metabolic_flux=sol.X;
Flux_index = find(metabolic_flux > 1.0e-15);

metabolic_Flux=Flux_index(Flux_index<=5738); %only metabolic reactions

% please note -2 for removing maintaince and GAM reactions
%metabolic_Flux=union(metabolic_Flux,24419:numel(model.rxns)); % -2
Flux_index(find(Flux_index==max(Flux_index)))=[]; %remove this solution. It is an internal variables during solutions
metabolic_Flux=Flux_index;
reactions=model.rxns(metabolic_Flux);
S=model.S(:,metabolic_Flux);
%mets=union(1:2220,12726); %for metabolic
mets=1:numel(model.mets);
S=S(mets,:);

% newModel.rxns=model.rxns(metabolic_Flux);
% newModel.rxnNames=model.rxnNames(metabolic_Flux)
% newModel.mets=model.mets(mets);
% newModel.metNames=model.metNames(mets);
% newModel.S=S;
% newModel.metComps=model.metComps(mets);

nullModel=initModel_yeast;
newModel=model;
metabolic_rxns=[1:5738 24419:numel(model.rxns) 5756 24293 24294 24295 24287]; %24287 for GDP transporter %24293 24294 X24295 for amino acids transporter

for i=metabolic_rxns
  %if strcmp(cell2mat(newModel.rxns(i)),'test_biomass')==0
    if(sol.X(i)>=1e-15 )
    p_i=find(newModel.S(:,i)>0);
    product='';
    for j=1:numel(p_i)
        if strcmp(product,'')==1
            if newModel.S(p_i(j),i)==1
            product=sprintf('%s[%s]', cell2mat(newModel.metNames(p_i(j))),cell2mat(newModel.metComps(p_i(j))));
            else
                            product=sprintf('%f %s[%s]', full(newModel.S(p_i(j),i)) ,cell2mat(newModel.metNames(p_i(j))),cell2mat(newModel.metComps(p_i(j))));
            end

        else
            if newModel.S(p_i(j),i)==1
               product=sprintf('%s + %s[%s]', product, cell2mat(newModel.metNames(p_i(j))),cell2mat(newModel.metComps(p_i(j))));
            else
               product=sprintf('%s + %f %s[%s]', product, full(newModel.S(p_i(j),i)), cell2mat(newModel.metNames(p_i(j))),cell2mat(newModel.metComps(p_i(j))));
           end
            
        end
    end
    
     s_i=find(newModel.S(:,i)<0);
    substrate='';
    for j=1:numel(s_i)
        if strcmp(substrate,'')==1
            if newModel.S(s_i(j),i)==-1
            substrate=sprintf('%s[%s]', cell2mat(newModel.metNames(s_i(j))),cell2mat(newModel.metComps(s_i(j))));
            else
              substrate=sprintf('%f %s[%s]', - full(newModel.S(s_i(j),i)), cell2mat(newModel.metNames(s_i(j))),cell2mat(newModel.metComps(s_i(j))));

            end
        else
            if newModel.S(s_i(j),i)==-1
            substrate=sprintf('%s + %s[%s]', substrate, cell2mat(newModel.metNames(s_i(j))),cell2mat(newModel.metComps(s_i(j))));
            else
            substrate=sprintf('%s + %f %s[%s]', substrate, -full(newModel.S(s_i(j),i)), cell2mat(newModel.metNames(s_i(j))),cell2mat(newModel.metComps(s_i(j))));
                
            end
            
        end
    end
    
    eq=sprintf('%s => %s',substrate,product);
    fprintf('%s\n',eq);

    nullModel=addYeastReaction(nullModel,eq,newModel.rxns(i),newModel.rxnNames(i),0,1000,0,{''});
    end
    
  %end % remove exporting biomass without proteins
end
% %addding protein
 [eq_cytosol eq_mito count_cytosplasm count_mito]=aminoAcidsComposition(model,fileName,mu);
rxns={'cytoplasm_protein'};
rxnsNames={'cytoplasm_protein'};
nullModel=addYeastReaction(nullModel,eq_cytosol,rxns,rxnsNames,0,1000,0,{''});

% eq_cytosol
rxns={'mito_protein'};
rxnsNames={'mito_protein'};
nullModel=addYeastReaction(nullModel,eq_mito,rxns,rxnsNames,0,1000,0,{''});
eq_mito

eq='protein [cytoplasm] => ';
rxns={'cytoplasm_protein_export'};
rxnsNames={'cytoplasm_protein_export'};
nullModel=addYeastReaction(nullModel,eq,rxns,rxnsNames,0,1000,0,{''});
eq



 eq='protein [mitochondrion] => '
 rxns={'mito_protein_export'};
 rxnsNames={'mito_protein_export'};
 nullModel=addYeastReaction(nullModel,eq,rxns,rxnsNames,0,1000,0,{''});
 eq

% eq='protein [cytoplasm] => Total Protein [cytoplasm]'
% rxns={'total_proteome'};
% rxnsNames={'total_proteome'};
% nullModel=addYeastReaction(nullModel,eq,rxns,rxnsNames,0,1000,0,{''});
% eq
% 
% 
% eq='Total Protein [cytoplasm] + biomass_test [cytoplasm] => Total Biomass [cytoplasm]'
% rxns={'total_biomass'};
% rxnsNames={'total_biomass'};
% nullModel=addYeastReaction(nullModel,eq,rxns,rxnsNames,0,1000,0,{''});
% eq
% 
%  eq='Total Biomass [cytoplasm] => '
%  rxns={'total_biomass_export'};
%  rxnsNames={'total_biomass_export'};
%  nullModel=addYeastReaction(nullModel,eq,rxns,rxnsNames,0,1000,0,{''});
%  eq

rr=rank(full(nullModel.S));
size(rr);
n=null(full(nullModel.S));
size(n)

%finding EFMs
%loop for the metabolites
fprintf('Check S\n');
model=nullModel;
for i=1:numel(model.S(:,1))
    row=full(model.S(i,:));
    n_reactions = length(find(row~=0));
    if n_reactions==1
        fprintf('%s[%s]\n',cell2mat(model.metNames(i)),cell2mat(model.metComps(i)));
        %add this reaction to the biomass
        model.S(i,find(ismember(model.rxns,'biomass_without_protein')))=-1;
    end
end

rr=rank(full(model.S));
size(rr);
n=null(full(model.S));
size(n)

struct.stoich=full(model.S);
struct.reversibilities=zeros(numel(model.rxns),1);

cd '/zhome/88/9/107870/efmtool'
%  exchange=[395 398 408 413 420 435 438 444 451 539 547 548];
% struct.reversibilities(exchange)=1;
efms=CalculateFluxModes(struct.stoich,struct.reversibilities)
 cd '/zhome/88/9/107870/Yeast-ME-GEM/ComplementaryScripts/Generation-Server'
