function [constrain1 constrain2 constrain3] = cytosolCondition(model, mu,kdeg)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');

[a b complex_list]=xlsread('TableS1.xlsx','Cytosolic_complex');

 constrain1='';
 constrain2='';
  constrain3='';

for i=1:numel(complex_list) %2
  
    c=cell2mat(strrep(complex_list(i),':','_'));
    rxnID=sprintf('%s_translation',c);
    rxnID=strrep(rxnID,'-','');
    s=find(ismember(model.rxns,rxnID));
    
    
    if mod(i,50)==0
        sep=char(10);
    else
        sep='';
    end
    
    MW=0;
    peptide=regexp(complex_list{i},':','split');
    for m1=1:numel(peptide)
        peptide_index=find(ismember(Proteins_MW(:,1),peptide(m1)));
        MW = MW + cell2mat(Proteins_MW(peptide_index,2));
    end
    
    area= peptideVolume(MW);
    
    one_copy_per_cell=(1/(13*6.6023e8));
    number = 1 / one_copy_per_cell;   %1e9 for convert nm^3 to mu m^3
    total_volume = number*area/(mu+kdeg);
    if i<=350
        if strcmp(constrain1,'')
            constrain1 = sprintf('%f X%d',total_volume,s);
        else
            for j=1:numel(s)
            constrain1 = sprintf('%s + %f X%d',constrain1, total_volume,s(j));
            end
        end
    elseif i<=600
        if strcmp(constrain2,'')
            constrain2 = sprintf('%f X%d',total_volume,s);
        else
            for j=1:numel(s)
                constrain2 = sprintf('%s + %f X%d',constrain2, total_volume,s(j));
            end
        end
    else
        if strcmp(constrain3,'')
            constrain3 = sprintf('%f X%d',total_volume,s);
        else
            for j=1:numel(s)
                constrain3 = sprintf('%s + %f X%d',constrain3, total_volume,s(j));
            end
        end
    end
    
end

