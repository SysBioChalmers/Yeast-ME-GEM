function [constrain1 constrain2] = cytosolCondition(model, mu,kdeg)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');

[a b complex_list]=xlsread('TableS1.xlsx','PlasmaMembrane');

 constrain1='';
 constrain2='';

for i=1:numel(complex_list) %2
    c=cell2mat(strrep(complex_list(i),':','_'));
    rxnID=sprintf('%s_complex_formation',c);
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
    
    area= peptideArea(MW);
    
    one_copy_per_cell=(1/(13*6.6023e8))*(mu+kdeg);
    number = 1 / one_copy_per_cell;   %1e9 for convert nm^3 to mu m^3
    total_volume = number* area;
    if i<=50
        if strcmp(constrain1,'')
            constrain1 = sprintf('%.15f X%d',total_volume,s);
        else
            constrain1 = sprintf('%s + %.15f X%d',constrain1, total_volume,s);
        end
    else
        if strcmp(constrain2,'')
            constrain2 = sprintf('%.15f X%d',total_volume,s);
        else
            constrain2 = sprintf('%s + %.15f X%d',constrain2, total_volume,s);
        end
    end
    
end

