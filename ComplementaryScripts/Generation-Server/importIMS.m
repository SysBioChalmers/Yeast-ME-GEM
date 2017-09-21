function [ve_mim ve_matrix] = importIMS(model)
[a b proteins]=xlsread('TableS1.xlsx','Mito_Intermembrane_proteins');
 MIM_proteins={''};
k=1;
for i=1:numel(proteins)
    peptide=regexp(cell2mat(proteins(i)),':','split');
    for j=1:numel(peptide);
        MIM_proteins(k)=peptide(j);
        k=k+1;
    end
end

k=0;
ve_mim='';
ve_matrix='';

for i=1:numel(model.rxns)
    name=regexp(cell2mat(model.rxns(i)),'importing_mitochondrion_Matrix');
    if numel(name)>0
        gene=strrep(cell2mat(model.rxns(i)),'_importing_mitochondrion_Matrix','');
        if numel(find(ismember(MIM_proteins,gene)))>0
            %if a gene is in MIM, it is added to MIM
            if strcmp(ve_mim,'')
                ve_mim = sprintf('X%d',i);
            else
                ve_mim=sprintf('%s + X%d%c',ve_mim,i,sep);
            end
        else
            if strcmp(ve_matrix,'')
                ve_matrix = sprintf('X%d',i);
            else
                ve_matrix = sprintf('%s + X%d%c',ve_matrix,i,sep);
            end
        end
        
        k=k+1;
        if mod(k,300)==0
            sep=char(10);
        else
            sep='';
        end
        
    end
    
end
