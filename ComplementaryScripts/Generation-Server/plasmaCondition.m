function constrain = plasmaCondition(model, mu,kdeg)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');

[a b complex_list]=xlsread('TableS1.xlsx','PlasmaMembrane');

 constrain='';

for i=1:numel(complex_list) %2
    c=cell2mat(strrep(complex_list(i),':','_'));
    rxnID=sprintf('%s_complex_formation_extracellular',c);
  
    rxnID=strrep(rxnID,'-','');
    s=find(ismember(model.rxns,rxnID));

    if numel(s)==0
         rxnID=sprintf('%s_complex_formation_cytoplasm',c);
         rxnID=strrep(rxnID,'-','');
         s=find(ismember(model.rxns,rxnID));
    end
    
    if numel(s)==0
        error(sprintf('%s is not in plasma membrane',c));
    end
     if mod(i,300)==0
        sep=char(10);
    else
        sep='';
     end
    
    MW=0;
    peptide=regexp(complex_list{i},':','split');
    for m1=1:numel(peptide)
        peptide_index=find(ismember(Proteins_MW(:,1),strrep(peptide(m1),'-','')));
        MW = MW + cell2mat(Proteins_MW(peptide_index,2));
    end
    
    area= peptideArea(MW);
    
    one_copy_per_cell=(1/(13*6.023e8))*(mu+kdeg);
    number = 1 / one_copy_per_cell ;
    total_area = number* area;
    
    if strcmp(constrain,'')
        constrain = sprintf('%.15f X%d',total_area,s);
    else
        constrain = sprintf('%s + %.15f X%d%s',constrain, total_area,s,sep);
    end
    
end

