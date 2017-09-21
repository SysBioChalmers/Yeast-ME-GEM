function constrain = GlycolysisCondition(model, mu,kdeg)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');

[a b complex_list1]=xlsread('TableS1.xlsx','kcat');
I=find(ismember(complex_list1(:,3),'G'));
complex_list=complex_list1(I,1);
 constrain='';

for i=1:numel(complex_list) %2
    c=cell2mat(strrep(complex_list(i),':','_'));
    c=strrep(c,'-','');

    rxnID=sprintf('%s_complex_formation_cytoplasm',c);
    
    s=find(ismember(model.rxns,rxnID));
    if numel(s)==0
         rxnID=sprintf('%s_complex_formation_mitochondrion',c);
         s=find(ismember(model.rxns,rxnID));   
    end
     if mod(i,300)==0
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

    constant=MW/1000.0/(mu+kdeg);
  
    
    if strcmp(constrain,'')
        constrain = sprintf('%.15f X%d',constant,s);
    else
        constrain = sprintf('%s + %.15f X%d%s',constrain, constant,s,sep);
    end
    
end

