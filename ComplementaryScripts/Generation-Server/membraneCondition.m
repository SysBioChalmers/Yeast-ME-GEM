function constrain = membraneConditionAngelica(model, mu,kdeg)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');
[a b complex_list]=xlsread('TableS1.xlsx','Mito_Intermembrane_proteins');
[a b subunit]=xlsread('TableS1.xlsx','Subunit');
 
constrain='';

complex_list(numel(complex_list)+1)={'TIM23_Complex_biosynthesis'};
complex_list(numel(complex_list)+1)={'TIM22_Complex_biosynthesis'};

for i=1:numel(complex_list) 
   
    
    c=cell2mat(strrep(complex_list(i),':','_'));
    I=find(ismember({'YDL119C','YDL022W'},c));
    if numel(I)==0
        rxnID=sprintf('%s_complex_formation_mitochondrion',c);
    else
       rxnID=sprintf('%s_complex_formation',c);

    end
    I=find(ismember({'YDL022W'},c));
    if numel(I)==1
          rxnID=sprintf('%s_complex_formation_cytoplasm',c);
    end
    rxnID=strrep(rxnID,'-','');
    s=find(ismember(model.rxns,rxnID));
    if numel(s)==0
        s=find(ismember(model.rxns,cell2mat(complex_list(i))));
        if numel(s)==0
          error(sprintf('%s reaction is not in the model',rxnID));
        end
    end
     if mod(i,300)==0
        sep=char(10);
    else
        sep='';
     end
    
     %%
     MW=0;
     if strcmp(complex_list{i},'TIM23_Complex_biosynthesis')
         MW=406270.3;
     elseif strcmp(complex_list{i},'TIM22_Complex_biosynthesis')
         MW=110299.1;
     else
         peptide=regexp(complex_list{i},':','split');
         
         %complex is dimer
         if numel(peptide)==1
             ns=1;
             I=find(ismember(subunit(:,1),peptide));
             if numel(I)==1
                 nsubunit=strrep(cell2mat(subunit(I,2)),'A','');
                 nsubunit=sprintf('%s ',nsubunit);
                 ns=str2num(nsubunit);
             end
         elseif (numel(peptide)>1)  %complex is a hetrodimer.
             I=find(ismember(subunit(:,1),complex_list{i}));
             ns=ones(numel(peptide),1);
             if numel(I)==1
                 subunits=regexp(cell2mat(subunit(I,2)),',','split');
                 for jj=1:numel(subunits)
                     ns(jj)=str2num(subunits{jj});
                 end
             end
         end
         for m1=1:numel(peptide)
             peptide_index=find(ismember(Proteins_MW(:,1),strrep(peptide(m1),'-','')));
             MW = MW + ns(m1) * cell2mat(Proteins_MW(peptide_index,2));
         end
     end
     %%
    area= peptideArea(MW);
    
    one_copy_per_cell=(1/(13*6.023e8));
    number = 13*6.023e8;
    total_area = number* area/(mu+kdeg);
    
    if strcmp(constrain,'')
        constrain = sprintf('%.15f X%d',total_area,s);
    else
        constrain = sprintf('%s + %.15f X%d%s',constrain, total_area,s,sep);
    end
    
end

