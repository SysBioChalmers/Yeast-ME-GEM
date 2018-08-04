function [vol vol1, molecules molecules1 variable count]= mitoVolume(model, fileName,mu,kdeg)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');

[a b complex_list]=xlsread('TableS1.xlsx','mitoProteins');
sol=readSoplex3_results(fileName,model);
%complex_list(778)= {'YDRUMD1'};

vol=0;
vol1=0;
molecules=0;
molecules1=0;

k=1;
variable={''};
count=0;
for i=1:numel(complex_list) %2
     flag=0;
    c=cell2mat(strrep(complex_list(i),':','_'));
    rxnID=sprintf('%s_translation',c);
    rxnID=strrep(rxnID,'-','');
    s=find(ismember(model.rxns,rxnID));
    if  numel(s)==0
        rxnID=sprintf('%s_translation_mitochondrion',c);
        rxnID=strrep(rxnID,'-','');
        s=find(ismember(model.rxns,rxnID));
        flag=1;
    end
    
    if mod(i,50)==0
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
    
    V= peptideVolume(MW);
    
    one_copy_per_cell=(1/(13*6.023e8));
    number = 1 / one_copy_per_cell;   %1e9 for convert nm^3 to mu m^3
    total_number = number;
    %total_volume1 = MW/1000;
    total_volume1 = number * V;
    for j=1:numel(s)
        vol  = vol + total_volume1*sol.X(s(j))/(mu+kdeg);
        molecules= molecules + total_number*sol.X(s(j))/(mu+kdeg);
        %fprintf('%s:%d %.15f\n',c,i,total_volume1 *sol.X(s(j)));
        variable{k,1}=sprintf('X%d',s(j));
        count(k,1)=total_volume1*sol.X(s(j))/(mu+kdeg);
        k=k+1;
        if flag==1
                    vol1  = vol1 + total_volume1 *sol.X(s(j));
                    molecules1 = molecules1 + total_number *sol.X(s(j));
        end

    end
    
end

        fprintf('Total mito %.15f\n',vol);
        fprintf('without Membrane %.15f\n',vol-vol1);
        fprintf('Membrane %.15f\n',vol1);
