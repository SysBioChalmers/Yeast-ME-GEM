%[newModel1,reaction_list complex_list compartment substrate_list product_list  lb ub type_list ] 
[a b Subunit]=xlsread('TableS1.xlsx','Subunit');
%get reaction ID
fptr=fopen('reaction_gene.txt','w');
r={''};
for i=1:numel(reaction_list)
    s=strrep(reaction_list(i),'_forward','');
    s=strrep(s,'_reverse','');
    r(i)=s;
end

% get reaction Complexes mapping
unique_complex=unique(complex_list);
unique_reaction=unique(r);
reaction_complex=zeros(numel(unique_complex),numel(r));
for i=2:numel(unique_complex)
    g=regexp(cell2mat(unique_complex(i)),':','split');
    if numel(g)==1
        I=find(ismember(complex_list,unique_complex(i)));
        reaction=r(I);
        for j=1:numel(reaction)
            J=find(ismember(r,reaction(j)));
            reaction_complex(i,J)=1;
        end
    end
    
    
end
s=0;
for i=1:numel(unique_reaction)
    K=0;
    I=find(reaction_complex(:,i));
    genes=unique_complex(I);
    for j=1:numel(genes)
        I=find(ismember(Subunit(:,1),genes(j)));
        type=Subunit(I,2);
        if numel(type)==0
        fprintf(fptr,'%s\t%s\tNA\n',cell2mat(unique_reaction(i)),cell2mat(genes(j)));
        else
        for k=1:numel(type)
            fprintf(fptr,'%s\t%s\t%s\n',cell2mat(unique_reaction(i)),cell2mat(genes(j)),cell2mat(type(k)));
        end
        end
    end
    s=s+K;
end
fclose(fptr);