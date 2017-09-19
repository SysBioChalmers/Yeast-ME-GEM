function  [eq_product number_of_aa]= degradeMito(seq)
[a b AA]=xlsread('tRNA_modification_position_1.xlsx','AA_coding');
aa=AA(1:20,1);
aa_names=AA(1:20,2);

% The aim of this function is to count how many codon in the RNA seq
codon=aacount(nt2aa(seq));
aa_count = zeros(20,1);
for i=1:20
    aa_count(i)=getfield(codon,cell2mat(aa(i)));
end

eq_product='';
number_of_aa=0;
for i=1:20
    if (aa_count(i) >0)
        number_of_aa = number_of_aa + aa_count(i);
        if (strcmp(eq_product,''))
            eq_product= sprintf('%d %s [mitochondrion]',aa_count(i),cell2mat(aa_names(i)));
        else
            eq_product= sprintf('%s + %d %s [mitochondrion]',eq_product,aa_count(i),cell2mat(aa_names(i)));
        end
    end

end

