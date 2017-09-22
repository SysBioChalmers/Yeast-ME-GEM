function [eq_substrates eq_product nATP nGTP]= codon_features(seq)
[a b ctRNA]=xlsread('tRNA_modification_position_1.xlsx','charged');
codon=ctRNA(2:65,1);
uncharged=ctRNA(2:65,2);
charged=ctRNA(2:65,3);


% The aim of this function is to count how many codon in the RNA seq
codon=codoncount(seq);
codon_names=fieldnames(codon);
codon_count = zeros(64,1);
for i=1:64
    codon_count(i)=getfield(codon,cell2mat(codon_names(i)));
end
eq_substrates='';
eq_product='';
number_of_tRNA=0;
for i=1:64
    if (codon_count(i) >0)
        number_of_tRNA = number_of_tRNA + codon_count(i);
        if (strcmp(eq_substrates,''))
            eq_substrates= sprintf('%d %s',codon_count(i),cell2mat(charged(i)));
            eq_product= sprintf('%d %s',codon_count(i),cell2mat(uncharged(i)));
        else
            eq_substrates= sprintf('%s + %d %s',eq_substrates,codon_count(i),cell2mat(charged(i)));
            eq_product= sprintf('%s + %d %s',eq_product,codon_count(i),cell2mat(uncharged(i)));
        end
    end

end
nATP = number_of_tRNA;
nGTP = number_of_tRNA*2;