function [nATP nGTP]= codon2tRNA(codon)
[a b ctRNA]=xlsread('tRNA_modification_position_1.xlsx','charged');
codon=ctRNA(2:65,1);
uncharged=ctRNA(2:65,2);
charged=ctRNA(2:65,3);
