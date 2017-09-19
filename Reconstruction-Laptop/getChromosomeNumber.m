function n=getChromosomeNumber(chr)
%this function convert the chromosome ID to number or index
if strcmp(chr,'chrI')
    n=1;
elseif strcmp(chr,'chrII')
    n=2;
elseif strcmp(chr,'chrIII')
    n=3;
elseif strcmp(chr,'chrIV')
    n=4;
elseif strcmp(chr,'chrV')
    n=5;
elseif strcmp(chr,'chrVI')
    n=6;
elseif strcmp(chr,'chrVII')
    n=7;
elseif strcmp(chr,'chrVIII')
    n=8;
elseif strcmp(chr,'chrIX')
    n=9;
elseif strcmp(chr,'chrX')
    n=10;
elseif strcmp(chr,'chrXI')
    n=11;
elseif strcmp(chr,'chrXII')
    n=12;
elseif strcmp(chr,'chrXIII')
    n=13;
elseif strcmp(chr,'chrXIV')
    n=14;
elseif strcmp(chr,'chrXV')
    n=15;
elseif strcmp(chr,'chrXVI')
    n=16;
elseif strcmp(chr,'chrMito')
    n=17;
elseif strcmp(chr,'chrM')
    n=17;
elseif strcmp(chr,'2-micron')
    n=17;
end
    
    
   