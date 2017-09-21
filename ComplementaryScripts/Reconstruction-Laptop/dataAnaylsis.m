% length of 5'UTR region
[a b Proteins_length]=xlsread('TableS1.xlsx','Length');
n=numel(Proteins_length(:,1));
L=cell2mat(Proteins_length(2:n,3));
hist(L,[1:10:1200]);