function [count aa_names] = aaConentration(con,seq)
aa=aacount(seq);
aa_names=fieldnames(aa);
count = zeros(20,1);
for i=1:20
    count(i)=getfield(aa,cell2mat(aa_names(i)));
    % * concentration
    count(i)=count(i)*con;
end