c=2;
[n a ctRNA]=xlsread('tRNA_modification_position_1.xlsx','charged');
for i=2:numel(ctRNA(:,c))
    met=sprintf('%s [cytoplasm]',cell2mat(ctRNA(i,c)));
    sol=testMetabolite(newModel2,met);
   fprintf('%s %f\n',met,sol.f);
end