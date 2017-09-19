c=2;
[n a proteins]=xlsread('TableS1.xlsx','Annotation');
for i=10:1100
    met=sprintf('%s_peptide [cytoplasm]',cell2mat(proteins(i,c)));
    sol=testMetabolite(newModel2,met);
    if sol.f==0
        fprintf('%s %f\n',met,sol.f);
    end
end