fptr_metabolic_proteome = fopen('metabolic_proteome.txt','w');
for j=1:836 %number of complex
  fprintf('%f\t',PredictedProteinMass(j,1));

for k=2:length(mu1)
fprintf('%.15\t',PredictedProteinMass(j,k));
end
fprintf('\n');
end