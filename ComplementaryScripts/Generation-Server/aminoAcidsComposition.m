function [eq_cytosol eq_mito count_cytosplasm count_mito]=aminoAcidsComposition(model,outputFile,mu)

%reading the solution file
sol=readSoplex3_results(outputFile,model);


[a b AAs]=xlsread('aminoAcidsNames.xlsx');

count_cytosplasm=zeros(20,1);
count_mito=zeros(20,1);


for i=1:numel(model.rxns)
    p=regexp(cell2mat(model.rxns(i)),'loading');
    if numel(p)>0
         p=regexp(cell2mat(model.rxns(i)),'_loading_mitochondrion');
         if numel(p)==0
             %amino acids in amino acids
             rxns=cell2mat(model.rxns(i));
             aa=rxns(1:3);
             
             index=find(ismember(AAs(:,1),aa));
             
             %fprintf('%s %f\n',aa,sol.X(i)/mu);

             count_cytosplasm(index)=count_cytosplasm(index) + sol.X(i)/mu;
             
             
         else
             rxns=cell2mat(model.rxns(i));
             aa=rxns(1:3);
             
             index=find(ismember(AAs(:,1),aa));
             count_mito(index)=count_mito(index) + sol.X(i)/mu;
  
             %fprintf('Mito: %s %f\n',aa,sol.X(i)/mu);
             
         end
    end
end

%estimating GTP in Mitochondira proteins

flux_sum_mito_protein=0;
for i=1:numel(model.rxns)
    p=regexp(cell2mat(model.rxns(i)),'_translation_mitochondrion');
    if numel(p)>0
        flux_sum_mito_protein=flux_sum_mito_protein + sol.X(i);
    end
end

nGTP=flux_sum_mito_protein*4440/mu;


eq_cytosol='';
eq_mito='';

for i=1:20
 if i==1
  eq_cytosol=sprintf('%f %s [cytoplasm]',count_cytosplasm(i), cell2mat(AAs(i,2)));

  eq_mito=sprintf('%f %s [mitochondrion]',count_mito(i), cell2mat(AAs(i,2)));

 else
    eq_cytosol=sprintf('%s + %f %s [cytoplasm]',eq_cytosol, count_cytosplasm(i), cell2mat(AAs(i,2)));

    eq_mito=sprintf('%s + %f %s [mitochondrion]',eq_mito, count_mito(i), cell2mat(AAs(i,2)));

 end
end
eq_cytosol=sprintf('%s => protein [cytoplasm]',eq_cytosol );
eq_mito=sprintf('%s + %f GTP [mitochondrion] + %f H2O [mitochondrion] => protein [mitochondrion] + %f H+ [mitochondrion] + %f phosphate [mitochondrion] + %f GDP [mitochondrion]',eq_mito,nGTP,nGTP,nGTP,nGTP,nGTP);
