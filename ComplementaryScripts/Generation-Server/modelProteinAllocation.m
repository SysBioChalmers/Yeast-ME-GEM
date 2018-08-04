function [geneList PredictedExpression ribosomeRatio] = modelProteinAllocation(model,solFile,mu,kdeg,SoplexVersion)
if strcmp(SoplexVersion,'Soplex-3.0')==1
  sol=readSoplex3_results(solFile,model);
else
sol=readSoplex_results(solFile,model);
end
flux=sol.X;
[a b Proteins]=xlsread('TableS1.xlsx','gene_seq');
[a b Proteins_mw]=xlsread('TableS1.xlsx','MW');
[a b Proteins_annotation]=xlsread('annotation.xlsx');

% get amino acids
geneName={''};
Expression=[];
MW=[];
k=1;
for r=1:numel(model.rxns)
    %fprintf('%d\n',r)
    name=regexp(cell2mat(model.rxns(r)),'_translation');
    if numel(name)>=1
        gene=regexp(cell2mat(model.rxns(r)),'\w*_translation','match');
        gene=strrep(gene,'_translation','');
        geneName(k)=gene;
        Expression(k)=flux(r)/(mu+kdeg); %unit is mmol/gDW
        
        indexMW=find(ismember(Proteins_mw(:,1),gene));
        MW(k)= cell2mat(Proteins_mw(indexMW,2));
        
        Expression(k)=flux(r)/(mu+kdeg)*MW(k)/1000; %unit is g/gDW

        
        k=k+1;
    end
end

geneList=unique(geneName);
PredictedExpression=[];
molecule_number=floor(proteinRatio(mu)*13/(300*128)*6.023e23/1e12); 

for i=1:numel(geneList)
    index=find(ismember(geneName,geneList(i)));
    PredictedExpression(i)=sum(Expression(index))/proteinRatio(mu); %the unit is g/g prpotein;*13*6.023e8/molecule_number; % the unit is g per gDW %%*13*6.023e8; %%
end

%print the data
fileName=strrep(solFile,'.out','.list');
fptr=fopen(fileName,'w');
for i=1:numel(geneList)
    fprintf(fptr,'%s\t%.15f\n', cell2mat(geneList(i)), PredictedExpression(i));
end
fclose(fptr);
fprintf('The data was stored in the file %s\n',fileName);

% %Estimate Metabolic enzyme weights
fileName=strrep(solFile,'.out','.allocation');
fptr=fopen(fileName,'w');
total_metabolic =0;
total_ribosome =0;
total_mitoribosome =0;
total_proteasome=0;
total_unmodel =0;
total_TranslationFactors=0;
total_Chaperone=0;
total_Ribosome_Assembly_Factors=0;
total_amino_acids=0;
total_carbon=0;
total_lipds=0;
total_NTs=0;
total_TOM=0;
for i=1:numel(geneList)
    index=find(ismember(Proteins_annotation(:,2),geneList(i)));
    if numel(index)==0
      total_unmodel =total_unmodel+PredictedExpression(i);

    end

    if numel(index)>=1
        class2='';
        class=cell2mat(Proteins_annotation(index(1),1));
        if numel(index)>1
            class2=cell2mat(Proteins_annotation(index(2),1));
        end
        f=0;
        if strcmp(class, 'Metabolic')==1
            total_metabolic=total_metabolic+ PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'Ribosome')==1
            total_ribosome=total_ribosome+PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'MitoRibosome')==1
            total_mitoribosome=total_mitoribosome+PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'TranslationFactors')==1
            total_TranslationFactors =total_TranslationFactors+PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'Chaperone')==1
            total_Chaperone =total_Chaperone+PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'Ribosome Assembly Factors')==1
            total_Ribosome_Assembly_Factors =total_Ribosome_Assembly_Factors+PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'TOM')==1
            total_TOM =total_TOM+PredictedExpression(i);
            f=1;
        end
        if strcmp(class, 'proteasome')==1
            total_proteasome =total_proteasome+PredictedExpression(i);
            f=1;
        end
        if strcmp(class2, 'C')==1
            total_carbon=total_carbon+PredictedExpression(i);
            f=1;
        end
        if strcmp(class2, 'AA')==1
            total_amino_acids=total_amino_acids+PredictedExpression(i);
            f=1;
        end
        if strcmp(class2, 'Lipid')==1
            total_lipds=total_lipds+PredictedExpression(i);
            f=1;
        end
        if strcmp(class2, 'NTs')==1
            total_NTs=total_NTs+PredictedExpression(i);
            f=1;
        end

        if f==0
            total_unmodel =total_unmodel+PredictedExpression(i);
        end
    end
end
fprintf('The allocation data was stored in the file %s\n',fileName);
fprintf(fptr,'total protein molecules\t%.15f\n',sum(PredictedExpression));
fprintf(fptr,'metabolic\t%.15f\n',total_metabolic/sum(PredictedExpression));
fprintf(fptr,'Ribosome\t%.15f\n',total_ribosome/sum(PredictedExpression));
fprintf(fptr,'MitoRibosome\t%.15f\n',total_mitoribosome/sum(PredictedExpression));
fprintf(fptr,'Translation factors\t%.15f\n',total_TranslationFactors/sum(PredictedExpression));
fprintf(fptr,'Chaperone\t%.15f\n',total_Chaperone/sum(PredictedExpression));
fprintf(fptr,'Ribosome Assembly Factors\t%.15f\n',total_Ribosome_Assembly_Factors/sum(PredictedExpression));
fprintf(fptr,'Proteasome\t%.15f\n',total_proteasome/sum(PredictedExpression));
fprintf(fptr,'TOM\t%.15f\n',total_TOM/sum(PredictedExpression));
fprintf(fptr,'Other\t%.15f\n',total_unmodel/sum(PredictedExpression));
fprintf(fptr,'Metabolism allocation');
fprintf(fptr,'Carbon\t%.15f\n',total_carbon/sum(PredictedExpression));
fprintf(fptr,'Amino Acids\t%.15f\n',total_amino_acids/sum(PredictedExpression));
fprintf(fptr,'Lipids\t%.15f\n',total_lipds/sum(PredictedExpression));
fprintf(fptr,'NTs\t%.15f\n',total_NTs/sum(PredictedExpression));

fclose(fptr)
ribosomeRatio= total_TranslationFactors/sum(PredictedExpression);
