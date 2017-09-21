[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');
constraint={'' '' '' '' '' ''};
k=1;
for i=1:numel(model.rxns)

    name=regexp(cell2mat(model.rxns(i)),'_translation');
    
    if numel(name)>0
       name=regexp(cell2mat(model.rxns(i)),'\w*_translation','match');
       gene=strrep(name,'_translation','');
       peptide_index=find(ismember(Proteins_MW(:,1),gene));
       kdeg= cell2mat(Proteins_MW(peptide_index,3));
       MW= cell2mat(Proteins_MW(peptide_index,3));
       
       if length(constraint{k}) >= 7500
           k=k+1;
       end
       if strcmp(constraint{k},'')==1
           constraint{k}=sprintf('%.15f X%d',MW/(mu + kdeg)/1000, i);
       else
           constraint{k}=sprintf('%s + %.15f X%d',constraint{k},MW/(mu + kdeg)/1000, i);
       end
    end
    
end

for i=1:k
    fprintf('fptr,%s -XTP%d =0\n',constraint{i},i);
end
XTP='';
for i=1:k
    if strcmp(XTP,'')
        XTP='XTP1';
    else
        XTP=sprintf('fptr,%s + XTP%d',XTP,i);
    end
end
fprintf(fptr,'%s  =  %f\n',XTP,proteinRatio(mu)*(1-upr1));

    