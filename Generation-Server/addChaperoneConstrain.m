function addChaperoneConstrain(fptr,model,c1,c2,c3,chaperone)
[a b chap]=xlsread('TableS1.xlsx','Folding');
J=find(ismember(chap(1,:),chaperone));
substrate=chap(:,J);
i=0;
rs='';
for i=1:numel(substrate)
    rxnID=sprintf('%s_folding_cytoplasm',cell2mat(substrate(i)));
    ss=find(ismember(model.rxns,rxnID));
    if numel(ss)==0
         rxnID=sprintf('%s_folding',cell2mat(substrate(i)));
         ss=find(ismember(model.rxns,rxnID));
    end
    if numel(ss)>0
    syn=ss(1);
    
     if strcmp(rs,'')
        rs=sprintf('X%d',syn);
    else
        rs=sprintf('%s + X%d',rs,syn);
    end
    end
end

rxnID=sprintf('%s_biosynthesis',chaperone);
ss=find(ismember(model.rxns,rxnID));
syn=ss(1);
fprintf(fptr,'%sSYN: %s - %.15f X%d <= 0\n',chaperone,rs,c1,syn );

rxnID=sprintf('%s_dilution',chaperone);
ss=find(ismember(model.rxns,rxnID));
dil=ss(1);
fprintf(fptr,'%sDIL: X%d - %.15f X%d = 0\n',chaperone,dil,c2,syn );

rxnID=sprintf('%s_degradation',chaperone);
ss=find(ismember(model.rxns,rxnID));
deg=ss(1);
fprintf(fptr,'%sDEG: X%d - %.15f X%d = 0\n',chaperone,deg,c3,syn );



  