function addTIM22Constraint(fptr,model,c1,c2,c3)
rxnID=sprintf('TIM22_Complex_biosynthesis');
syn=find(ismember(model.rxns,rxnID));

rxnID=sprintf('TIM22_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('TIM22_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

[ve_mim ve_matrix] = importIMS(model);
fprintf(fptr,'CTIM22_syn: %s - %.15f X%d = 0\n',ve_mim,c1,syn );
fprintf(fptr,'CTIM22_dil: X%d - %.15f X%d = 0\n',syn_dil,c2,syn );
fprintf(fptr,'CTIM22_deg: X%d - %.15f X%d = 0\n',syn_deg,c3,syn );