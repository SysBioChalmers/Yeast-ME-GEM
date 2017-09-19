function addTIM23Constraint(fptr,model,c1,c2,c3)
rxnID=sprintf('TIM23_Complex_biosynthesis');
syn=find(ismember(model.rxns,rxnID));

rxnID=sprintf('TIM23_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('TIM23_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

[ve ve_matrix] = importIMS(model);
fprintf(fptr,'CTIM23_syn: %s - %.15f X%d = 0\n',ve_matrix,c1,syn );
fprintf(fptr,'CTIM23_dil: X%d - %.15f X%d = 0\n',syn_dil,c2,syn );
fprintf(fptr,'CTIM23_deg: X%d - %.15f X%d = 0\n',syn_deg,c3,syn );
