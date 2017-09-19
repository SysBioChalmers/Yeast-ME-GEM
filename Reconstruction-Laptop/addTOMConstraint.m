function addTOMConstraint(fptr,model,c1,c2,c3)
rxnID=sprintf('Tom_Complex_biosynthesis');
syn=find(ismember(model.rxns,rxnID));

rxnID=sprintf('TOM_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('TOM_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

[ve ve_matrix] = importIMS(model);
fprintf(fptr,'CTOM_syn: %s - %.15f X%d = 0\n',ve_matrix,c1,syn );
fprintf(fptr,'CTOM_dil: X%d - %.15f X%d = 0\n',syn_dil,c2,syn );
fprintf(fptr,'CTOM_deg: X%d - %.15f X%d = 0\n',syn_deg,c3,syn );
