function printeETF(model,fptr,genes, c1,c2,c3,ve)


rs='';
syn=zeros(numel(genes),1);
for gi=1:numel(genes)
    rxnID=sprintf('%s_complex_formation',cell2mat(genes(gi)));
    ss=find(ismember(model.rxns,rxnID));
    syn(gi)=ss(1);
    if strcmp(rs,'')
        rs=sprintf('- %.15f X%d',c1,syn(gi));
    else
        rs=sprintf('%s - %.15f X%d',rs,c1, syn(gi));
    end
end

fprintf(fptr,'CeEF31%d: %s %s <= 0\n',1,ve,rs);

for gi=1:numel(genes)
  rxnID=sprintf('%s_dilution',cell2mat(genes(gi)));
  dil=find(ismember(model.rxns,rxnID));
  fprintf(fptr,'CeEF32: X%d - %.15f X%d = 0\n',dil,c2,syn(gi) );

  rxnID=sprintf('%s_subunit_degradation',cell2mat(genes(gi)));
  deg=find(ismember(model.rxns,rxnID));
  fprintf(fptr,'CeEF32: X%d - %.15f X%d = 0\n',deg(1),c3,syn(gi) );
end
