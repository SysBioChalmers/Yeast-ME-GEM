function [syn dil deg nsubunits]=printeETF(model,genes)


rs='';
syn=zeros(numel(genes),1);
for gi=1:numel(genes)
    rxnID=sprintf('%s_complex_formation',cell2mat(genes(gi)));
    ss=find(ismember(model.rxns,rxnID));
    syn(gi)=ss(1);
    
end

nsubunits = length(find(full(model.S(:,syn(1)))<0));


  rxnID=sprintf('%s_dilution',cell2mat(genes(gi)));
  dil=find(ismember(model.rxns,rxnID));
  if numel(find(ismember({'Intialtion_Factors','eEF1B', 'Relase_Factors'},genes)))>=1
      rxnID=sprintf('%s_subunit_degradation',cell2mat(genes(gi)));
  else
      rxnID=sprintf('%s_degradation',cell2mat(genes(gi)));
  end
  deg=find(ismember(model.rxns,rxnID));

