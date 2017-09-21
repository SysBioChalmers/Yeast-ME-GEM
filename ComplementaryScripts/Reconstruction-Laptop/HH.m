clear all
HH1 = readCbModel('yeast_7.6_cobra.xml');

HH =HH1;
HH.S(:,11:numel(HH.rxns))=[];
HH.rev(11:numel(HH.rxns))=[];
HH.c(11:numel(HH.rxns))=[];

HH.lb(11:numel(HH.rxns))=[];
HH.ub(11:numel(HH.rxns))=[];
HH.rules(11:numel(HH.rxns))=[];
HH.grRules(11:numel(HH.rxns))=[];
HH.subSystems(11:numel(HH.rxns))=[];
HH.confidenceScores(11:numel(HH.rxns))=[];
HH.rxnReferences(11:numel(HH.rxns))=[];
HH.rxnECNumbers(11:numel(HH.rxns))=[];
HH.rxnNotes=[];
HH.rxnNames(11:numel(HH.rxns))=[];
HH.rxns(11:numel(HH.rxns))=[];

HH2=HH;
HH=HH2;
newmets={''};
k=1;
index=0;
for i=1:numel(HH.mets)

    Index=find(full(HH.S(i,:)==0));
     if numel(Index)==10
         index(k)=i;
         k=k+1;
     end
end


       HH.mets(index)=[];
         HH.metNames(index)=[];
         HH.metFormulas(index)=[];
         HH.metCharge(index)=[];
         HH.metChEBIID(index)=[];
         HH.metKEGGID(index)=[];
         HH.metPubChemID(index)=[];
         HH.metInChIString(index)=[];
         HH.b(index)=[];
         HH.S(index,:)=[];
         
         HH3=HH;
         HH=HH3;
         k=1; 
         index=0;
for i=1:909
    i
      Index=find(full(full(HH.rxnGeneMat(:,i))==0));
     if numel(Index)==10
         index(k)=i;
         k=k+1;
     end
end
HH.genes(index)=[];
HH.rxnGeneMat(:,index)=[];
HH.rxnGeneMat(index)=[];