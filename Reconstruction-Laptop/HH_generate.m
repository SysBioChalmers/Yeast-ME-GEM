%clear all
HH1=importModel('iTO977_v1.00_raven.xml',true);

HH =HH1;
HH.S(:,11:numel(HH.rxns))=[];
HH.rev(11:numel(HH.rxns))=[];
HH.c(11:numel(HH.rxns))=[];

HH.lb(11:numel(HH.rxns))=[];
HH.ub(11:numel(HH.rxns))=[];
HH.grRules(11:numel(HH.rxns))=[];
HH.subSystems(11:numel(HH.rxns))=[];
HH.eccodes(11:numel(HH.rxns))=[];
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
         HH.metMiriams(index)=[];
         HH.inchis(index)=[];
         HH.metComps(index)=[];
         HH.b(index)=[];
         HH.S(index,:)=[];
         
         HH3=HH;
         HH=HH3;
         k=1; 
         index=0;
         HH.rxnGeneMat(11:numel(HH1.rxns),:)=[];
for i=1:numel(HH1.genes)
    i
      aa=full(HH.rxnGeneMat(:,i));
      Index=find(aa==0);
     if numel(Index)==10
         index(k)=i;
         k=k+1;
     end
ends
HH.genes(index)=[];
HH.rxnGeneMat(:,index)=[];
HH.geneComps(index)=[];
