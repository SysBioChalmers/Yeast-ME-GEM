clear mEq nmEq
k=1;
for i=1:numel(model.rxns)
    i
rxn=sprintf('%s_forward',cell2mat(model.rxns(i)));
ii=find(ismember(newModel.rxns,rxn));
ms=model.metNames(find(model.S(:,i)<0));
mp=model.metNames(find(model.S(:,i)>0));
s=newModel.metNames(find(newModel.S(:,ii(1))<0));
p=newModel.metNames(find(newModel.S(:,ii(1))>0));
cs=newModel.metComps(find(newModel.S(:,ii(1))<0));
cp=newModel.metComps(find(newModel.S(:,ii(1))>0));
eq1=sprintf('%s',cell2mat(ms(1)));
eq2=sprintf('%s[%s]',cell2mat(s(1)),cell2mat(cs(1)));
for j=2:numel(ms)
    eq1=sprintf('%s + %s',eq1,cell2mat(ms(j)));
    eq2=sprintf('%s + %s[%s]',eq2,cell2mat(s(j)),cell2mat(cs(j)));
end

eq11=sprintf('%s',cell2mat(mp(1)));
eq22=sprintf('%s[%s]',cell2mat(p(1)),cell2mat(cp(1)));

for j=2:numel(mp)
    eq11=sprintf('%s + %s',eq11,cell2mat(mp(j)));
    eq22=sprintf('%s + %s[%s]',eq22,cell2mat(p(j)),cell2mat(cp(j)));
end

if strcmp(eq1,eq2) && strcmp(eq11,eq22)
    comp=1;
   
else
    comp=0;
    mEq(k,1)={sprintf('%s => %s',eq1,eq11)};
    nmEq(k,1)={sprintf('%s => %s',eq2,eq22)};
    k=k+1;
end
end