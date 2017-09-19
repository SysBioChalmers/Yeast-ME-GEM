function [comp S1 P1 lb ub C]=getCompartment(model,rxn_index)
mets_index=find(model.S(:,rxn_index));
mets=model.metNames(mets_index);
lb=model.lb(rxn_index);
ub=model.ub(rxn_index);
C = model.c(rxn_index);
S =find(model.S(:,rxn_index)<0);
P =find(model.S(:,rxn_index)>0);

aS =full(model.S(S,rxn_index));
aP =full(model.S(P,rxn_index));

eq= sprintf('%f %s',abs(aS(1)),cell2mat(model.metNames(S(1))));

for i=2:numel(S)
     eq=sprintf('%s + %f %s',eq,abs(aS(i)),cell2mat(model.metNames(S(i))));
    
    
   
end
S1=eq;
if numel(P)>0
    eq= sprintf('%f %s',aP(1),cell2mat(model.metNames(P(1))));;
      
    for i=2:numel(P)
        eq=sprintf('%s + %f %s',eq,aP(i),cell2mat(model.metNames(P(i))));
    end
    P1=eq;
else
    eq=sprintf('%s =>',eq);
    P1='';
end
for i=1:numel(mets)
    c=regexp(mets{i},'extracellular','match');
    if isempty(c)==0
        comp='extracellular';
    end
end


for i=1:numel(mets)
    c=regexp(mets{i},'cytoplasm','match');
    if isempty(c)==0
        comp='cytoplasm';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'mitochondrial membrane','match');
    if isempty(c)==0
        comp='mitochondrial membrane';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'mitochondrion','match');
    if isempty(c)==0
        comp='mitochondrion';
    end
end


for i=1:numel(mets)
    c=regexp(mets{i},'endoplasmic reticulum','match');
    if isempty(c)==0
        comp='endoplasmic reticulum';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'peroxisome','match');
    if isempty(c)==0
        comp='peroxisome';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'nucleus','match');
    if isempty(c)==0
        comp='nucleus';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'Golgi','match');
    if isempty(c)==0
        comp='Golgi';
    end
end



for i=1:numel(mets)
    c=regexp(mets{i},'vacuole','match');
    if isempty(c)==0
        comp='vacuole';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'lipid particle','match');
    if isempty(c)==0
        comp='lipid particle';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'cell envelope','match');
    if isempty(c)==0
        comp='cell envelope';
    end
end

for i=1:numel(mets)
    c=regexp(mets{i},'vacuolar membrane','match');
    if isempty(c)==0
        comp='vacuolar membrane';
    end
end
