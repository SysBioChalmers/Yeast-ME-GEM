function M=getModificationMito(aa,codon,nucleus_compartment,cytosol_compartment,Position,Reaction)

aa_index=find(ismember(Position(:,1),aa));
codon_index=find(ismember(Position(aa_index,2),codon));

aa_modification=Position(aa_index(codon_index),4);

comp={};
k=1;
M=aa_index(codon_index);
% for k=1:numel(index)
%     %k
%     i=index(k);
%     modification_index=find(ismember(Reaction(1:58,4),Position(i,4)));
%     modification_position=cell2mat(Reaction(modification_index,3));
%     J=find(modification_position==cell2mat(Position(i,3)));
%     %reaction index
%     r=modification_index(J);
%     
%     comp{k,1}=changeCompartmentFormat(cell2mat(Reaction(r,2)));
% end
% 
% N=index(find(ismember(comp,nucleus_compartment)));
% C=index(find(ismember(comp,cytosol_compartment)));
