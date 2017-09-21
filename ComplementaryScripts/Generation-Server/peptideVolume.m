function vol = peptideVolume(M)
%In this function we compute the number peptide in mitochondrial intermembrane
% We used The redius of sphere containing the protein molecular weight
% http://biologicalproceduresonline.biomedcentral.com/articles/10.1007/s12575-009-9008-x

r_min=0.066*M^(1/3)/1000; %convert nm to mu m
vol=(4/3)*pi*r_min^3;
