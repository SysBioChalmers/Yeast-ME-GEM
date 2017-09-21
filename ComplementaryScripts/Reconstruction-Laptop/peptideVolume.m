function v = peptideVolume(M)
%In this function we compute the number peptide in mitochondrial intermembrane
% We used The redius of sphere containing the protein molecular weight
% http://biologicalproceduresonline.biomedcentral.com/articles/10.1007/s12575-009-9008-x

v= (1.212e-3 * M)/1e9; %the unit is nm^3
