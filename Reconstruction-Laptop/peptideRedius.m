function number = peptideRedius(V,cell_mito_ratio,M)
%In this function we compute the number peptide in mitochondrial intermembrane
% We used The redius of sphere containing the protein molecular weight
% http://biologicalproceduresonline.biomedcentral.com/articles/10.1007/s12575-009-9008-x

V=100; % as micro meter power 3 \mum^3
V=V*cell_mito_ratio;
sphere_redius = ((3*V)/(4*pi))^(1/3);
mito_surface=4*pi*sphere_redius^2;

r_min=0.066*M^(1/3)/1000; %convert nm to mu m
area_ATP9=pi*r_min^2;
total_atp9=area_ATP9;%*10;
number=mito_surface*.15/total_atp;