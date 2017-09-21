function number = peptideRedius()
%In this function we compute the number peptide in mitochondrial intermembrane
% We used The redius of sphere containing the protein molecular weight
cell_mito_ratio=.015;
V=100; % as micro meter power 3 \mum^3
V=V*cell_mito_ratio;
sphere_redius = ((3*V)/(4*pi))^(1/3);
mito_surface=10%4*pi*sphere_redius^2;

M=29100.1+12*7760.5; %5825.3; %Da
r_min=0.066*M^(1/3)/1000;
area_ATP9=pi*r_min^2;
total_atp9=area_ATP9;%*10;
number_f1_f0=mito_surface*.15/total_atp8