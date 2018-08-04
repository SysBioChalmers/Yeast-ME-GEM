%model=MyModel;
figure
i=1;
key='D37M_YML120C';
clear vol_mito vol_mim vol_cyto Molecules_mim Molecules_mito Molecules_cyto
mu1=[0.25:0.01:0.38];
mu1=sort(mu1);
k=1;
model=MyModel;
for i=1:numel(mu1)
    mu=mu1(i);

fileName=sprintf('Max_Mu_GEKO_Glucose_%f_T33_Copies_0_k_1.lp.out',mu);

[vol_mito(i) vol_mim(i) Molecules_mito(i) Molecules_mim(i)] = mitoVolume(model, fileName,mu,0.05);
[vol_cyto(i) Molecules_cyto(i)]= cytosolVolume(model, fileName,mu,0.05);

end
mu=mu1;
subplot(2,2,1)
plot(mu,vol_cyto);
title('Volume Cytosol');

subplot(2,2,2)
plot(mu,vol_mito);
title('Volume mitochondrion');
subplot(2,2,3)
plot(mu,vol_mim);
title('Volume IMM');

subplot(2,2,4)
plot(mu,vol_mito+vol_cyto);
title('Volume C+M');
%molecules
subplot(2,4,5)
plot(mu,Molecules_cyto);
title('# molecules in Cytosol');
subplot(2,4,6)
plot(mu,Molecules_mito);
title('# molecules in mitochondrion');

subplot(2,4,7)
plot(mu,Molecules_mim);
title('# molecules in MIM');
subplot(2,4,8)
plot(mu,Molecules_mito+Molecules_cyto);
title('# molecules in C+M');