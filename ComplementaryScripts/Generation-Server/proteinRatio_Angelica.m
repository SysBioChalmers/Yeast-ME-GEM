function r=proteinRatio_Angelica(gr)
pr=[43.67732545 35.41934296 51.96100385 67.03996373 62.83782307 48.81608478];
mu=[0.2 0.23 0.27 0.3 0.32 0.34];
index=find(mu==gr);

r=pr(index)/100;

if gr>0.23 && gr <0.3
    r= (452.63*gr - 69.079)/100;
end

