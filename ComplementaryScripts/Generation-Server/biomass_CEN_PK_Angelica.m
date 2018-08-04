function [eq total_gram free_aa_ratio]= biomass_CEN_PK_Angelica(mu)
%This function returns the biomass reaction for the strain CEN.PK 113-7D
% Ibrahim Elsemman

[a b lipidComposition]=xlsread('TableS1.xlsx','Lipid');
[a b AAComposition]=xlsread('TableS1.xlsx','AminoAcids');

%lipid 
total_gram=0;
lipid_ratio =(-4.5576*mu + 7.4549)/100;
total_gram = total_gram + lipid_ratio;
lipid='';
    for i=2:numel(lipidComposition(:,1))-1
        c= cell2mat(lipidComposition(i,5))*lipid_ratio*1000/cell2mat(lipidComposition(i,3));
        if i==2
           lipid=sprintf('%.15f %s',c,cell2mat(lipidComposition(i,1)));
        else
             lipid=sprintf('%s + %.15f %s',lipid,c,cell2mat(lipidComposition(i,1)));
        end
    end

 
 % H2O
    H2O_ratio= (0)/100;
%    H2O= sprintf('%.15f H2O [cytoplasm]', H2O_ratio*1000/18.01528);
    H2O= sprintf('%.15f H2O [cytoplasm]', 0);
    
 %RNA
    RNA_ratio =(144.38*mu + 62.075)/1000;
    total_gram = total_gram + RNA_ratio;
    factor=RNA_ratio/0.063;
    
    AMP=0.051*factor; 
    GMP=0.051*factor;
    UMP=0.067*factor;
    CMP=0.05*factor;
    RNA=sprintf('%f AMP [cytoplasm] + %f GMP [cytoplasm] + %f CMP [cytoplasm] + %f UMP [cytoplasm]',AMP,GMP,UMP,CMP); 
    
 %DNA
    DNA='0.00359 dAMP [cytoplasm] + 0.00243 dCMP [cytoplasm] + 0.00243 dGMP [cytoplasm] + 0.00359 dTMP [cytoplasm]';
    DNA_ratio=0.005;
    total_gram = total_gram + DNA_ratio;
    
 %phosphate
    phosphate_ratio=0;
     phosphate =sprintf('%f phosphate [cytoplasm]', phosphate_ratio*1000/94.9714);
      total_gram = total_gram + phosphate_ratio;
% sulfate
    sulfate_ratio =  0.0019;
    total_gram = total_gram + sulfate_ratio;

    sulfate_MW = 96.063; %gram / mol
    sulfate=sprintf('%.6f sulphate [cytoplasm]', 1000/sulfate_MW *sulfate_ratio);

% riboflavin
    riboflavin_ratio  = 3.3873e-04;
    total_gram = total_gram + riboflavin_ratio;

    riboflavin_MW = 376.3639; %gram /mol
    riboflavin=sprintf('%.6f riboflavin [cytoplasm]', riboflavin_ratio*1000/riboflavin_MW);
 % heme a
    heme_ratio = 0;
    total_gram = total_gram + heme_ratio;

    heme_MW = 852.837;
    heme= sprintf('%.6f heme a [cytoplasm]',1000/heme_MW *heme_ratio);
    
 % chitin
    chitin_ratio = 2.03E-7;
    total_gram = total_gram + chitin_ratio;

    chitin_MW = 221.2;
    chitin= sprintf('%.15f chitin [cytoplasm]',1000/chitin_MW *chitin_ratio);
    
 % free amino acids
 % carbohydrates
    glucan_MW   = 162;
    glycogen_MW = 162;
    mannan_MW   = 162;
    trehalose_MW   = 342.3;
    
    glucan_ratio=(540117*mu^6 - 745285*mu^5 + 401048*mu^4 - 104290*mu^3 + 12921*mu^2 - 611.55*mu + 19.633)/100;
    glycogen_ratio= (-581.6*mu^3 + 488.33*mu^2 - 148.74*mu + 18.095)/100;
    mannan_ratio= (-1.8223*mu + 10.367)/100;
    trehalose_ratio = (-43.601*mu + 9.2581)/100;
    if trehalose_ratio <0
        trehalose_ratio =0;
    end
          deference= (1- (proteinRatio_Angelica(mu) + total_gram + glucan_ratio + glycogen_ratio + mannan_ratio + trehalose_ratio));
          if deference<0
              glucan_ratio = glucan_ratio - 0.5*abs(deference);
              mannan_ratio = mannan_ratio - 0.5*abs(deference);
          end
          
          total_gram = total_gram + glucan_ratio + glycogen_ratio + mannan_ratio + trehalose_ratio;
    
    glucan = glucan_ratio *1000/glucan_MW;
    glycogen = glycogen_ratio *1000/glycogen_MW;
    mannan = mannan_ratio *1000/mannan_MW;
    trehalose=trehalose_ratio *1000/trehalose_MW;
    
    
      
    
    
    carbohydrates =sprintf('%.6f (1-_GT_3)-beta-D-glucan [cytoplasm] + %.6f mannan [cytoplasm] + %.6f glycogen [cytoplasm] + %.6f trehalose [cytoplasm]',glucan,mannan,glycogen, trehalose);

    
free_aa_ratio =1-(proteinRatio_Angelica(mu) + total_gram);
if free_aa_ratio<1e-6
    free_aa_ratio=0;
end
total_gram = total_gram + free_aa_ratio;
free_aa='';
    for i=2:numel(AAComposition(:,1))-1
        c= cell2mat(AAComposition(i,3))*free_aa_ratio*1000/cell2mat(AAComposition(i,4));
        if i==2
           free_aa=sprintf('%.15f %s',c,cell2mat(AAComposition(i,1)));
        else
             free_aa=sprintf('%s + %.15f %s',free_aa,c,cell2mat(AAComposition(i,1)));
        end
    end
    
    eq=sprintf('%s + %s + %s + %s + %s + %s + %s + %s + %s + %s => biomass_test [cytoplasm]',free_aa, lipid, carbohydrates,RNA,DNA, phosphate, sulfate, riboflavin,heme, chitin);
    
   total_gram = total_gram + proteinRatio_Angelica(mu);