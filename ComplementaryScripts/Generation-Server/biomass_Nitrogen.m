function eq= biomass_Nitrogen(mu)
%This function returns the biomass reaction for the strain CEN.PK 113-7D
% Ibrahim Elsemman

[a b lipidComposition]=xlsread('TableS1.xlsx','Lipid');
[a b AAComposition]=xlsread('TableS1.xlsx','AminoAcids');

%lipid 

lipid_ratio =(4.8)/100;
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
    H2O= sprintf('%.15f H2O [cytoplasm]', H2O_ratio*1000/18.01528);
    
 %RNA
    RNA_ratio =(3.7)/100;%(-14.883*mu*mu + 21.131*mu + 4.6712)/100;
    factor=RNA_ratio/0.063;
    
    AMP=0.051*factor; 
    GMP=0.051*factor;
    UMP=0.067*factor;
    CMP=0.05*factor;
    RNA=sprintf('%f AMP [cytoplasm] + %f GMP [cytoplasm] + %f CMP [cytoplasm] + %f UMP [cytoplasm]',AMP,GMP,UMP,CMP); 
    
 %DNA
    DNA='0.00359 dAMP [cytoplasm] + 0.00243 dCMP [cytoplasm] + 0.00243 dGMP [cytoplasm] + 0.00359 dTMP [cytoplasm]';
    DNA_ratio=0.005;
 %phosphate
    phosphate_ratio=0;
     phosphate =sprintf('%f phosphate [cytoplasm]', phosphate_ratio*1000/94.9714);
% carbohydrates      
   
    glucan_MW   = 162;
    glycogen_MW = 162;
    mannan_MW   = 162;
    trehalose_MW   = 342.3;
    
    glucan_ratio=(59*.28)/100;
    glycogen_ratio= (59*.1952)/100;
    mannan_ratio= (59*.28)/100;
    trehalose_ratio =  (59*.2428)/100;
    
    glucan = glucan_ratio *1000/glucan_MW;
    glycogen = glycogen_ratio *1000/glycogen_MW;
    mannan = mannan_ratio *1000/mannan_MW;
    trehalose=trehalose_ratio *1000/trehalose_MW;
    
    carbohydrates =sprintf('%.6f (1-_GT_3)-beta-D-glucan [cytoplasm] + %.6f mannan [cytoplasm] + %.6f glycogen [cytoplasm] + %.6f trehalose [cytoplasm]',glucan,mannan,glycogen, trehalose);
% sulfate
    sulfate_ratio = 0.000;
    sulfate_MW = 96.063; %gram / mol
    sulfate=sprintf('%.6f sulphate [cytoplasm]', 1000/sulfate_MW *sulfate_ratio);

% riboflavin
    riboflavin_ratio  = 0;
    riboflavin_MW = 376.3639; %gram /mol
    riboflavin=sprintf('%.6f riboflavin [cytoplasm]', 1000/riboflavin_MW *riboflavin_ratio);
 % heme a
    heme_ratio = 0.00;
    heme_MW = 852.837;
    heme= sprintf('%.6f heme a [cytoplasm]',1000/heme_MW *heme_ratio);
    
 % chitin
    chitin_ratio = 0.00;
    chitin_MW = 683.311282535;
    chitin= sprintf('%.6f chitin [cytoplasm]',2.03E-7);
    
    % free amino acids
free_aa_ratio =0.04;
%total_gram = total_gram + free_aa_ratio;
free_aa='';
    for i=2:numel(AAComposition(:,1))-1
        c= cell2mat(AAComposition(i,3))*free_aa_ratio*1000/cell2mat(AAComposition(i,4));
        if i==2
           free_aa=sprintf('%.15f %s',c,cell2mat(AAComposition(i,1)));
        else
             free_aa=sprintf('%s + %.15f %s',free_aa,c,cell2mat(AAComposition(i,1)));
        end
    end

    
    eq=sprintf('%s + %s + %s + %s + %s + %s + %s + %s + %s + %s + %s => biomass_test [cytoplasm]',free_aa, H2O, lipid, carbohydrates,RNA,DNA, phosphate, sulfate, riboflavin,heme, chitin);
    
