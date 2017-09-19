clc
 fptr=fopen('tRNA_reactions_output.txt','w');
 fptr1=fopen('charged_tRNA_names.txt','w');
 [n a tRNAsequence]=xlsread('tRNA_modification_position_1.xlsx','Sequence');
 [n a Position]=xlsread('tRNA_modification_position_1.xlsx','Position');
 [n a Reaction]=xlsread('tRNA_modification_position_1.xlsx','Reaction');


% compartment in the model
nucleus_compartment = '[n]';
cytosol_compartment ='[c]';

for i=2:43 % loop for tRNA
   i
   aa=cell2mat(tRNAsequence(i,3));
   codon=cell2mat(tRNAsequence(i,4));
   [N C]=getModificationCompartment(aa,codon,nucleus_compartment,cytosol_compartment,Position,Reaction);
   
   
   seq=tRNAsequence(i,10);
   bc=basecount(dna2rna(cell2mat(seq)));
   
   codon_names=fieldnames(bc);
   
   fprintf(fptr,'tRNA_%s%s%s%s\tRNA_Pol_III[n] + TFIIIB[n] + TFIIIC[n] + %d%s%d%s%d%s%d%s%s%s%s%s%d%s\n',cell2mat(tRNAsequence(i,3)),'_',cell2mat(tRNAsequence(i,4)),'_transcription',bc.A,' ATP[n] + ',bc.C,' CTP[n] + ',bc.G,' GTP[n] + ',bc.T,' UTP[n] => tRNA_',cell2mat(tRNAsequence(2,3)),'_',cell2mat(tRNAsequence(i,4)),'_transcript[n] + ',(bc.A + bc.C + bc.G + bc.T), ' PPi + tRNA_Pol_III_inactive[n] + TFIIIB_inactive[n] + TFIIIC_inactive[n]');
   base=sprintf('tRNA_%s_%s_transcript',cell2mat(tRNAsequence(i,3)),cell2mat(tRNAsequence(i,4)));
   
   [rxnID rxnNames eq modification]=tRNAreactions(base,N,Reaction,Position);
   %print equation into file
   for j=1:numel(eq)
       fprintf(fptr,'\t%s\n',cell2mat(eq(j)));
   end
   
   fprintf(fptr,'\t%s%s + Ran-GTP[n] + Karyopherins-tRNA[n] + Nuclear_Pore_Complex_tRNA + ATP[n] + H2O[n] => ADP[c] + 2 Pi[c] + 2 H[c] + Ran-GDP[c] + Karyopherins-tRNA[c] + Nuclear_Pore_Complex + Rna1p[c] + %s%s\n',modification_base,nucleus_compartment,modification_base,cytosol_compartment);
   base=modification_base;
   
   %add CCA to tRNA
   fprintf(fptr,'\t%s%s + ATP[c] + 2 CTP[c] => %s_CCA%s + 3 diphosphate\n',modification_base,cytosol_compartment,modification_base,cytosol_compartment);
   
   base=sprintf('%s_CCA',base);
   [rxnID rxnNames eq modification]=tRNAreactions(base,C,Reaction,Position);
   
   for j=1:numel(eq)
       fprintf(fptr,'\t%s\n',cell2mat(eq(j)));
   end
   
   %splicing and amino-acyl
   AA=cell2mat(tRNAsequence(i,13));
   
   if cell2mat(tRNAsequence(i,10))==1
       seq=cell2mat(tRNAsequence(i,9));
       from=cell2mat(tRNAsequence(i,11));
       to=cell2mat(tRNAsequence(i,12));
       intron=seq(from:to);
       bc=basecount(intron);
       
       fprintf(fptr,'\t%s%s + tRNA_splicing_complex[m] => %s_splicing%s + PPi[n] + AMP[n] + Arr_GT_p[n] + %d ATP[c] + %d CTP[c] + %d GTP[c] + %d UTP[c]\n',modification_base,cytosol_compartment,modification_base,cytosol_compartment,bc.A,bc.C,bc.G,bc.T);
       fprintf(fptr,'\t%s_splicing%s + H2O[c] + ATP[c] + %s%s => AMP[c] + diphosphate[c] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,codon,cytosol_compartment);
       fprintf(fptr1,'%s_splicing\t%s-tRNA(%s)\n',modification_base,aa,codon);
   else
       fprintf(fptr,'\t%s%s + H2O[c] + ATP[c] + %s%s => AMP[c] + diphosphate[c] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,codon,cytosol_compartment)
       fprintf(fptr1,'%s\t%s-tRNA(%s)\n',modification_base,aa,codon);
   end
   
   
end

fclose(fptr);
fclose(fptr1);