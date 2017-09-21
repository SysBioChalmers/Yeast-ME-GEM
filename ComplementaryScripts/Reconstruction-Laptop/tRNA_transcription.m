clc
fptr=fopen('tRNA_reactions_output.txt','w');
[n a data]=xlsread('tRNA_modification_position.xlsx','Sequence');
for i=2:42 % loop for tRNA
   seq=data(i,9);
   bc=basecount(dna2rna(cell2mat(seq)));
   if cell2mat(data(i,10))==1
       fprintf(fptr,'tRNA_%s%s%s%s\tRNA_Pol_III[n] + TFIIIB[n] + TFIIIC[n] + %d%s%d%s%d%s%d%s%s%s%s%s%d%s\n',cell2mat(data(i,3)),'_',cell2mat(data(i,4)),'_transcription',bc.A,' ATP[n] + ',bc.C,' CTP[n] + ',bc.G,' GTP[n] + ',bc.T,' UTP[n] => tRNA_',cell2mat(data(2,3)),'_',cell2mat(data(i,4)),'_pre_mature_transcribt[n] + ',(bc.A + bc.C + bc.G + bc.T), ' PPi + tRNA_Pol_III_inactive[n] + TFIIIB_inactive[n] + TFIIIC_inactive[n]');
       seq=cell2mat(data(i,9));
       from=cell2mat(data(i,11));
       to=cell2mat(data(i,12));
       intron=seq(from:to);
       bc=basecount(intron);
       fprintf(fptr,'tRNA_%s_%s%s\ttRNA_%s_%s%stRNA_%s_%s%s + %d ATP[n] + %d CTP[n] + %d GTP[n] + %d UTP[n] + tRNA_splicing_complex_inactive[n]\n',cell2mat(data(i,3)),cell2mat(data(i,4)),'_splicing',cell2mat(data(i,3)),cell2mat(data(i,4)),'_pre_mature_transcript[n] + tRNA_splicing_complex_inactive[n] + ATP[n] + H2O[n] + NAD+[n] + tRNA_splicing => PPi[n] + AMP[n] + Arr>p[n] + ',cell2mat(data(2,3)),cell2mat(data(2,4)),'_mature_transcribt[n]',bc.A,bc.C,bc.G,bc.T);
   else
      fprintf(fptr,'tRNA_%s%s%s%s\tRNA_Pol_III[n] + TFIIIB[n] + TFIIIC[n] + %d%s%d%s%d%s%d%s%s%s%s%s%d%s\n',cell2mat(data(i,3)),'_',cell2mat(data(i,4)),'_transcription',bc.A,' ATP[n] + ',bc.C,' CTP[n] + ',bc.G,' GTP[n] + ',bc.T,' UTP[n] => tRNA_',cell2mat(data(2,3)),'_',cell2mat(data(i,4)),'_mature_transcript[n] + ',(bc.A + bc.C + bc.G + bc.T), ' PPi + tRNA_Pol_III_inactive[n] + TFIIIB_inactive[n] + TFIIIC_inactive[n]');
      
   end
   %tRNA transporting
   fprintf(fptr,'tRNA_%s_%s_export\ttRNA_%s_%s_mature_transcript[n] + Ran-GTP[n] + Karyopherins-tRNA[n] + Nuclear_Pore_Complex_tRNA => Ran-GTP[c] + Karyopherins-tRNA[c] + Pi[c] + Nuclear_Pore_Complex + Rna1p[c] + tRNA_%s_%s_mature_transcript[c]\n',cell2mat(data(i,3)),cell2mat(data(i,4)),cell2mat(data(i,3)),cell2mat(data(i,4)),cell2mat(data(i,3)),cell2mat(data(i,4)));
   

end
fclose(fptr);