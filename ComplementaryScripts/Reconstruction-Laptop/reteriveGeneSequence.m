% The aim of this script is to identify DNA sequence from Transcription
% Start Side (TSS) until the end of gene
% The reference genome is R64-1-1 March 2011, because the TSS data was
% computed based of this reference. The genome reference was download from 
% http://downloads.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
% We then uncompress this file and we got the gene boundiress from file gff
% file. After that we filter file with gene feature and mapped TSS from NAR
% paper which has Pubmed id xxxxxxx. The final processed file is stored in
% file geneLocation.xls.
%
% The reference genome is  S288C_reference_sequence_R64-1-1_20110203.fsa
% Here we try to compile the both files geneLocation.xls and 
% S288C_reference_sequence_R64-1-1_20110203.fsa to extract the 5'UTR
% sequence.
%
% Written by Ibrahim Elsemman
% DTU Denmark
% 19 - 10 - 2015
clear all

% [a b ctRNA]=xlsread('tRNA_modification_position_1.xlsx','charged');
% codon=ctRNA(2:65,1);
% uncharged=ctRNA(2:65,2);
% charged=ctRNA(2:65,3);

% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start Extracting\n');

fptr=fopen('utr_seq.txt','w');
fptr1=fopen('reactions.txt','w');

[a b geneData]=xlsread('TableS1.xlsx','Annotation');
n=numel(geneData(:,7));

for i=2:1504 %n
    i
   gene=cell2mat(geneData(i,2));
   Start=cell2mat(geneData(i,12));
   End = cell2mat(geneData(i,13));
   direction = cell2mat(geneData(i,5));
   TSS=cell2mat(geneData(i,11));
   if ~(strcmp(cell2mat(geneData(i,4)),'chrM'))
      chr= getChromosomeNumber(cell2mat(geneData(i,4)));

   if direction=='+'
       if TSS==0
          utr_length=50;
       else
          if TSS<Start % TSS mab be in CDS
             utr_length =Start-TSS;
          else
            utr_length = 0;
          end
       end
       DNAsequence{i-1,1}=yeastgenome(chr).Sequence(Start:End);
   else
        if TSS==0
          utr_length=50;
        else
           if TSS > End % TSS mab be in CDS
               utr_length = TSS-End;
          else
            utr_length = 0;
          end
       end
       DNAsequence{i-1,1}=seqcomplement(yeastgenome(chr).Sequence(End:-1:Start));
   end
 
 %intial 
 step_average = randomScanning(utr_length,.9,0.002);
 nATP_intial=1+step_average;
 nGTP_intial=1;
 %elongation
 seq=DNAsequence{i-1,1};
 elnogation_seq=seq(4:numel(seq)-3);
 [eq_substrates eq_product nATP_elnogation nGTP_elnogation] = elnogation(elnogation_seq);       
 
 %termination
 nATP_termination=1;
 nGTP_termination=2;
 total_atp=nATP_intial + nATP_elnogation + nATP_termination;
 total_gtp= nGTP_intial + nGTP_elnogation + nGTP_termination;
 total_h2o = total_atp + total_gtp;
 total_pi = total_atp + total_gtp;
 total_h = total_atp + total_gtp;
 eq=sprintf('\tbiosynthesis of %s Peptide\t%d H2O [cytoplasm] + %d ATP [cytoplasm] + %d GTP [cytoplasm] + Met1-tRNA(AUG) [cytoplasm] + %s => %s_peptide [cytoplasm] + %d ADP [cytoplasm] + %d GDP [cytoplasm] + %d phosphate [cytoplasm] + %d H+ [cytoplasm] + 1-tRNA(AUG) [cytoplasm] + %s',gene, total_h2o, total_atp,total_gtp,eq_substrates,gene,total_atp,total_gtp,total_pi,total_h,eq_product);
 fprintf(fptr1,'%s\n',eq);
 fprintf(fptr1,'\tfolding of %s Peptide\t%s_peptide [cytoplasm] => %s_folding [cytoplasm]\n',gene,gene,gene);
 
   end % end if chrM
     
   end

fclose(fptr1);
fprintf('finishing\n');

