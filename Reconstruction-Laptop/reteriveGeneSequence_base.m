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
% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start Extracting\n');

fptr=fopen('utr_seq.txt','w');
fptr1=fopen('reactions.txt','w');

sheets={'Metabolic_genes','NPC'};
titles={'Metabolic','Nuclear_pore_complex'};
for k=1:2
    [a b geneData]=xlsread('geneLocation.xlsx',cell2mat(sheets(k)));
    n=numel(geneData(:,7));

for i=2:n
   gene=cell2mat(geneData(i,1));
   Start=cell2mat(geneData(i,3));
   End = cell2mat(geneData(i,4));
   direction = cell2mat(geneData(i,6));
   TSS= cell2mat(geneData(i,2));
   PAS= cell2mat(geneData(i,5));
   TSS_length= cell2mat(geneData(i,7));
   PAS_length= cell2mat(geneData(i,8));
   chr= getChromosomeNumber(cell2mat(geneData(i,9)));
   

   if direction=='+'
              
          DNAsequence{i-1,1}=yeastgenome(chr).Sequence(Start-TSS_length:End+PAS_length);
        % utrSequence{i-1,1}=yeastgenome(chr).Sequence(Start-TSS_length:Start-1);
       
   else
      
       DNAsequence{i-1,1}=seqcomplement(yeastgenome(chr).Sequence(End+TSS_length:-1:Start-PAS_length));
    %   utrSequence{i-1,1}= seqcomplement(yeastgenome(chr).Sequence(TSS:-1:End-1));
       
   end
  
       r=basecount(DNAsequence{i-1,1});
       
       fprintf(fptr1,'R_%s_transcript_%s\tTranscription for %s gene %s\t%d ATP [nucleus] + %d CTP [nucleus] + %d GTP [nucleus] + %d UTP [nucleus] => mRNA_transcript_%s [nucleus] + %d diphosphate [nucleus]\t\t\t\t\t\t\t\t\t\n',titles{k},gene,titles{k},gene, r.A,r.C,r.G,r.T,gene,r.A+r.C+r.G+r.T);

       fprintf(fptr,'%s\t%s\t%d\t%d\t%d\t%d\n',gene,DNAsequence{i-1,1},Start,End,TSS_length,PAS_length);
     
     
   end
end
fclose(fptr);
fprintf('finishing\n');

