function fiveUTRLength(ng,mg)

% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start Extracting\n');

fptr=fopen('utr_Length.txt','w');

[a b geneData]=xlsread('TableS1.xlsx','Annotation');
n=numel(geneData(:,7));

for i=ng:mg
   i
   gene  = cell2mat(geneData(i,2));
   Start = cell2mat(geneData(i,12));
   End   = cell2mat(geneData(i,13));
   direction = cell2mat(geneData(i,5));
   
   TSS=cell2mat(geneData(i,11));
   PAL=cell2mat(geneData(i,14));
   
   intron=cell2mat(geneData(i,10));
   protein_seq=cell2mat(geneData(i,9));
   if ~(strcmp(cell2mat(geneData(i,4)),'chrM'))
      chr= getChromosomeNumber(cell2mat(geneData(i,4)));

   if direction=='+'
       if TSS==0
          utr_length=50;
          three_utr_length=50;
       else
          if TSS<Start % TSS may be in CDS
             utr_length =Start-TSS;
             three_utr_length=PAL-End;
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
               three_utr_length=Start-PAL;
          else
            utr_length = 0;
          end
       end
       DNAsequence{i-1,1}=seqcomplement(yeastgenome(chr).Sequence(End:-1:Start));
   end
 % Remove intron
 if strcmp(gene,'YDRUMD1')==0
      seq=DNAsequence{i-1,1};
 else
     seq=  cell2mat(geneData(i,8));
 end
 
 if intron==1
     %search for inton information in the Intron sheet
     [a b intronData]=xlsread('TableS1.xlsx','Intron');
     index=find(ismember(intronData(:,2),gene));
     a1=0;
     for k=1:numel(index)
         intronStart = cell2mat(intronData(index(k),14));
         intronEnd = cell2mat(intronData(index(k),15));
          if direction=='+'
              a=intronStart - Start-a1;
              b=intronEnd - Start-a1;
              intron_seq=seq(a+1:b+1);
              a1=numel(a+1:b+1);
              seq(a+1:b+1)=[];
          else
              a= End - intronStart-a1;
              b= End - intronEnd-a1;
              intron_seq=seq(b+1:a+1);
               a1=numel(b+1:a+1);
              seq(b+1:a+1)=[];
          end
         
          if(strcmp(intron_seq,cell2mat(intronData(index(k),12)) )==0)
              msg=sprintf('the gene %s has error in intron',gene);
              error(msg);
          end
     end
 end
 if strcmp(nt2aa(seq),protein_seq)==0
              fprintf('the gene %s has error in sequence\n',gene);
 end
 %intial
 fprintf(fptr,'%s\t%d\t%d\n',gene,utr_length,three_utr_length);
   else
       fprintf(fptr,'%s\t%d\n',gene,0); 
   end 
end


fclose(fptr);
fprintf('finishing\n');

