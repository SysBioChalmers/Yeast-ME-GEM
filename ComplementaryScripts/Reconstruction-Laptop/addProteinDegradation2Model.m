function newModel= addProteinDegradation2Model(newModel,n,m)

% step 1: read the reference genome
yeastgenome=fastaread('S288C_reference_sequence_R64-1-1_20110203.fsa');

% step 2: read the gene location excel file

% please note a and b are unimportant
fprintf('Start Extracting\n');

fptr=fopen('utr_seq.txt','w');
fptr1=fopen('reactions.txt','w');

[a b geneData]=xlsread('TableS1.xlsx','Annotation');


for i=n:m
    
   gene  = cell2mat(geneData(i,2));
   Start = cell2mat(geneData(i,12));
   End   = cell2mat(geneData(i,13));
   direction = cell2mat(geneData(i,5));
   TSS=cell2mat(geneData(i,11));
   intron=cell2mat(geneData(i,10));
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
 % Remove intron
 seq=DNAsequence{i-1,1};
 
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
 
[eq_product number_of_aa]= degrade(seq);       
 
 %termination
nATP =floor(1.3*number_of_aa);
eq=sprintf('%d H2O [cytoplasm] + %d ATP [cytoplasm] + %s_subunit [cytoplasm] => %s + %d ADP [cytoplasm] + %d phosphate [cytoplasm] + %d H+ [cytoplasm]',nATP, nATP,gene,eq_product,nATP,nATP,nATP);

rxnID=sprintf('%s_subunit_degradation',gene);
 rxnID=strrep(rxnID,'-','');
 rxnName=sprintf('degradation of %s Peptide',gene);

 newModel=addYeastReaction(newModel,eq,{rxnID},{rxnName},0,1000,0,{'Proteasome'});

 end % end if chrM
     
end

fclose(fptr1);
fprintf('finishing\n');

