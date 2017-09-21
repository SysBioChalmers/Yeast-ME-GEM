function model = addtRNA2Model(model) 
 fptr=fopen('tRNA_reactions_output.txt','w');
 fptr1=fopen('charged_tRNA_names.txt','w');
 [n a tRNAsequence]=xlsread('tRNA_modification_position_1.xlsx','Sequence');
 [n a Position]=xlsread('tRNA_modification_position_1.xlsx','Position');
 [n a Reaction]=xlsread('tRNA_modification_position_1.xlsx','Reaction');
 [n a RS]=xlsread('tRNA_modification_position_1.xlsx','AARS_Cytosol');


% compartment in the model
nucleus_compartment = '[n]';
cytosol_compartment ='[c]';

for i=2:43 % loop for tRNA
   aa=cell2mat(tRNAsequence(i,3));
   codon=cell2mat(tRNAsequence(i,4));
   antiCodons=regexp(cell2mat(tRNAsequence(i,5)),';','split');
   for ii=1:numel(antiCodons)
   antiCodon = cell2mat( antiCodons(ii));   
   [N C]=getModificationCompartment(aa,codon,nucleus_compartment,cytosol_compartment,Position,Reaction);
   
   %finding AARS
   RS_index=find(ismember(RS(:,1),aa));
   RS_gene=cell2mat(RS(RS_index,3));
   
   seq=tRNAsequence(i,10);
   bc=basecount(dna2rna(cell2mat(seq)));
   
   codon_names=fieldnames(bc);
   
   eq=sprintf('%d ATP [nucleus] + %d CTP [nucleus] + %d GTP [nucleus] + %d UTP [nucleus] => tRNA_%s_%s_%s_transcript [nucleus] + %d diphosphate [nucleus]',bc.A , bc.C , bc.G , bc.T,aa,codon,antiCodon,(bc.A + bc.C + bc.G + bc.T));
   base=sprintf('tRNA_%s_%s_%s_transcript',cell2mat(tRNAsequence(i,3)),cell2mat(tRNAsequence(i,4)),antiCodon);
  
   rxnID=sprintf('tRNA_%s_%s_%s_transcription',aa,codon,antiCodon);
   rxnNames=sprintf('tRNA_%s_%s_%s transcription',aa,codon,antiCodon);
   model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
      
   [rxnID rxnNames eq modification_base]=tRNAreactions_normal(base,antiCodon,N,Reaction,Position);
   %print equation into file
   for j=1:numel(eq)
       %fprintf(fptr,'\t%s\n',cell2mat(eq(j)));
       rxnID(j)=strrep(rxnID(j),'-','_');
       %fprintf('%s\n',cell2mat(eq(j)));
       model=addYeastReaction(model,cell2mat(eq(j)),rxnID(j),rxnNames(j),0,100,0,{''});
   end
   
   eq=sprintf('\t%s%s + GTP[n] + ATP[n] + 2 H2O[n] => ADP[c] + 2 phosphate[c] + 2 H[c] + %s%s\n',modification_base,nucleus_compartment,modification_base,cytosol_compartment);
   eq=strrep(eq,'[n]',' [nucleus]');
   eq=strrep(eq,'[c]',' [cytoplasm]');
   eq=strrep(eq,'+ H [','+ H+ [');
   eq=strrep(eq,'+ 2 H [','+ 2 H+ [');

   rxnID=sprintf('tRNA_%s_%s_%s_exporting',aa,codon,antiCodon);
   rxnNames=sprintf('tRNA_%s_%s_%s exporting',aa,codon,antiCodon);
   model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
   
   base=modification_base;
   
   %add CCA to tRNA
   eq=sprintf('\t%s%s + ATP[c] + 2 CTP[c] => %s_CCA%s + 3 diphosphate[c]\n',modification_base,cytosol_compartment,modification_base,cytosol_compartment);
       eq=strrep(eq,'[n]',' [nucleus]');
       eq=strrep(eq,'[c]',' [cytoplasm]');
       eq=strrep(eq,'+ H [','+ H+ [');
       eq=strrep(eq,'+ 2 H [','+ 2 H+ [');
       
   rxnID=sprintf('tRNA_%s_%s_adding_CCA',aa,codon);
   rxnNames=sprintf('tRNA_%s_%s adding CCA to the tRNA',aa,codon);
   model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
   
   
   base=sprintf('%s_CCA',base);
   [rxnID rxnNames eq modification_base]=tRNAreactions_normal(base,antiCodon,C,Reaction,Position);
   
   for j=1:numel(eq)
       %fprintf(fptr,'\t%s\n',cell2mat(eq(j)));
       eq(j)=strrep(eq(j),'[n]',' [nucleus]');
       eq(j)=strrep(eq(j),'[c]',' [cytoplasm]');
       eq(j)=strrep(eq(j),'+ H [','+ H+ [');
       eq(j)=strrep(eq(j),'2 H [','2 H+ [');
        eq(j)=strrep(eq(j),'5 H [','5 H+ [');
       eq(j)=strrep(eq(j),'SAM [','S-adenosyl-L-methionine [');
       eq(j)=strrep(eq(j),'SAH [','S-adenosyl-L-homocysteine [');
       
       rxnID(j)=strrep(rxnID(j),'-','_');
       
       %fprintf('%s\n',cell2mat(eq(j)));
       model=addYeastReaction(model,cell2mat(eq(j)),rxnID(j),rxnNames(j),0,100,0,{''});
     
   end
   
   %splicing and amino-acyl
   AA=cell2mat(tRNAsequence(i,14));
   %check for intron
   if cell2mat(tRNAsequence(i,11))==1
       seq=cell2mat(tRNAsequence(i,10));
       from=cell2mat(tRNAsequence(i,12));
       to=cell2mat(tRNAsequence(i,13));
       intron=seq(from:to);
       bc=basecount(intron);
       if bc.A >=1
          eq = sprintf('\t%s [cytoplasm] + NAD [cytoplasm] => %s_splicing [cytoplasm] + diphosphate [cytoplasm] + AMP [cytoplasm] + Arr_GT_p [cytoplasm] + %d ATP [cytoplasm] + %d CTP [cytoplasm] + %d GTP [cytoplasm] + %d UTP [cytoplasm]\n',modification_base,modification_base,(bc.A -1),bc.C,bc.G,bc.T);
       else
          eq = sprintf('\t%s [cytoplasm] + ATP [cytoplasm] + NAD [cytoplasm] => %s_splicing [cytoplasm] + diphosphate [cytoplasm] + AMP [nucleus] + Arr_GT_p [nucleus] + %d CTP [cytoplasm] + %d GTP [cytoplasm] + %d UTP [cytoplasm]\n',modification_base,modification_base,bc.C,bc.G,bc.T);
       end
       
       rxnID=sprintf('tRNA_%s_%s_splicing',aa,codon);
       rxnNames=sprintf('tRNA_%s_%s splicing',aa,codon);
       model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
           
       if strcmp(RS_gene,'YDR023W;YHR011W')==0
           eq= sprintf('\t%s_splicing%s + ATP[c] + %s%s => AMP[c] + diphosphate[c] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
           eq=strrep(eq,'[c]',' [cytoplasm]');
           rxnID=sprintf('%s_%s_loading',aa,codon);
           rxnNames=sprintf('%s_%s_loading',aa,codon);
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{RS_gene});
           
       else
           eq=sprintf('\t%s_splicing%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
           eq=strrep(eq,'[c]',' [cytoplasm]');
           rxnID=sprintf('%s_%s_loading_YDR023W',aa,codon);
           rxnNames=sprintf('%s_%s_loading',aa,codon);
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YDR023W'});
           
           eq=sprintf('\t%s_splicing%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
           eq=strrep(eq,'[c]',' [cytoplasm]');
           rxnID=sprintf('%s_%s_loading_YHR011W',aa,codon);
           rxnNames=sprintf('%s_%s_loading',aa,codon);
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YHR011W'});
       end
       
       
               eq=sprintf('tRNA(%s) [cytoplasm] => %s_splicing [cytoplasm]\n',antiCodon,modification_base);
               rxnID=sprintf('codon_synonym_%s_%s_relase',codon,antiCodon)
               rxnNames=sprintf('codon %s has %s synonym relase',codon,antiCodon)
               model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
               
               eq=sprintf('%s_splicing [cytoplasm] => \n',modification_base);
               rxnID=sprintf('codon_synonym_%s_%s_dilution',codon,antiCodon)
               rxnNames=sprintf('codon %s has %s synonym dilution',codon,antiCodon)
               model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
               
               
               fprintf(fptr1,'%s-tRNA(%s) [cytoplasm]\t%s_splicing [cytoplasm]\n',aa,antiCodon,modification_base);
               
                 
       
       
   else
       i
       if strcmp(RS_gene,'YBR121C;YPR081C')==0 && strcmp(RS_gene,'YDR023W;YHR011W')==0
           
       eq=sprintf('\t%s%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
       eq=strrep(eq,'[c]',' [cytoplasm]');
       rxnID=sprintf('%s_%s_%s_loading',aa,codon,antiCodon);
       rxnNames=sprintf('%s_%s_%s_loading',aa,codon,antiCodon);
       model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{RS_gene});
       %fprintf(fptr1,'%s\t%s-tRNA(%s)\n',modification_base,aa,codon);
       elseif strcmp(RS_gene,'YDR023W;YHR011W')==1
           eq=sprintf('\t%s%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
           eq=strrep(eq,'[c]',' [cytoplasm]');
           rxnID=sprintf('%s_%s_%s_loading_YDR023W',aa,codon,antiCodon);
           rxnNames=sprintf('%s_%s_%s_loading',aa,codon,antiCodon);
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YDR023W'});
           
           eq=sprintf('\t%s%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
           eq=strrep(eq,'[c]',' [cytoplasm]');
           rxnID=sprintf('%s_%s_%s_loading_YHR011W',aa,codon,codon,antiCodon);
           rxnNames=sprintf('%s_%s_%s_loading',aa,codon,antiCodon);
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YHR011W'});
           
       else
       eq=sprintf('\t%s%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
       eq=strrep(eq,'[c]',' [cytoplasm]');
       rxnID=sprintf('%s_%s_%s_loading_YBR121C',aa,codon,codon,antiCodon);
       rxnNames=sprintf('%s_%s_%s_loading',aa,codon,codon,antiCodon);
       model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YBR121C'});
       
       eq=sprintf('\t%s%s + ATP [cytoplasm] + %s%s => AMP [cytoplasm] + diphosphate [cytoplasm] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
       eq=strrep(eq,'[c]',' [cytoplasm]');
       rxnID=sprintf('%s_%s_%s_loading_YPR081C',aa,codon,codon,antiCodon);
       rxnNames=sprintf('%s_%s_%s_loading',aa,codon,codon,antiCodon);
       model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YPR081C'});
       end
       
                
           eq=sprintf('tRNA(%s) [cytoplasm] => %s [cytoplasm]\n',antiCodon,modification_base);
           rxnID=sprintf('codon_synonym_%s_%s_relase',codon,antiCodon)
           rxnNames=sprintf('codon %s has %s synonym relase',codon,antiCodon)
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});    
           
           eq=sprintf('%s [cytoplasm] => \n',modification_base);
           rxnID=sprintf('tRNA_%s_%s_dilution',codon,antiCodon)
           rxnNames=sprintf('tRNA_%s_%s_dilution',codon,antiCodon)
           model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});    
           
         
           
           fprintf(fptr1,'%s-tRNA(%s) [cytoplasm]\t%s [cytoplasm]\n',aa,antiCodon,modification_base);
      
      
   end
   
   
   end
end
fclose(fptr);
fclose(fptr1);