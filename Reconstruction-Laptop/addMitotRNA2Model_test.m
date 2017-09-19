function model = addMitotRNA2Model_test(model) 
 fptr=fopen('tRNA_reactions_output.txt','w');
 fptr1=fopen('charged_tRNA_names.txt','w');
 [n a tRNAsequence]=xlsread('tRNA_modification_position_1.xlsx','MitoSequence');
 [n a Position]=xlsread('tRNA_modification_position_1.xlsx','Position_mito');
 [n a Reaction]=xlsread('tRNA_modification_position_1.xlsx','Reaction');
 [n a RS]=xlsread('tRNA_modification_position_1.xlsx','AARS_Mito');

% compartment in the model
nucleus_compartment = '[n]';
cytosol_compartment ='[c]';

for i=2:25 % loop for tRNA
    
    aa=cell2mat(tRNAsequence(i,3));
    codon=cell2mat(tRNAsequence(i,4));
    antiCodons=regexp(cell2mat(tRNAsequence(i,5)),';','split');
    codon1=codon;
    for ii=1:numel(antiCodons)
        
        antiCodon = cell2mat( antiCodons(ii));
        
        codon=sprintf('%s_%s',codon,antiCodon);
        M =getModificationMito(aa,codon1,nucleus_compartment,cytosol_compartment,Position,Reaction);
        
        %finding AARS
        RS_index=find(ismember(RS(:,1),aa));
        RS_gene=cell2mat(RS(RS_index,3));
        if isnan(RS_gene)
            RS_gene='';
        end
        
        seq=tRNAsequence(i,9);
        bc=basecount(dna2rna(cell2mat(seq)));
        
        codon_names=fieldnames(bc);
        
        eq=sprintf('%d ATP [mitochondrion] + %d CTP [mitochondrion] + %d GTP [mitochondrion] + %d UTP [mitochondrion] => tRNA_%s_%s_transcript [mitochondrion] + %d diphosphate [mitochondrion]',bc.A , bc.C , bc.G , bc.T,aa,codon,(bc.A + bc.C + bc.G + bc.T));
        base=sprintf('tRNA_%s_%s_%s_transcript',cell2mat(tRNAsequence(i,3)),cell2mat(tRNAsequence(i,4)),antiCodon);
        
        rxnID=sprintf('tRNA_%s_%s_transcription_mitochondrion',aa,codon);
        rxnNames=sprintf('tRNA_%s_%s transcription in mitochondrion',aa,codon);
        model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
        
        [rxnID rxnNames eq modification_base]=tRNAreactionsMito(base,M,Reaction,Position,antiCodon);
        %print equation into file
        for j=1:numel(eq)
            %fprintf(fptr,'\t%s\n',cell2mat(eq(j)));
            rxnID(j)=strrep(rxnID(j),'-','_');
            %fprintf('%s\n',cell2mat(eq(j)));
            model=addYeastReaction(model,cell2mat(eq(j)),rxnID(j),rxnNames(j),0,100,0,{''});
        end
        
        %    eq=sprintf('\t%s%s + GTP[n] + ATP[n] + 2 H2O[n] => ADP[c] + 2 phosphate[c] + 2 H[c] + %s%s\n',modification_base,nucleus_compartment,modification_base,cytosol_compartment);
        %    eq=strrep(eq,'[n]',' [mitochondrion]');
        %    eq=strrep(eq,'[c]',' [mitochondrion]');
        %    eq=strrep(eq,'+ H [','+ H+ [');
        %    eq=strrep(eq,'+ 2 H [','+ 2 H+ [');
        %
        %    rxnID=sprintf('tRNA_%s_%s_exporting',aa,codon);
        %    rxnNames=sprintf('tRNA_%s_%s exporting',aa,codon);
        %    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
        
        base=modification_base;
        
        %add CCA to tRNA
        eq=sprintf('\t%s%s + ATP [mitochondrion] + 2 CTP [mitochondrion] => %s_CCA%s + 3 diphosphate [mitochondrion]\n',modification_base,cytosol_compartment,modification_base,cytosol_compartment);
        eq=strrep(eq,'[n]',' [mitochondrion]');
        eq=strrep(eq,'[c]',' [mitochondrion]');
        eq=strrep(eq,'+ H [','+ H+ [');
        eq=strrep(eq,'+ 2 H [','+ 2 H+ [');
        
        rxnID=sprintf('tRNA_%s_%s_adding_CCA_mitochondrion',aa,codon);
        rxnNames=sprintf('tRNA_%s_%s adding CCA to the tRNA into mitochondrion',aa,codon);
        model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
        
        
        modification_base=sprintf('%s_CCA',base);
        %[rxnID rxnNames eq modification_base]=tRNAreactions(base,C,Reaction,Position);
        
        %    for j=1:numel(eq)
        %        %fprintf(fptr,'\t%s\n',cell2mat(eq(j)));
        %        eq(j)=strrep(eq(j),'[n]',' [mitochondrion]');
        %        eq(j)=strrep(eq(j),'[c]',' [mitochondrion]');
        %        eq(j)=strrep(eq(j),'+ H [','+ H+ [');
        %        eq(j)=strrep(eq(j),'2 H [','2 H+ [');
        %         eq(j)=strrep(eq(j),'5 H [','5 H+ [');
        %        eq(j)=strrep(eq(j),'SAM [','S-adenosyl-L-methionine [');
        %        eq(j)=strrep(eq(j),'SAH [','S-adenosyl-L-homocysteine [');
        %
        %        rxnID(j)=strrep(rxnID(j),'-','_');
        %
        %        %fprintf('%s\n',cell2mat(eq(j)));
        %        model=addYeastReaction(model,cell2mat(eq(j)),rxnID(j),rxnNames(j),0,100,0,{''});
        %
        %    end
        %
        %splicing and amino-acyl
        AA=cell2mat(tRNAsequence(i,13));
        %check for intron
        if cell2mat(tRNAsequence(i,10))==1
            seq=cell2mat(tRNAsequence(i,10));
            from=cell2mat(tRNAsequence(i,12));
            to=cell2mat(tRNAsequence(i,13));
            intron=seq(from:to);
            bc=basecount(intron);
            if bc.A >=1
                eq = sprintf('\t%s [mitochondrion] + NAD [mitochondrion] => %s_splicing [mitochondrion] + diphosphate [mitochondrion] + AMP [mitochondrion] + Arr_GT_p [mitochondrion] + %d ATP [mitochondrion] + %d CTP [mitochondrion] + %d GTP [mitochondrion] + %d UTP [mitochondrion]\n',modification_base,modification_base,(bc.A -1),bc.C,bc.G,bc.T);
            else
                eq = sprintf('\t%s [mitochondrion] + ATP [mitochondrion] + NAD [mitochondrion] => %s_splicing [mitochondrion] + diphosphate [mitochondrion] + AMP [mitochondrion] + Arr_GT_p [mitochondrion] + %d CTP [mitochondrion] + %d GTP [mitochondrion] + %d UTP [mitochondrion]\n',modification_base,modification_base,bc.C,bc.G,bc.T);
            end
            
            rxnID=sprintf('tRNA_%s_%s_splicing',aa,codon);
            rxnNames=sprintf('tRNA_%s_%s splicing',aa,codon);
            model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
            
            if strcmp(RS_gene,'YDR023W;YHR011W')==0
                eq= sprintf('\t%s_splicing%s + ATP[c] + %s%s => AMP[c] + diphosphate[c] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{RS_gene});
                
            else
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_YDR023W',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YDR023W'});
                
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_YHR011W',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YHR011W'});
            end
            

                              
             
            
        else
            
            if strcmp(RS_gene,'YBR121C;YPR081C')==0 && strcmp(RS_gene,'YDR023W;YHR011W')==0
                
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_mitochondrion',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{RS_gene});
                
                                
                
                % fprintf(fptr1,'%s\t%s-tRNA(%s)\n',modification_base,aa,codon);
            elseif strcmp(RS_gene,'YDR023W;YHR011W')==1
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_YDR023W',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YDR023W'});
                
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_YHR011W',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YHR011W'});
                
                
                
            else
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_YBR121C',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YBR121C'});
                
                eq=sprintf('\t%s%s + ATP [mitochondrion] + %s%s => AMP [mitochondrion] + diphosphate [mitochondrion] + %s-tRNA(%s)%s\n',modification_base,cytosol_compartment,AA,cytosol_compartment,aa,antiCodon,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_loading_YPR081C',aa,codon);
                rxnNames=sprintf('%s_%s_loading',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{'YPR081C'});
            end
            
                
               %add dilution of tRNA
                eq=sprintf('\t%s%s => ',modification_base,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_dilution_mitochondrion',aa,codon);
                rxnNames=sprintf('%s_%s_mitochondrion',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''});
                
                %relase a codon
                eq=sprintf('tRNA(%s) [mitochondrion] => %s%s',antiCodon,modification_base,cytosol_compartment);
                eq=strrep(eq,'[c]',' [mitochondrion]');
                rxnID=sprintf('%s_%s_relase_mitochondrion',aa,codon);
                rxnNames=sprintf('%s_%s relase mitochondrion',aa,codon);
                model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,100,0,{''}); 
            
            
        end
    end
    
end
fclose(fptr);
fclose(fptr1);