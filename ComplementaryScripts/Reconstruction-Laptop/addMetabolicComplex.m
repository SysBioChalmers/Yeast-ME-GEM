function model=addMetabolicComplex(model,complex_list,compartment,type_list)
[a b subunit]=xlsread('TableS1.xlsx','Subunit');
complex=unique(complex_list);
k1=1;
reaction_complex={};
mitoProteins={'Q0045','Q0080','Q0085','Q0105','Q0130','Q0250','Q0275'};
for i= 1:numel(complex_list)
             peptide=regexp(complex_list{i},':','split');
             nsubunit={''};
            for jj=1:numel(peptide)
                nsubunit{jj}='';
            end
            if strcmp(peptide,'')==0
                if numel(peptide)==1
                    I=find(ismember(subunit(:,1),peptide));
                    if numel(I)==1
                      nsubunit{1}=strrep(cell2mat(subunit(I,2)),'A','');
                      nsubunit{1}=sprintf('%s ',nsubunit{1});
                      if strcmp(nsubunit{1},'1 ')==1
                          nsubunit{1}='';
                      end
                    end
                      
                else
                    I=find(ismember(subunit(:,1),complex_list{i}));
                    if numel(I)==1
                      subunits=regexp(cell2mat(subunit(I,2)),',','split');
                      nsubunit={''};
                      for jj=1:numel(subunits)
                          nsubunit{jj}=sprintf('%s ',subunits{jj});
                      if strcmp(nsubunit{jj},'1 ')==1
                          nsubunit{jj}='';
                      end
                      end
                      
                    end
                    
                end
            eq= sprintf('%s%s_folding [%s]',nsubunit{1},cell2mat(peptide(1)),cell2mat(compartment(i)));
            if strcmp(compartment(i),'cytoplasm')==0
                %
                if strcmp(compartment(i),'mitochondrion')==1
                        index_mim=find(ismember(mitoProteins,cell2mat(peptide(1))));
                        if numel(index_mim)==0
                            r=sprintf('%s_peptide [cytoplasm] => %s_peptide [mitochondrial membrane]',cell2mat(peptide(1)),cell2mat(peptide(1)));
                            rxnID=sprintf('%s_importing_%s_IMS',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('Importing %s into IMS %s',cell2mat(peptide(1)),cell2mat(compartment(i)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            r=sprintf('%s_peptide [mitochondrial membrane] => %s_peptide [mitochondrion]',cell2mat(peptide(1)),cell2mat(peptide(1)));
                            rxnID=sprintf('%s_importing_%s_Matrix',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('Importing %s into Matrix %s',cell2mat(peptide(1)),cell2mat(compartment(i)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            
                        end
                else
                % other compartments
                r=sprintf('%s_peptide [cytoplasm] => %s_peptide [%s]',cell2mat(peptide(1)),cell2mat(peptide(1)),cell2mat(compartment(i)));
                rxnID=sprintf('%s_importing_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('Importing %s into %s',cell2mat(peptide(1)),cell2mat(compartment(i)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                end

                
                %
                r=sprintf('%s_peptide [%s] => %s_folding [%s]',cell2mat(peptide(1)),cell2mat(compartment(i)),cell2mat(peptide(1)),cell2mat(compartment(i)));
                rxnID=sprintf('%s_folding_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s folding in %s',cell2mat(peptide(1)),cell2mat(compartment(i)));

                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                %misfolding
                r=sprintf('%s_folding [%s] => %s_misfolding [%s]',cell2mat(peptide(1)),cell2mat(compartment(i)),cell2mat(peptide(1)),cell2mat(compartment(i)));
                rxnID=sprintf('%s_misfolding_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

                
                %refolding
                r=sprintf('%s_misfolding [%s] => %s_folding [%s]',cell2mat(peptide(1)),cell2mat(compartment(i)),cell2mat(peptide(1)),cell2mat(compartment(i)));
                rxnID=sprintf('%s_refolding_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s refolding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                %misfolding degradation
                r=sprintf('%s_misfolding [%s] => %s_subunit [%s]',cell2mat(peptide(1)),cell2mat(compartment(i)),cell2mat(peptide(1)),cell2mat(compartment(i)));
                rxnID=sprintf('%s_degradation_misfolding_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                %misfolding dilution
                r=sprintf('%s_misfolding [%s] => ',cell2mat(peptide(1)),cell2mat(compartment(i)));
                rxnID=sprintf('%s_dilution_misfolding_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding dilution',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                
                 index_mim=find(ismember(mitoProteins,cell2mat(peptide(1))));
                 if numel(index_mim)==0
                     r=sprintf('%s_subunit [%s] => %s_subunit [cytoplasm]',cell2mat(peptide(1)),cell2mat(compartment(i)),cell2mat(peptide(1)));
                     rxnID=sprintf('%s_subunit_export_%s',cell2mat(peptide(1)),strrep(cell2mat(compartment(i)),' ','_'));
                     rxnID=strrep(rxnID,'-','');
                     rxnName=sprintf('Exporting %s from %s',cell2mat(peptide(1)),cell2mat(compartment(i)));
                       model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                 end
            else % complex is in the cytoplasm
                %folding
                r=sprintf('%s_peptide [cytoplasm] => %s_folding [cytoplasm]',cell2mat(peptide(1)),cell2mat(peptide(1)));
                rxnID=sprintf('%s_folding_cytoplasm',cell2mat(peptide(1)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s folding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                %misfolding
                r=sprintf('%s_folding [cytoplasm] => %s_misfolding [cytoplasm]',cell2mat(peptide(1)),cell2mat(peptide(1)));
                rxnID=sprintf('%s_misfolding_cytoplasm',cell2mat(peptide(1)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

                %refolding
                r=sprintf('%s_misfolding [cytoplasm] => %s_folding [cytoplasm]',cell2mat(peptide(1)),cell2mat(peptide(1)));
                rxnID=sprintf('%s_refolding_cytoplasm',cell2mat(peptide(1)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s refolding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                %misfolding degradation
                r=sprintf('%s_misfolding [cytoplasm] => %s_subunit [cytoplasm]',cell2mat(peptide(1)),cell2mat(peptide(1)));
                rxnID=sprintf('%s_degradation_misfolding_cytoplasm',cell2mat(peptide(1)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding degradation',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

                %misfolding dilution
                r=sprintf('%s_misfolding [cytoplasm] => ',cell2mat(peptide(1)));
                rxnID=sprintf('%s_dilution_misfolding_cytoplasm',cell2mat(peptide(1)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding',cell2mat(peptide(1)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

                
            end
            
            for k=2:numel(peptide)
                if strcmp(compartment(i),'cytoplasm')==0
                    %
                    if strcmp(compartment(i),'mitochondrion')==1
                        index_mim=find(ismember(mitoProteins,cell2mat(peptide(k))));
                        if numel(index_mim)==0
                            r=sprintf('%s_peptide [cytoplasm] => %s_peptide [mitochondrial membrane]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                            rxnID=sprintf('%s_importing_%s_IMS',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('Importing %s into IMS %s',cell2mat(peptide(k)),cell2mat(compartment(i)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            r=sprintf('%s_peptide [mitochondrial membrane] => %s_peptide [mitochondrion]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                            rxnID=sprintf('%s_importing_%s_Matrix',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('Importing %s into Matrix %s',cell2mat(peptide(k)),cell2mat(compartment(i)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            %misfolding
                            r=sprintf('%s_folding [mitochondrion] => %s_misfolding [mitochondrion]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                            rxnID=sprintf('%s_misfolding_mitochondrion',cell2mat(peptide(k)));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('%s Misfolding',cell2mat(peptide(k)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            %refolding
                            r=sprintf('%s_misfolding [mitochondrion] => %s_folding [mitochondrion]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                            rxnID=sprintf('%s_refolding_mitochondrion',cell2mat(peptide(k)));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('%s refolding',cell2mat(peptide(k)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            %misfolding degradation
                            r=sprintf('%s_misfolding [mitochondrion] => %s_subunit [mitochondrion]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                            rxnID=sprintf('%s_degradation_misfolding_mitochondrion',cell2mat(peptide(k)));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('%s Misfolding degradation',cell2mat(peptide(k)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                            % misfolding dilution
                            r=sprintf('%s_misfolding [mitochondrion] => ',cell2mat(peptide(k)));
                            rxnID=sprintf('%s_dilution_misfolding_mitochondrion',cell2mat(peptide(k)));
                            rxnID=strrep(rxnID,'-','');
                            rxnName=sprintf('%s Misfolding dilution',cell2mat(peptide(k)));
                            model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                            
                        end
                    else
                        %
                        r=sprintf('%s_peptide [cytoplasm] => %s_peptide [%s]',cell2mat(peptide(k)),cell2mat(peptide(k)),cell2mat(compartment(i)));
                        rxnID=sprintf('%s_importing_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                        rxnID=strrep(rxnID,'-','');
                        rxnName=sprintf('Importing %s into %s',cell2mat(peptide(k)),cell2mat(compartment(i)));
                        model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                        %
                        
                    end
                    
                    r=sprintf('%s_peptide [%s] => %s_folding [%s]',cell2mat(peptide(k)),cell2mat(compartment(i)),cell2mat(peptide(k)),cell2mat(compartment(i)));
                    rxnID=sprintf('%s_folding_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                    rxnID=strrep(rxnID,'-','');
                    rxnName=sprintf('%s folding into %s',cell2mat(peptide(k)),cell2mat(compartment(i)));
                    model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                    
                    %misfolding
                        r=sprintf('%s_folding [%s] => %s_misfolding [%s]',cell2mat(peptide(k)),cell2mat(compartment(i)),cell2mat(peptide(k)),cell2mat(compartment(i)));
                        rxnID=sprintf('%s_misfolding_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                        rxnID=strrep(rxnID,'-','');
                        rxnName=sprintf('%s Misfolding',cell2mat(peptide(k)));
                        model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                 
                     %refolding
                        r=sprintf('%s_misfolding [%s] => %s_folding [%s]',cell2mat(peptide(k)),cell2mat(compartment(i)),cell2mat(peptide(k)),cell2mat(compartment(i)));
                        rxnID=sprintf('%s_refolding_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                        rxnID=strrep(rxnID,'-','');
                        rxnName=sprintf('%s refolding',cell2mat(peptide(k)));
                        model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                        
                        
                        %misfoldeing degradation
                        r=sprintf('%s_misfolding [%s] => %s_subunit [%s]',cell2mat(peptide(k)),cell2mat(compartment(i)),cell2mat(peptide(k)),cell2mat(compartment(i)));
                        rxnID=sprintf('%s_degradation_misfolding_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                        rxnID=strrep(rxnID,'-','');
                        rxnName=sprintf('%s Misfolding',cell2mat(peptide(k)));
                        model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                    
                    %misfolding dilution
                        r=sprintf('%s_misfolding [%s] => ',cell2mat(peptide(k)),cell2mat(compartment(i)));
                        rxnID=sprintf('%s_dilution_misfolding_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                        rxnID=strrep(rxnID,'-','');
                        rxnName=sprintf('%s Misfolding dilution',cell2mat(peptide(k)));
                        model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});

                    index_mim=find(ismember(mitoProteins,cell2mat(peptide(k))));
                    if numel(index_mim)==0
                        r=sprintf('%s_subunit [%s] => %s_subunit [cytoplasm]',cell2mat(peptide(k)), cell2mat(compartment(i)), cell2mat(peptide(k)));
                        rxnID=sprintf('%s_subunit_export_%s',cell2mat(peptide(k)),strrep(cell2mat(compartment(i)),' ','_'));
                        rxnID=strrep(rxnID,'-','');
                        rxnName=sprintf('export %s subunit from %s',cell2mat(peptide(k)),cell2mat(compartment(i)));
                        
                        model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                    end
           else % complex is in the cytoplasm
                %folding
                r=sprintf('%s_peptide [cytoplasm] => %s_folding [cytoplasm]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                rxnID=sprintf('%s_folding_cytoplasm',cell2mat(peptide(k)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s folding',cell2mat(peptide(k)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
               %misfolding
                r=sprintf('%s_folding [cytoplasm] => %s_misfolding [cytoplasm]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                rxnID=sprintf('%s_misfolding_cytoplasm',cell2mat(peptide(k)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding',cell2mat(peptide(k)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                %misfolding
                r=sprintf('%s_misfolding [cytoplasm] => %s_folding [cytoplasm]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                rxnID=sprintf('%s_refolding_cytoplasm',cell2mat(peptide(k)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s refolding',cell2mat(peptide(k)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                
                %misfolding degradation
                r=sprintf('%s_misfolding [cytoplasm] => %s_subunit [cytoplasm]',cell2mat(peptide(k)),cell2mat(peptide(k)));
                rxnID=sprintf('%s_degradation_misfolding_cytoplasm',cell2mat(peptide(k)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s Misfolding degradation',cell2mat(peptide(k)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                %misfolding dilution
                r=sprintf('%s_misfolding [cytoplasm] => ',cell2mat(peptide(k)));
                rxnID=sprintf('%s_dilution_misfolding_cytoplasm',cell2mat(peptide(k)));
                rxnID=strrep(rxnID,'-','');
                rxnName=sprintf('%s misfolding dilution',cell2mat(peptide(k)));
                model=addYeastReaction_check(model,r,{rxnID},{rxnName},0,1000,0,{''});
                end
                
                eq=sprintf('%s + %s%s_folding [%s]',eq,nsubunit{k},cell2mat(peptide(k)),cell2mat(compartment(i)));
            end
            
             c=cell2mat(strrep(complex_list(i),':','_'));
             complexName= eq;
             
             eq1=strrep(eq,'folding','subunit');
             
             eq = sprintf('%s => %s_complex [%s]',eq,c,cell2mat(compartment(i)));

             
              rxnID=sprintf('%s_complex_formation_%s',c,strrep(cell2mat(compartment(i)),' ',''));
              rxnID=strrep(rxnID,'-','');
              rxnName=sprintf('%s complex formation from its foldings',c);
              model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
                   
              % complex dilution
              eq=sprintf('%s_complex [%s] => ',c,cell2mat(compartment(i)));
              rxnID=sprintf('%s_complex_dilution_%s',c,strrep(cell2mat(compartment(i)),' ',''));
              rxnID=strrep(rxnID,'-','');
              rxnName=sprintf('%s complex dilution',c);
              model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});
              
              %complex degradation
              eq=sprintf('%s_complex [%s] => %s',c,cell2mat(compartment(i)),eq1);

              rxnID=sprintf('%s_complex_degradation_%s',c,strrep(cell2mat(compartment(i)),' ',''));
              rxnID=strrep(rxnID,'-','');
              rxnName=sprintf('%s complex degradation_%s',c,cell2mat(compartment(i)));
              model=addYeastReaction_check(model,eq,{rxnID},{rxnName},0,1000,0,{''});


            else
                %Reactions that are without complex
%              reaction_complex(k1,1)={'NO Complex'};
%              k1=k1+1;
            end
end
 