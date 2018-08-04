function  Sfactor= Estimate_abundance(model,growth_rate, complexFileName,outputLPFile, ExpressionFileName,total_protein_molecules_model,SoplexVersion,Iteration)

[a b subunit]=xlsread('TableS1.xlsx','Subunit');


%read solution file
if strcmp(SoplexVersion,'Soplex-3.0')==1
  sol=readSoplex3_results(outputLPFile,model);
else
sol=readSoplex_results(outputLPFile,model);
end

%read the expression file

[a b Protein]=xlsread(ExpressionFileName);
   %total protein molecules in the expression data
   Expression_col=4;
   total_protein_molecules_data=0;
   for kk1=2:numel(Protein(:,Expression_col))
       total_protein_molecules_data = total_protein_molecules_data + cell2mat(Protein(kk1,Expression_col));
   end
   
%read complex file
[s complex_list compartment] =textread(complexFileName,'%s%s%s');

%

Sfactor.s=s;
Sfactor.complex_list=complex_list;
Sfactor.compartment=compartment;
expression=zeros(numel(complex_list),1);
ModelExpression=zeros(numel(complex_list),1);

complex_list_Expressed={''};
kk=1;
for i=1:numel(complex_list)
   complex=complex_list{i};
   %set the dilution reaction
   rxns=sprintf('%s_complex_dilution_%s',strrep(strrep(complex,':','_'),'-',''),strrep(cell2mat(compartment(i)),'_',''));
   rxn_index = find(ismember(model.rxns,rxns));
   ModelExpression(i)=sol.X(rxn_index)/growth_rate*(13*6.023e8);
   
   %Counting number of subunits in the complex
      rxns=sprintf('%s_complex_formation_%s',strrep(strrep(complex,':','_'),'-',''),strrep(cell2mat(compartment(i)),'_',''));
      rxn_index = find(ismember(model.rxns,rxns));
      ns_index=find(model.S(:,rxn_index)<0);
      nsubunits= -sum(full(model.S(ns_index,rxn_index)));
   %Compute the ratio of number of molecules in this ratio to the the total
   %number of protein molecules
      ModelExpression(i)=ModelExpression(i)*nsubunits/total_protein_molecules_model;
    
   
   peptide=regexp(complex,':','split');
   %number of subunits
        if numel(peptide)==1
                        ns=1;
                        I=find(ismember(subunit(:,1),peptide));
                        if numel(I)==1
                            nsubunit=strrep(cell2mat(subunit(I,2)),'A','');
                            nsubunit=sprintf('%s ',nsubunit);
                            ns=str2num(nsubunit);
                        end
                    elseif (numel(peptide)>1)  %complex is a hetrodimer.
                        I=find(ismember(subunit(:,1),complex_list{i}));
                        ns=ones(numel(peptide),1);
                        if numel(I)==1
                            subunits=regexp(cell2mat(subunit(I,2)),',','split');
                            for jj=1:numel(subunits)
                               ns(jj)=str2num(subunits{jj});
                             end
                        end
        end
                    
        
   abundance=0;
      k=0;
      for j=1:numel(peptide)
          %find abundance for each subunit
          index=find(ismember(Protein(:,1),peptide{j}));
          if numel(index)>=1
              k=k+1; %number of subunits is found in the data
              abundance=abundance+ cell2mat(Protein(index,Expression_col))/ns(j);
          end
      end
          if k==0 %no subunit is in the data
              %expression(i)=1000000; % all expressed
          else
              expression(i)=abundance/k*nsubunits/total_protein_molecules_data; %the expression is the aveage of the subunits
          end
      
end

Sfactor.expression=expression;
Sfactor.ModelExpression=ModelExpression;

%
if Iteration==1
    S_factor=zeros(numel(complex_list),1);
else
    [a b c express_data express_model S_factor]=textread(strrep(sprintf('%s.sfactor',outputLPFile),sprintf('k_%d',Iteration),sprintf('k_%d',Iteration-1)),'%s%s%s%f%f%f');
end

fptr_sfactor=fopen(sprintf('%s.sfactor',outputLPFile),'w');
for i=1:numel(complex_list)
    
    if ModelExpression(i)==0 ||expression(i)==0
        fprintf(fptr_sfactor,'%s\t%s\t%s\t%f\t%f\t%f\n',cell2mat(Sfactor.s(i)),cell2mat(Sfactor.complex_list(i)),cell2mat(compartment(i)),expression(i),ModelExpression(i),1);
        
    else
        if S_factor(i) ~=1 %save the old saturation factors
            if Iteration==1
              fprintf(fptr_sfactor,'%s\t%s\t%s\t%f\t%f\t%f\n',cell2mat(Sfactor.s(i)),cell2mat(Sfactor.complex_list(i)),cell2mat(compartment(i)),expression(i),ModelExpression(i),ModelExpression(i)/expression(i));

            else
            fprintf(fptr_sfactor,'%s\t%s\t%s\t%f\t%f\t%f\n',cell2mat(Sfactor.s(i)),cell2mat(Sfactor.complex_list(i)),cell2mat(compartment(i)),expression(i),ModelExpression(i),S_factor(i));
            end
        else
           fprintf(fptr_sfactor,'%s\t%s\t%s\t%f\t%f\t%f\n',cell2mat(Sfactor.s(i)),cell2mat(Sfactor.complex_list(i)),cell2mat(compartment(i)),expression(i),ModelExpression(i),ModelExpression(i)/expression(i));

        end
    end
end
fclose(fptr_sfactor);