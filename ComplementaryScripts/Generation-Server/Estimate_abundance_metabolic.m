function  [ ModelExpression]= Estimate_abundance_metabolic(model,growth_rate, complexFileName,outputLPFile, SoplexVersion)

[a b subunit]=xlsread('TableS1.xlsx','Subunit');
[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');


%read solution file
if strcmp(SoplexVersion,'Soplex-3.0')==1
  sol=readSoplex3_results(outputLPFile,model);
else
sol=readSoplex_results(outputLPFile,model);
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
   
   MW=0;
   
   %complex is dimer
   peptide=regexp(complex_list{i},':','split');
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
   for m1=1:numel(peptide)
       peptide_index=find(ismember(Proteins_MW(:,1),strrep(peptide(m1),'-','')));
       MW = MW + ns(m1) * cell2mat(Proteins_MW(peptide_index,2));
   end
   
   rxns=sprintf('%s_complex_dilution_%s',strrep(strrep(complex,':','_'),'-',''),strrep(cell2mat(compartment(i)),'_',''));
   rxn_index = find(ismember(model.rxns,rxns));
   ModelExpression(i)=sol.X(rxn_index)/growth_rate*MW/1000/proteinRatio(growth_rate);% the unit is g/g protein
   
   peptide=regexp(complex,':','split');
   
end

fptr_mass =fopen(sprintf('%s.metabolic_mass',outputLPFile),'w');
for i=1:numel(complex_list)
    fprintf(fptr_mass,'%s\t%s\t%s\t%f\t%f\n',cell2mat(Sfactor.s(i)),cell2mat(Sfactor.complex_list(i)),cell2mat(compartment(i)),ModelExpression(i),1);
end

fclose(fptr_mass);