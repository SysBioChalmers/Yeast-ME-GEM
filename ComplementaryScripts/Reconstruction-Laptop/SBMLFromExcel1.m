function SBMLFromExcel1(fileName, outputFileName,COBRAFormat,printWarnings)
% SBMLFromExcel
%   Converts a model in the Excel format to SBML.
%
%   fileName        the Excel file
%   outputFileName  the SBML file
%   COBRAFormat     true if the model should be saved in COBRA Toolbox
%                   format. Only limited support at the moment (opt,
%                   default false)
%   printWarnings   true if warnings about model issues should be reported
%                   (opt, default true)
%
%   For a detailed description of the file format, see the supplied manual.
%
%   Usage: SBMLFromExcel(fileName,outputFileName,COBRAFormat,printWarnings)
%
%   Rasmus Agren, 2012-05-24
%

if nargin<3
    COBRAFormat=false;
end
if nargin<4
    printWarnings=true;
end

try
    [crap,crap,raw]=xlsread(fileName,'MODEL');
catch
    throw(MException('','Could not load the MODEL sheet'));
end
raw=cleanImported(raw);

%It is assumed that the first line is labels and that the second one is
%info
allLabels={'MODELID';'MODELNAME';'DEFAULT LOWER';'DEFAULT UPPER';'CONTACT GIVEN NAME';'CONTACT FAMILY NAME';'CONTACT EMAIL';'ORGANIZATION';'TAXONOMY';'NOTES'};
modelAnnotation=[];
modelAnnotation.ID=[];
modelAnnotation.name=[];
modelAnnotation.givenName=[];
modelAnnotation.familyName=[];
modelAnnotation.email=[];
modelAnnotation.organization=[];
modelAnnotation.taxonomy=[];
modelAnnotation.note=[];
defaultLower=-1000;
defaultUpper=1000;

%Loop through the labels
[I J]=ismember(upper(raw(1,:)),allLabels);
I=find(I);
for i=1:numel(I)
    switch J(I(i))
        case 1
            if any(raw{I(i),2})
                modelAnnotation.ID=num2str(raw{2,I(i)}); %Should be string already
                if ~isempty(regexp(modelAnnotation.ID,'[^a-z_A-Z0-9]', 'once'))
                    throw(MException('','Illegal character(s) in model id')); 
                end
            else
                throw(MException('','No model ID supplied'));
            end
        case 2
            if any(raw{2,I(i)})
                modelAnnotation.name=num2str(raw{2,I(i)}); %Should be string already
            else
                throw(MException('','No model name supplied'));
            end
        case 3
            if isnumeric(raw{2,I(i)})
                if ~isnan(raw{2,I(i)})
                    defaultLower=raw{2,I(i)};
                else
                    fprintf('NOTE: DEFAULT LOWER not supplied. Uses -1000.');
                    defaultLower=-1000;
                end
            else
                %Try to convert string to number
                if isnan(str2double(raw{2,I(i)}))
                    throw(MException('','DEFAULT LOWER must be numeric'));
                else
                    defaultLower=raw{2,I(i)};
                end
            end
        case 4
            if isnumeric(raw{2,I(i)})
                if ~isnan(raw{2,I(i)})
                    defaultUpper=raw{2,I(i)};
                else
                    fprintf('NOTE: DEFAULT UPPER not supplied. Uses 1000.');
                    defaultUpper=1000;
                end
            else
                %Try to convert string to number
                if isnan(str2double(raw{2,I(i)}))
                    throw(MException('','DEFAULT UPPER must be numeric'));
                else
                    defaultUpper=raw{2,I(i)};
                end
            end
        case 5
            if any(raw{2,I(i)})
                modelAnnotation.givenName=num2str(raw{2,I(i)}); %Should be string already 
            end
        case 6
            if any(raw{2,I(i)})
                modelAnnotation.familyName=num2str(raw{2,I(i)}); %Should be string already 
            end
        case 7
            if any(raw{2,I(i)})
                modelAnnotation.email=num2str(raw{2,I(i)}); %Should be string already 
            end
        case 8
            if any(raw{2,I(i)})
                modelAnnotation.organization=num2str(raw{2,I(i)}); %Should be string already 
            end    
        case 9
            if any(raw{2,I(i)})
                modelAnnotation.taxonomy=num2str(raw{2,I(i)}); %Should be string already 
            end 
        case 10
            if any(raw{2,I(i)})
                modelAnnotation.note=num2str(raw{2,I(i)}); %Should be string already 
            end     
    end  
end

%Check some needed stuff
if isempty(modelAnnotation.ID)
    throw(MException('','There must be a column named MODELID in the MODEL sheet'));   
end
if isempty(modelAnnotation.name)
    throw(MException('','There must be a column named MODELNAME in the MODEL sheet'));   
end

%Get compartment information
try
    [crap,crap,raw]=xlsread(fileName,'COMPS');
catch
    throw(MException('','Could not load the COMPS sheet'));
end	
raw=cleanImported(raw);

allLabels={'COMPABBREV';'COMPNAME';'INSIDE';'GO TERM'};
compAbbrev={};
compName={};
compOutside={};
compGO={};

%Loop through the labels
[I J]=ismember(upper(raw(1,:)),allLabels);
I=find(I);
for i=1:numel(I)
    switch J(I(i))
        case 1
           compAbbrev=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 2
           compName=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 3
           compOutside=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 4
            if COBRAFormat==false
                compGO=raw(2:end,I(i));
            end
    end
end

%Check that necessary fields are loaded
if isempty(compAbbrev)
    throw(MException('','There must be a column named COMPABBREV in the COMPS sheet'));   
end
if isempty(compName)
    throw(MException('','There must be a column named COMPNAME in the COMPS sheet'));   
end
if isempty(compOutside)
    compOutside=cell(numel(compAbbrev),1);
    compOutside(:)={'NaN'};
end
if isempty(compGO) && COBRAFormat==false
    compGO=cell(numel(compAbbrev),1);
end

%Check that the abbreviated form only contains one character
%Ibrahim Modification
for i=1:length(compAbbrev)
    %if length(compAbbrev{i})==0 || ~all(isstrprop(compAbbrev{i},'alphanum'))
      %  throw(MException('',['The abbreviation of compartment ' compName{i} ' does not follow the form of one alphanumeric']));
    %else
        compAbbrev{i}=lower(compAbbrev{i});
   % end
end

%Check to see that all the OUTSIDE compartments are defined
for i=1:length(compOutside)
   if ~strcmp('NaN',compOutside{i}) %A little weird, but easier
      index=find(strcmp(compOutside{i},compAbbrev),1);
      if isempty(index)
          throw(MException('',['The outside compartment for ' compName{i} ' does not have a corresponding compartment'])); 
      else
          compOutside{i}=int2str(index);
      end
   else
       compOutside{i}=NaN;
   end
end

%Gene info is not needed to load the model, so check for this here
geneNames={};
geneID1={};
geneID2={};
geneShortNames={};
geneCompartments={};
geneKEGG={};
geneComps={};
%Get all the genes and info about them
foundGenes=true;
try
    [crap,crap,raw]=xlsread(fileName,'GENES');
catch
    foundGenes=false;
    if printWarnings==true
        fprintf('WARNING: There is no spreadsheet named GENES\n')
    end
end
if foundGenes==true
    raw=cleanImported(raw);

    allLabels={'GENE NAME';'GENE ID 1';'GENE ID 2';'SHORT NAME';'COMPARTMENT';'KEGG MAPS'};
    
    %Loop through the labels
    [I J]=ismember(upper(raw(1,:)),allLabels);
    I=find(I);
    foundGenes=false;
    for i=1:numel(I)
        switch J(I(i))
            case 1
                geneNames=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
                foundGenes=true;
            case 2
                geneID1=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
            case 3
                geneID2=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
            case 4
                geneShortNames=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
            case 5
                geneCompartments=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
            case 6
                geneKEGG=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        end
    end

    if foundGenes==false
        throw(MException('','There must be a column named GENE NAME in the GENES sheet'));
    end

    if isempty(geneID1)
        geneID1=cell(numel(geneNames));
    end
    if isempty(geneID2)
        geneID2=cell(numel(geneNames));
    end
    if isempty(geneShortNames)
        geneShortNames=cell(numel(geneNames));
    end
    if isempty(geneKEGG)
        geneKEGG=cell(numel(geneNames));
    end
    if isempty(geneCompartments)
       geneCompartments=cell(numel(geneNames),1);
       geneCompartments(:)=compAbbrev(1);
       if printWarnings==true
            fprintf('WARNING: There is no column named COMPARTMENT in the GENES sheet. All genes will be assigned to the first compartment in COMPS. This is merely for annotation and has no effect on the functionality of the model\n');
       end
    end

    %Check that geneName contain only strings and no empty strings
    if ~iscellstr(geneNames)
        throw(MException('','All gene names have to be strings'));
    else
        if any(strcmp('',geneNames)) || any(strcmp('NaN',geneNames))
            throw(MException('','There can be no empty strings in gene names'));
        end
    end

    %Check that geneComp contain only strings and no empty string
    if ~iscellstr(geneCompartments)
        throw(MException('','All gene compartments have to be strings'));
    else
        if ~isempty(find(strcmp('',geneCompartments),1))
            throw(MException('','There can be no empty strings in gene compartments'));
        end
    end

    %Check that all gene compartments correspond to a compartment
    for i=1:length(geneNames)
        index=find(strcmp(geneCompartments{i},compAbbrev));
        if length(index)==1
            geneComps{i}=int2str(index);
        else
            throw(MException('',['The gene ' geneNames{i} ' has a compartment abbreviation that could not be found']));
        end
    end

    %Check that all gene names are unique
    if length(geneNames)~=length(unique(geneNames))
        throw(MException('','Not all gene names are unique'));
    end

    %Check that geneNames contain no weird characters
    illegalCells=regexp(geneNames,'[();:]', 'once'); %Should check for ';' and ':' too
    if ~isempty(cell2mat(illegalCells))
        errorText='Illegal character(s) in gene names:\n';
        for i=1:length(illegalCells)
            if ~isempty(illegalCells{i})
                errorText=[errorText geneNames{i} '\n'];
            end
        end
        throw(MException('',errorText));
    end

    %To fit with other code
    geneID1(strcmp('NaN',geneID1))={NaN};
    geneID2(strcmp('NaN',geneID2))={NaN};
    geneShortNames(strcmp('NaN',geneShortNames))={NaN};
    geneKEGG(strcmp('NaN',geneKEGG))={NaN};
end

%Loads the reaction data
try
    [crap,crap,raw]=xlsread(fileName,'RXNS');
catch
    throw(MException('','Could not load the RXNS sheet'));
end
raw=cleanImported(raw);

allLabels={'RXNID';'NAME';'EQUATION';'EC-NUMBER';'GENE ASSOCIATION';'LOWER BOUND';'UPPER BOUND';'OBJECTIVE';'COMPARTMENT';'SUBSYSTEM';'SBO TERM';'REPLACEMENT ID'};
reactionIDs={};
reactionNames={};
equations={};
ecNumbers={};
geneAssociations={};
lowerBounds=[];
upperBounds=[];
objectives=[];
reactionCompartments={};
reactionSubsystem={};
reactionSBO={};
reactionReplacement={};

%Loop through the labels
[I J]=ismember(upper(raw(1,:)),allLabels);
I=find(I);
for i=1:numel(I)
    switch J(I(i))
        case 1
           reactionIDs=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 2
           reactionNames=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 3
           equations=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 4
           ecNumbers=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 5
           geneAssociations=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 6
           %Check here if they are all numeric
           if all(cellfun(@isnumeric,raw(2:end,I(i))))
                lowerBounds=cell2mat(raw(2:end,I(i)));
           else
                throw(MException('','The lower bounds must be numerical values')); 
           end
        case 7
           %Check here if they are all numeric
           if all(cellfun(@isnumeric,raw(2:end,I(i))))
                upperBounds=cell2mat(raw(2:end,I(i)));
           else
                throw(MException('','The upper bounds must be numerical values')); 
           end
        case 8
           %Check here if they are all numeric
           if all(cellfun(@isnumeric,raw(2:end,I(i))))
                objectives=cell2mat(raw(2:end,I(i)));
           else
                throw(MException('','The objectives must be numerical values')); 
           end
        case 9
        	reactionCompartments=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 10
        	reactionSubsystem=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 11
            if COBRAFormat==false
                reactionSBO=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
            end
        case 12
        	reactionReplacement=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);    
    end
end

%Check that all necessary reaction info has been loaded (RXNID and
%EQUATION)
if isempty(reactionIDs)
     throw(MException('','There must be a column named RXNID in the RXNS sheet'));   
end
if isempty(equations)
     throw(MException('','There must be a column named EQUATION in the RXNS sheet'));   
end

%Check if some other stuff is loaded and populate with default values
%otherwise
if isempty(reactionNames)
   reactionNames=cell(numel(reactionIDs),1);
   if printWarnings==true
        fprintf('WARNING: There is no column named NAME in the RXNS sheet. Empty strings will be used as reaction names\n');
   end
end
if isempty(lowerBounds)
   lowerBounds=nan(numel(reactionIDs),1);
   if printWarnings==true
        fprintf('WARNING: There is no column named LOWER BOUND in the RXNS sheet. Default bounds will be used\n');
   end
end
if isempty(upperBounds)
   upperBounds=nan(numel(reactionIDs),1);
   if printWarnings==true
        fprintf('WARNING: There is no column named UPPER BOUND in the RXNS sheet. Default bounds will be used\n');
   end
end
if isempty(objectives)
   objectives=nan(numel(reactionIDs),1);
   if printWarnings==true
        fprintf('WARNING: There is no column named OBJECTIVE in the RXNS sheet\n');
   end
end
if isempty(reactionCompartments)
   reactionCompartments=cell(numel(reactionIDs),1);
   reactionCompartments(:)=compAbbrev(1);
   if printWarnings==true
        fprintf('WARNING: There is no column named COMPARTMENT in the RXNS sheet. All reactions will be assigned to the first compartment in COMPS. This is merely for annotation and has no effect on the functionality of the model\n');
   end
end

%To fit with other code
reactionNames(strcmp('NaN',reactionNames))={NaN};
ecNumbers(strcmp('NaN',ecNumbers))={NaN};
geneAssociations(strcmp('NaN',geneAssociations))={NaN};
reactionSubsystem(strcmp('NaN',reactionSubsystem))={NaN};
reactionSBO(strcmp('NaN',reactionSBO))={NaN};
reactionReplacement(strcmp('NaN',reactionReplacement))={NaN};

%Check that an SBO-term is associated with each reaction
if COBRAFormat==false
    if ~isempty(reactionSBO)
        if ~iscellstr(reactionSBO)
            reactionSBO=[];
        else
            if ~isempty(find(strcmp('',reactionSBO),1))
                reactionSBO=[];
                if printWarnings==true
                    fprintf('WARNING: Not all reactions have associated SBO-terms. SBO-terms will not be used.\n');
                end
            end
        end
    end
end

if isempty(ecNumbers)
   ecNumbers=cell(numel(reactionIDs),1);
end
if isempty(geneAssociations)
   geneAssociations=cell(numel(reactionIDs),1);
end
if isempty(reactionSubsystem)
   reactionSubsystem=cell(numel(reactionIDs),1);
end

%Replace the reaction IDs for those IDs that have a corresponding 
%replacement name.
I=cellfun(@any,reactionReplacement);
reactionIDs(I)=reactionReplacement(I);

%Check that all reaction IDs are unique
if length(reactionIDs)~=length(unique(reactionIDs))
    for i=1:length(reactionIDs)
    for j=i+1:length(reactionIDs)
        if strmatch(reactionIDs(i),reactionIDs(j),'exact')
            fprintf('The reaction %s is repeated in row %d and %d\n',cell2mat(reactionIDs(i)),i+1,j+1);
        end
    end
end
     throw(MException('','Not all reaction IDs are unique'));
end

%Check that there are no empty strings in reactionIDs or equations
if any(strcmp('NaN',reactionIDs)) || any(strcmp('',reactionIDs))
    throw(MException('','There are empty reaction IDs')); 
end

if any(strcmp('NaN',equations)) || any(strcmp('',equations))
    throw(MException('','There are empty equations')); 
end

%Check that reactionIDs contain no weird characters
illegalCells=regexp(reactionIDs,'[^a-z_A-Z0-9]', 'once');
if ~isempty(cell2mat(illegalCells))
	errorText='Illegal character(s) in reaction IDs:\n';
    for i=1:length(illegalCells)
        if ~isempty(illegalCells{i})
            errorText=[errorText reactionIDs{i} '\n'];
        end
    end
    throw(MException('',errorText));
end

%Check that all reactions have compartments defined
if any(strcmp('NaN',reactionCompartments)) || any(strcmp('',reactionCompartments))
    throw(MException('','All reactions must have an associated compartment string')); 
end

%Fix empty reaction names
I=~cellfun(@any,reactionNames);
reactionNames(I)={''};
    
%Check gene association and compartment for each reaction
for i=1:length(reactionNames)    
    %Check that all gene associations have a match in the gene list
    if ischar(geneAssociations{i}) && length(geneAssociations{i})>0
        indexes=strfind(geneAssociations{i},':'); %Genes are separated by ":" for AND and ";" for OR
        indexes=unique([indexes strfind(geneAssociations{i},';')]);
        if isempty(indexes)
            %See if you have a match (it can't have more than one since the
            %names are unique)
            if isempty(find(strcmp(geneAssociations{i},geneNames),1))
                throw(MException('',['The gene association in reaction ' reactionIDs{i} ' (' geneAssociations{i} ') is not present in the gene list']));
            end   
        else
            temp=[0 indexes numel(geneAssociations{i})+1];
            for j=1:numel(indexes)+1;
                %The reaction has several associated genes
                geneName=geneAssociations{i}(temp(j)+1:temp(j+1)-1);
                if isempty(find(strcmp(geneName,geneNames),1))
                    throw(MException('',['The gene association in reaction ' reactionIDs{i} ' (' geneName ') is not present in the gene list']));
                end
            end
        end
    end
    
    %Check that the compartment for each reaction can be found and save the
    %position in compAbbrev
    index=find(strcmp(reactionCompartments{i},compAbbrev));
    if length(index)==1
        reactionComps{i}=int2str(index);
    else
        throw(MException('',['The reaction ' reactionNames{i} ' has a compartment abbreviation that could not be found'])); 
    end
end

%Check that the reaction names don't contain any forbidden characters
illegalCells=regexp(reactionNames,'[%"<>\\]', 'once');
if ~isempty(cell2mat(illegalCells))
    errorText='Illegal character(s) in reaction names:\n';
    for i=1:length(illegalCells)
        if ~isempty(illegalCells{i})
            errorText=[errorText reactionNames{i} '\n'];
        end
    end
    throw(MException('',errorText));
end

%Get all the metabolites and info about them
try
    [crap,crap,raw]=xlsread(fileName,'METS');
catch
    throw(MException('','Could not load the METS sheet'));
end
raw=cleanImported(raw);			

allLabels={'METID';'METNAME';'UNCONSTRAINED';'MIRIAM';'COMPOSITION';'INCHI';'COMPARTMENT';'REPLACEMENT ID'};

%Load the metabolite information
metaboliteIDs={};
metaboliteNames={};
metConstrained={};
metMiriam={};
metComposition={};
metInchi={};
metCompartments={};
metReplacement={};

%Loop through the labels
[I J]=ismember(upper(raw(1,:)),allLabels);
I=find(I);
for i=1:numel(I)
    switch J(I(i))
        case 1
           metaboliteIDs=raw(2:end,I(i));
        case 2
           metaboliteNames=raw(2:end,I(i)); 
        case 3
           metConstrained=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);
        case 4
           metMiriam=raw(2:end,I(i));  
        case 5
           metComposition=raw(2:end,I(i));
        case 6
           metInchi=raw(2:end,I(i)); 
        case 7
            metCompartments=cellfun(@num2str,raw(2:end,I(i)),'UniformOutput',false);

            %Check that all metabolites have compartments defined
            if any(strcmp('NaN',metCompartments)) || any(strcmp('',metCompartments))
                throw(MException('','All metabolites must have an associated compartment string')); 
            end
        case 8
            metReplacement=raw(2:end,I(i));
    end
end

%Check that necessary fields are loaded (METID)
if isempty(metaboliteIDs)
     throw(MException('','There must be a column named METID in the METS sheet'));   
end

%Check that some other stuff is loaded and use default values otherwise
if isempty(metaboliteNames)
   metaboliteNames=cell(numel(metaboliteIDs),1);
   if printWarnings==true
        fprintf('WARNING: There is no column named METNAME in the METS sheet. METID will be used as name\n');
   end
end
if isempty(metConstrained)
   metConstrained=cell(numel(metaboliteIDs),1);
   metConstrained(:)={'false'};
   if printWarnings==true
        fprintf('WARNING: There is no column named UNCONSTRAINED in the METS sheet. All metabolites will be constrained\n');
   end
end
if isempty(metMiriam)
   metMiriam=cell(numel(metaboliteIDs),1);
end
if isempty(metComposition)
   metComposition=cell(numel(metaboliteIDs),1);
end
if isempty(metInchi)
   metInchi=cell(numel(metaboliteIDs),1);
end
if isempty(metCompartments)
   metCompartments=cell(numel(metaboliteIDs),1);
   metCompartments(:)=compAbbrev(1);
   if printWarnings==true
        fprintf('WARNING: There is no column named COMPARTMENT in the METS sheet. All metabolites will be assigned to the first compartment in COMPS. Note that RAVEN makes extensive use of metabolite names and compartments. Some features will therefore not function correctly if metabolite compartments are not correctly assigned\n');
   end
end

metConstrained(strcmp('NaN',metConstrained))={NaN};
metConstrained(strcmp('1',metConstrained))={'true'};
metConstrained(strcmp('0',metConstrained))={'false'};

[I J]=unique(upper(metaboliteIDs));
if numel(metaboliteIDs)~=numel(I)
    K=1:numel(metaboliteIDs);
    K(J)=[];
    toPrint=metaboliteIDs{K(1)};
    for i=2:numel(K)
        toPrint=[toPrint '\n' metaboliteIDs{K(i)}];
    end
    throw(MException('',['The following metabolites are duplicates:\n' toPrint]));
end

%Replace the metabolite IDs for those IDs that have a corresponding 
%replacement metabolite. This is not used for matching, but will be checked
%for consistency with SBML naming conventions
finalMetIDs=metaboliteIDs; 
I=cellfun(@any,metReplacement);
finalMetIDs(I)=metReplacement(I);

%Check that the metaboliteIDs are strings
if ~iscellstr(finalMetIDs)
    throw(MException('','All metabolite IDs must be strings'));
end

%Check that all metabolites IDs are unique after the mapping as well.
if length(finalMetIDs)~=length(unique(finalMetIDs))
    throw(MException('','Not all metabolite IDs are unique'));
end

%Check that metaboliteIDs contain no weird characters.
illegalCells=regexp(finalMetIDs,'[^a-z_A-Z0-9]', 'once');
if ~isempty(cell2mat(illegalCells))
    errorText='Illegal character(s) in metabolite IDs:\n';
    for i=1:length(illegalCells)
        if ~isempty(illegalCells{i})
            errorText=[errorText finalMetIDs{i} '\n'];
        end
    end
    throw(MException('',errorText));
end

for i=1:length(finalMetIDs)
   %Check that the compartment for each metabolite can be found and save to
   %POSITION IN COMPNAME!!!
   %
   index=find(strcmp(metCompartments{i},compAbbrev));
   if length(index)==1
      metComps{i}=index;
   else
      throw(MException('',['The metabolite "' finalMetIDs{i} '" has a compartment abbreviation that could not be found'])); 
   end

   %Check that the "constrained" fields are "true", "false", or NaN
   if iscellstr(metConstrained(i))
        if ~strcmpi('false',metConstrained(i)) && ~strcmpi('true',metConstrained(i))
            throw(MException('',['The UNCONSTRAINED property for metabolite "' finalMetIDs{i} '" must be "true", "false", or not set']));
        else
            metConstrained{i}=lower(metConstrained{i});
        end
    else
       if ~isnan(metConstrained{i})
           throw(MException('',['The UNCONSTRAINED property for metabolite "' finalMetIDs{i} '" must be "true", "false", or not set']));
       else
           metConstrained{i}='false';
       end
   end
    
   %If the metabolite name isn't set, replace it with the metabolite id
   if ~ischar(metaboliteNames{i}) || isempty(metaboliteNames{i})
       metaboliteNames(i)=finalMetIDs(i);
   end
end

%Check that it doesn't contain any forbidden characters
illegalCells=regexp(metaboliteNames,'["%<>\\]', 'once');
if ~isempty(cell2mat(illegalCells))
    errorText='Illegal character(s) in metabolite names:\n';
    for i=1:length(illegalCells)
        if ~isempty(illegalCells{i})
            errorText=[errorText metaboliteNames{i} '\n'];
        end
    end
    throw(MException('',errorText));
end

%Everything seems fine with the metabolite IDs, compartments, genes, and
%reactions
revIndexes=strfind(equations,' <=> ');
irrevIndexes=strfind(equations,' => ');

if any(cellfun(@isempty,revIndexes) & cellfun(@isempty,irrevIndexes))
    throw(MException('',['The reaction ' reactionIDs{find(cellfun(@isempty,revIndexes) & cellfun(@isempty,irrevIndexes),1)} ' does not have reversibility data']));
end

%Split the reactions in left and right side
lhs=cell(numel(equations),1);
rhs=cell(numel(equations),1);
for i=1:numel(equations)
    stop=[revIndexes{i} irrevIndexes{i}];
    lhs{i}=equations{i}(1:stop(1)-1);
    if isempty(revIndexes{i})
        rhs{i}=equations{i}(stop(1)+4:end);
    else
        rhs{i}=equations{i}(stop(1)+5:end);
    end
end

leftPlus=strfind(lhs,' + ');
rightPlus=strfind(rhs,' + ');
leftSpace=strfind(lhs,' ');
rightSpace=strfind(rhs,' ');

%Preallocate a METxRXNS stoichiometric matrix
S=sparse(length(metaboliteIDs),length(reactionIDs));

%Loop through each of the equations
for i=1:numel(equations)
    for k=-1:2:1 %-1 is for reactant side
        if k==-1
           currentPlus=leftPlus;
           currentSpace=leftSpace;
           currentSide=lhs;
        else
           currentPlus=rightPlus;
           currentSpace=rightSpace;
           currentSide=rhs; 
        end
        
        starts=[1 currentPlus{i}+3 numel(currentSide{i})+4];
        for j=1:numel(starts)-1
            %Check to see if it starts with a coefficient
            coeff=str2double(currentSide{i}(starts(j):currentSpace{i}(find(currentSpace{i}>starts(j),1))-1));
            if ~isnan(coeff)
                %If it starts with a coefficient
                metName=currentSide{i}(currentSpace{i}(find(currentSpace{i}>starts(j),1))+1:starts(j+1)-4);
            else
                coeff=1;
                metName=currentSide{i}(starts(j):starts(j+1)-4);
            end

            %Check to see that it was found in the list
            %Check if the metabolite is present
            metID=find(strcmpi(metName,metaboliteIDs),1);

            %If it didn't find the metabolite
            if isempty(metID)
                throw(MException('',['The metabolite "' metName '" in reaction ' reactionIDs{i} ' was not found in the metabolite list']));
            end

            %If the metabolite was present in more than one copy
            if length(metID)>1
                throw(MException('',['The metabolite "' metName '" in reaction ' reactionIDs{i} ' was found in more than one copy in the metabolite list']));
            end

            %If the metabolites don't match cases
            if strcmp(metName,metaboliteIDs(metID))~=1
                fprintf('WARNING: The metabolite "%s" in reaction %s differs in upper/lower case compared to the metabolite list\n',metName,reactionIDs{i});
            end

            %Check to see that the metabolite isn't already present in
            %the reaction. This means that the reaction is on the form
            %A => A
            if printWarnings==true
                if S(metID,i)~=0
                    fprintf(['WARNING: The reaction ' reactionIDs{i} ' has one or more metabolites both as reactants and as products. Only the net reaction will be exported\n']);
                end
            end

            S(metID,i)=S(metID,i)+coeff*k;
        end
    end
end

%Check that there are no reactions with only reactants or only products
oneSided=find(~(any(S>0) & any(S<0)),1);
if any(oneSided)
    throw(MException('',['The reaction ' reactionIDs{oneSided} ' does not include both reactants and products. If you need such reactions, then make use of unconstrained metabolites instead']));
end

reversibility=zeros(length(equations),1);

reversibility(~cellfun(@isempty,revIndexes))=1;

%Check that there are no conflicting bounds
I=find(lowerBounds>upperBounds);
if any(I)
    errorText='The following reactions(s) have contradicting bounds:\n';
    for i=1:numel(I)
        errorText=[errorText reactionIDs{I(i)} '\n'];
    end
    throw(MException('',errorText));
end

I=find(reversibility==0 & lowerBounds<0);
if any(I)
    errorText='The following reactions(s) have bounds that contradict their directionality:\n';
    for i=1:numel(I)
        errorText=[errorText reactionIDs{I(i)} '\n'];
    end
    throw(MException('',errorText));
end

if printWarnings==true
    %Check that all the metabolites are being used
    involvedMat=S;
    involvedMat(involvedMat~=0)=1;
    usage=sum(involvedMat,2);
    notPresent=find(usage==0);
    unbalanced=find(usage==1);

    if ~isempty(notPresent)
        errorText='WARNING: The following metabolite(s) are never used:\n';
        for i=1:length(notPresent)
            errorText=[errorText '(' finalMetIDs{notPresent(i)} ') ' metaboliteNames{notPresent(i)} '\n'];
        end
        fprintf([errorText '\n']);
    end

    %Note: This should take reactants/products into account. If a
    %metabolite is only a product in all (irreversible) reactions, then is
    %is it unbalanced
    badMets=unbalanced(ismember(metConstrained(unbalanced),'false'));
    if any(badMets)
        errorText='WARNING: The following internal metabolite(s) are only used in one reaction (zero flux is the only solution):\n';
        for i=1:length(badMets)
            errorText=[errorText '(' finalMetIDs{badMets(i)} ' [' compAbbrev{metComps{badMets(i)}} ']) ' metaboliteNames{badMets(i)} '\n'];
        end
        fprintf([errorText '\n']);
    end

    %Check that all metabolites are balanced for C, N, S, P and the number
    %of R-groups
    
    [nC nN nS nP foundComp]=getComposition(metComposition, metInchi);

    %Loop through the reactions and look at those where all the metabolites
    %have composition data. Just count the others
    cantbalanceRxns=[];
    for i=1:size(S,2)
     
        foundMets=find(S(:,i));

        if sum(foundComp(foundMets))==numel(foundMets)
            cBalance=sum(S(foundMets,i).*nC(foundMets));
            nBalance=sum(S(foundMets,i).*nN(foundMets));
            sBalance=sum(S(foundMets,i).*nS(foundMets));
            pBalance=sum(S(foundMets,i).*nP(foundMets));
            %Arbitatary small number
            if abs(cBalance(1,1))>0.00000001
                fprintf('WARNING: The reaction %s is not balanced with respect to carbon\n',reactionIDs{i});
            end
            if abs(nBalance(1,1))>0.00000001
                fprintf('WARNING: The reaction %s is not balanced with respect to nitrogen\n',reactionIDs{i});
            end
            if abs(sBalance(1,1))>0.00000001
                fprintf('WARNING: The reaction %s is not balanced with respect to sulfur\n',reactionIDs{i});
            end
            if abs(pBalance(1,1))>0.00000001
                fprintf('WARNING: The reaction %s is not balanced with respect to phosphorus\n',reactionIDs{i});
            end
        else
            cantbalanceRxns=[cantbalanceRxns i];
        end
    end
    if any(cantbalanceRxns)>0
        fprintf('\nWARNING: The following %s reactions could not be checked for mass balancing\n',num2str(numel(cantbalanceRxns)));
        fprintf('%s\n',reactionIDs{cantbalanceRxns});
    end
end

%All the information has been collected. Time to write SBML!
exportSBML(outputFileName,modelAnnotation,finalMetIDs,metaboliteNames,...
metMiriam, metComposition, metInchi, metConstrained, metComps, reactionIDs, reactionNames, reactionComps, S, reversibility,...
compName, compOutside, compAbbrev, compGO, lowerBounds, upperBounds, objectives, ecNumbers, geneAssociations, reactionSubsystem,...
reactionSBO, defaultLower,defaultUpper, geneNames, geneID1, geneID2, geneShortNames, geneComps, geneKEGG, COBRAFormat);
end

%Cleans up the structure that is imported from using xlsread
function raw=cleanImported(raw)
    %First check that it's a cell array. If a sheet is completely empty,
    %then raw=NaN
    if iscell(raw)
        %Remove columns that aren't strings. If you cut and paste a lot in the sheet 
        %there tends to be columns that are NaN
        I=cellfun(@isstr,raw(1,:));
        raw=raw(:,I);

        %Find the lines that are not commented
        keepers=[1; find(strcmp('',raw(:,1)) | cellfun(@wrapperNAN,raw(:,1)))];
        raw=raw(keepers,:);

        %Check if there are any rows that are all NaN. This could happen if
        %xlsread reads too far. Remove any such rows.
        nans=cellfun(@wrapperNAN,raw);
        I=all(nans,2);
        raw(I,:)=[];

        %Also check if there are any lines that contain only NaNs or white
        %spaces. This could happen if you accidentaly inserted a space
        %somewhere
        whites=cellfun(@wrapperWS,raw);
        I=all(whites,2);
        raw(I,:)=[];
    else
        raw={''};
    end
    
    %Checks if something is NaN. Can't use isnan with cellfun as it does it
    %character by character for strings
    function I=wrapperNAN(A)
       I=any(isnan(A)); 
    end
    
    %Checks if something is all white spaces or NaN
    function I=wrapperWS(A)
        if isnan(A)
            I=true;
        else
            %isstrprob gives an error if boolean
            if islogical(A)
                I=false;
            else
                I=all(isstrprop(A,'wspace'));
            end
        end
    end
end

function [C N S P foundComp]=getComposition(metComposition, metInchi)
    %Assume that they are of the same length
    C=zeros(size(metComposition));
    N=zeros(size(metComposition));
    S=zeros(size(metComposition));
    P=zeros(size(metComposition));
    foundComp=zeros(size(metComposition));
    
    for i=1:length(metComposition)
        inchiError=0;
        if ischar(metInchi{i})
            if length(metInchi{i})>0
                %Find the formula in the Inchi string. Assume that it is
                %everything between the first and second "\"
                indexes=strfind(metInchi{i},'/');
                if length(indexes)==0
                   inchiError=1; 
                else
                   %For some simple molecules such as salts only the first "\" is present
                   if length(indexes)==1
                        formula=metInchi{i}(indexes(1)+1:length(metInchi{i}));
                   else
                        formula=metInchi{i}(indexes(1)+1:indexes(2)-1);
                   end
                   [nC, nN, nS, nP, errorFlag] = compFromFormula(formula);
                   if errorFlag==0
                        C(i)=nC;
                        N(i)=nN;
                        S(i)=nS;
                        P(i)=nP;
                        foundComp(i)=1;
                   end
                end
            else
                inchiError=1;
            end
        else
            inchiError=1;
        end          
            
        if inchiError==1
            %If no InChI could be found
            if ischar(metComposition{i})
                if length(metComposition{i})>0
                    [nC, nN, nS, nP, errorFlag] = compFromFormula(metComposition{i});
                    if errorFlag==0
                        C(i)=nC;
                        N(i)=nN;
                        S(i)=nS;
                        P(i)=nP;
                        foundComp(i)=1;
                    end
                end
            end
        end
    end
end

function [nC, nN, nS, nP, errorFlag] = compFromFormula(formula)
    %IMPORTANT! This does not work if elements can have more than one
    %character. Look in to that!
    errorFlag=0;
    
    % simple modification for Co and Fe we remove it from formula
    % modified by Ibrahim El-Semman
    formula=regexprep(formula,'Co','');
    formula=regexprep(formula,'Fe','');
    %
    %
    
   
    nonNumeric=regexp(formula,'[^0-9]');
    abbrevs=['C';'N';'S';'P'];
    comp=zeros(size(abbrevs));
    for i=1:size(abbrevs)
        index=strfind(formula,abbrevs(i));
        if length(index)==1
           nextNN=find(nonNumeric>index);
           if ~isempty(nextNN)
                comp(i)=str2double(formula(index+1:nonNumeric(nextNN(1))-1));
           else
                comp(i)=str2double(formula(index+1:length(formula)));
           end
           %Might be temporary!! Assumes that there is 1 atom if the
           %str2double thing didn't work
           if isnan(comp(i))
               comp(i)=1;
           end
        else
           if length(index)>1
               %THIS HAS TO BE FIXED! THIS CAN HAPPEN FOR POLYMERS FOR
               %EXAMPLE!!
                errorFlag=1;
           end
        end
    end
    nC=comp(1);
    nN=comp(2);
    nS=comp(3);
    nP=comp(4);
end
    
function [errorFlag] = exportSBML(outputFile,modelAnnotation,metaboliteIDs,...
    metaboliteNames, metMiriam, metComposition, metInchi, metConstrained, metComps,reactionIDs, reactionNames,...
    reactionComps, stochiometricMatrix, reversibility, compName, compOutside, compAbbrev, compGO, lowerBounds,...
    upperBounds, objectives, ecNumbers, geneAssociations, reactionSubsystem, reactionSBO, defaultLower,defaultUpper,...
    geneNames, geneID1, geneID2, geneShortNames, geneComps, geneKEGG, COBRAFormat)

%Generate temporary name
tempFile=tempname;

fid = fopen(tempFile,'w');
if COBRAFormat==false
	intro=['<?xml version="1.0" encoding="UTF-8" ?>'...
    '<sbml xmlns="http://www.sbml.org/sbml/level2/version3" level="2" version="3">'...
    '<model metaid="metaid_' modelAnnotation.ID '" id="' modelAnnotation.ID '" name="' modelAnnotation.name '">'];
    if any(modelAnnotation.note)
        intro=[intro '<notes><body xmlns="http://www.w3.org/1999/xhtml">' modelAnnotation.note '</body></notes>\n'];
    end
    intro=[intro '<annotation>'...
    '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
    '<rdf:Description rdf:about="#metaid_' modelAnnotation.ID '">'];
    if any(modelAnnotation.givenName) && any(modelAnnotation.familyName)
        intro=[intro '<dc:creator rdf:parseType="Resource"><rdf:Bag><rdf:li rdf:parseType="Resource">'];
        if any(modelAnnotation.givenName) && any(modelAnnotation.familyName)
           intro=[intro '<vCard:N rdf:parseType="Resource"><vCard:Family>' modelAnnotation.familyName '</vCard:Family><vCard:Given>' modelAnnotation.givenName '</vCard:Given></vCard:N>'];
        end
        if any(modelAnnotation.email)
           intro=[intro '<vCard:EMAIL>' modelAnnotation.email '</vCard:EMAIL>'];
        end
        if any(modelAnnotation.organization)
           intro=[intro '<vCard:ORG><vCard:Orgname>' modelAnnotation.organization '</vCard:Orgname></vCard:ORG>'];
        end
        
        intro=[intro '</rdf:li></rdf:Bag></dc:creator>'];
        intro=[intro '<dcterms:created rdf:parseType="Resource">'...
        '<dcterms:W3CDTF>' datestr(now,'yyyy-mm-ddTHH:MM:SSZ') '</dcterms:W3CDTF>'... 
        '</dcterms:created>'...
        '<dcterms:modified rdf:parseType="Resource">'...
        '<dcterms:W3CDTF>' datestr(now,'yyyy-mm-ddTHH:MM:SSZ') '</dcterms:W3CDTF>'... 
        '</dcterms:modified>'];
    end
    if any(modelAnnotation.taxonomy)
           intro=[intro '<bqbiol:is><rdf:Bag><rdf:li rdf:resource="urn:miriam:' modelAnnotation.taxonomy '" /></rdf:Bag></bqbiol:is>'];
    end
    intro=[intro '</rdf:Description></rdf:RDF></annotation>\n'];
else
    intro=['<?xml version="1.0" encoding="UTF-8"?><sbml xmlns="http://www.sbml.org/sbml/level2" level="2" version="1" xmlns:html="http://www.w3.org/1999/xhtml"><model id="' modelAnnotation.ID '" name="' modelAnnotation.name '">'];
    if any(modelAnnotation.note)
        intro=[intro '<notes><body xmlns="http://www.w3.org/1999/xhtml">' modelAnnotation.note '</body></notes>\n'];
    end
end

intro=[intro '<listOfUnitDefinitions>'...
    '<unitDefinition id="mmol_per_gDW_per_hr">'...
    '<listOfUnits>'...
    '<unit kind="mole" scale="-3"/>'...
    '<unit kind="second" multiplier="0.00027778" exponent="-1"/>'...
    '</listOfUnits>'...
    '</unitDefinition>'...
    '</listOfUnitDefinitions>\n'...
    '<listOfCompartments>'];

%Write intro
fprintf(fid,intro);

for i=1:length(compName)
    %Check if it's outside anything
    if ~isnan(compOutside{i})
        if COBRAFormat==false
            append=[' outside="C_' compOutside{i} '" spatialDimensions="3"'];
        else
            append=[' outside="C_' compOutside{i} '" spatialDimensions="3"'];
        end
    else
        append=' spatialDimensions="3"';
    end
    if COBRAFormat==false
        fprintf(fid,['<compartment metaid="metaid_C_' int2str(i) '" id="C_' int2str(i) '" name="' compName{i} '"' append ' sboTerm="SBO:0000290">']);
        
        %Check if there is an associated GO-term
        if isstr(compGO{i}) && length(compGO{i})>0
            compinfo=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms='...
                    '"http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" '...
                    'xmlns:bqmodel="http://biomodels.net/model-qualifiers/"><rdf:Description rdf:about="#metaid_C_' int2str(i) '">'...
                    '<bqbiol:is><rdf:Bag><rdf:li rdf:resource="urn:miriam:obo.go:' compGO{i} '"/></rdf:Bag></bqbiol:is>'...
                    '</rdf:Description></rdf:RDF></annotation></compartment>'];
        else
            compinfo='</compartment>';
        end
        fprintf(fid,compinfo);
    else
        fprintf(fid,['<compartment id="C_' int2str(i) '" name="' compName{i} '"' append '></compartment>']);
    end
end

intro='</listOfCompartments>\n<listOfSpecies>';
fprintf(fid,intro);

%Write metabolites
for i=1:length(metaboliteIDs)
    if COBRAFormat==false
        toprint=['<species metaid="metaid_M_' metaboliteIDs{i} '" id="M_' metaboliteIDs{i} '" name="' metaboliteNames{i} '" compartment="C_' int2str(metComps{i}) '" boundaryCondition="' metConstrained{i} '" sboTerm="SBO:0000299">'];
    else
        %Get the formula for the compound
        if ischar(metComposition{i}) && length(metComposition{i})>0
            if ~isempty(regexp(metComposition{i},'[^a-zA-Z0-9]', 'once'))
                formula='_';
            else
                formula=['_' metComposition{i}];
            end
        else
            %Find the formula in the Inchi string. Assume that it is
            %everything between the first and second "\"
            indexes=strfind(metInchi{i},'/');
            if length(indexes)<2
            	formula='_'; 
            else
                if ~isempty(regexp(metInchi{i}(indexes(1)+1:indexes(2)-1),'[^a-zA-Z0-9]', 'once'))
                    formula='_';
                else
                    formula=['_' metInchi{i}(indexes(1)+1:indexes(2)-1)];
                end
            end
        end
        toprint=['<species id="M_' metaboliteIDs{i} '" name="' metaboliteNames{i} formula '" compartment="C_' int2str(metComps{i}) '" boundaryCondition="' metConstrained{i} '">'];
    end
    %Print some stuff if there is a formula for the compound
    if COBRAFormat==false
        if ischar(metComposition{i})
            if length(metComposition{i})>0
                toprint=[toprint '<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>FORMULA: '  metComposition{i} '</p></body></notes>'];
            end
        end
    end
    %Only print annotations for metabolites with some miriam link. This is because I don't know how "unknown" 
    %metabolites should be presented, and it seems unlikely that you will
    %have InChI without a database link. This might be temporary....
    if COBRAFormat==false
        if ischar(metMiriam{i}) && length(metMiriam{i})>0
            toprint=[toprint '<annotation>'];
            %Print InChI if available
            if ischar(metInchi{i})
                if length(metInchi{i})>0
                    toprint=[toprint '<in:inchi xmlns:in="http://biomodels.net/inchi" metaid="metaid_M_' metaboliteIDs{i} '_inchi">InChI=' metInchi{i} '</in:inchi>']; 
                    isInchi=1;
                else
                    isInchi=0;
                end
            else
               isInchi=0; 
            end
            %Print some more annotation stuff
            toprint=[toprint '<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '...
                'xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
                '<rdf:Description rdf:about="#metaid_M_' metaboliteIDs{i} '">'...
                '<bqbiol:is>'...
                '<rdf:Bag>'];
            %If InchI
            if isInchi==1
                toprint=[toprint '<rdf:li rdf:resource="#metaid_M_' metaboliteIDs{i} '_inchi" />'];
            end
            %Print miriam
            toprint=[toprint '<rdf:li rdf:resource="urn:miriam:' metMiriam{i} '"/>'];

            %Finish up
            toprint=[toprint '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>'];
        end
    end
    toprint=[toprint '</species>\n'];
    fprintf(fid,toprint);
end    

%Add information on all modifiers (that is the genes)
%Loop through to replace empty cells with ''. NOT PRETTY!
for i=1:length(geneAssociations)
   if ~ischar(geneAssociations{i})
       geneAssociations{i}='';
   end
end

if COBRAFormat==false
    %First add all the genes
    for i=1:length(geneNames)
       toprint=['<species metaid="metaid_E_' int2str(i) '" id="E_' int2str(i) '" name="' geneNames{i} '" compartment="C_' num2str(geneComps{i}) '" sboTerm="SBO:0000014">'];
       %Print gene name if present
       if ischar(geneShortNames{i}) && length(geneShortNames{i})>0
            toprint=[toprint '<notes><body xmlns="http://www.w3.org/1999/xhtml"><p>SHORT NAME: '  geneShortNames{i} '</p></body></notes>'];
       end
       
       %Print annotation info if present
       if (ischar(geneKEGG{i}) && length(geneKEGG{i})>1) || (ischar(geneID1{i}) && length(geneID1{i})>1) || (ischar(geneID2{i}) && length(geneID2{i})>1)
            toprint=[toprint '<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dc="http://purl.org/dc/elements/1.1/" '...
                'xmlns:dcterms="http://purl.org/dc/terms/" xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" '...
                'xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
                '<rdf:Description rdf:about="#metaid_E_' int2str(i) '"><bqbiol:is><rdf:Bag>'];
                if ischar(geneID1{i}) && length(geneID1{i})>0
                    toprint=[toprint '<rdf:li rdf:resource="urn:miriam:' geneID1{i} '" />'];
                end
                if ischar(geneID2{i}) && length(geneID2{i})>0
                    toprint=[toprint '<rdf:li rdf:resource="urn:miriam:' geneID2{i} '" />'];
                end
                
                if ischar(geneKEGG{i}) && length(geneKEGG{i})>0
                    %KEGG maps are separated by ":"
                    [crap crap crap crap crap crap keggMaps]=regexp(geneKEGG{i},'[:]');
                    
                    for j=1:length(keggMaps)
                        toprint=[toprint '<rdf:li rdf:resource="urn:miriam:kegg.pathway:' keggMaps{j} '" />'];
                    end
                end
                
                toprint=[toprint '</rdf:Bag></bqbiol:is></rdf:Description></rdf:RDF></annotation>'];
       end
       toprint=[toprint '</species>\n'];
       fprintf(fid,toprint);
    end
    
    %Loop through all reactions and find gene associations which contain
    %":", which means that they are governed by several genes
    uniqueGenes=unique(geneAssociations);
    [crap reactions crap crap crap crap crap]=regexp(uniqueGenes,'[:]');
    reactions=find(cellfun('length',reactions));
    [crap crap crap crap crap crap mods]=regexp(uniqueGenes(reactions),'[;]');
    
    geneComplexes={};
    complexGenes={};
    %Loop through each modifier and add the ones that are complexes
    for i=1:numel(reactions)
        %Check to see if it's a complex
        [crap complex crap crap crap crap comGenes]=regexp(mods{i},'[:]');
        complex=find(cellfun('length',complex));
        for j=1:numel(complex)
        geneComplexes=[geneComplexes;mods{i}{complex(j)}];
        complexGenes=[complexGenes;comGenes(complex(j))];
        end
    end
    
    %Remove duplicate complexes
    [geneComplexes,I]=unique(geneComplexes);
    complexGenes=complexGenes(I);
    
    %The SBO term for the complex is set to be "protein". Might not be
    %correct.
    for i=1:length(geneComplexes)
    	toprint=['<species metaid="metaid_Cx_' int2str(i) '" id="Cx_' int2str(i) '" name="' geneComplexes{i} '" compartment="C_2" sboTerm="SBO:0000297">'];
    	toprint=[toprint '</species>\n'];
    	fprintf(fid,toprint);
    end
end

%Finish metbolites
fprintf(fid,'</listOfSpecies>');

%Add reactions
fprintf(fid,'<listOfReactions>');
for i=1:length(reactionIDs)
%Get reversibility
    reversible='false';
    if reversibility(i)==1
        reversible='true';
    end

    if COBRAFormat==false
        if ~isempty(reactionSBO)
             SBO=[' sboTerm="' reactionSBO{i} '"'];
        else
            SBO='';
        end
        fprintf(fid,['<reaction metaid="metaid_R_' reactionIDs{i} '" id="R_' reactionIDs{i} '" name="' reactionNames{i} '" reversible="' reversible '"' SBO '>']);
    else
        fprintf(fid,['<reaction id="R_' reactionIDs{i} '" name="' reactionNames{i} '" reversible="' reversible '">']);
    end
    if COBRAFormat==false
        if ischar(reactionSubsystem{i}) && length(reactionSubsystem{i})>0
            fprintf(fid,'<notes>');
            fprintf(fid,['<body xmlns="http://www.w3.org/1999/xhtml"><p>SUBSYSTEM: ' reactionSubsystem{i} '</p></body>']);
            fprintf(fid,'</notes>');
        end
    else
        fprintf(fid,'<notes>');
        if ~isnan(geneAssociations{i})
            %In order to adhere to the COBRA standards it should be like
            %this:
            %-If only one gene then no parentheses
            %-If only "and" or only "or" there should only be one set of
            %parentheses
            %-If both "and" and "or", then split on "or". This is not
            %complete, but it's the type of relationship supported by the
            %Excel formulation
            aSign=strfind(geneAssociations{i},':');
            oSign=strfind(geneAssociations{i},';');
            if isempty(aSign) && isempty(oSign)
                geneString=geneAssociations{i};
            else
                if isempty(aSign)
                    geneString=['( ' strrep(geneAssociations{i},';',' or ') ' )'];
                else
                    if isempty(oSign)
                        geneString=['( ' strrep(geneAssociations{i},':',' and ') ' )'];
                    else
                        geneString=['(( ' strrep(geneAssociations{i},';',' ) or ( ') ' ))'];
                        geneString=strrep(geneString,':',' and ');
                    end
                end
            end
            
            %geneString=strrep(geneAssociations{i},':',' and ');
            %geneString=strrep(geneString,';',' or ');
            fprintf(fid,['<html:p>GENE_ASSOCIATION: ' geneString '</html:p>']);
        end
        if ischar(reactionSubsystem{i}) && length(reactionSubsystem{i})>0
            fprintf(fid,['<html:p>SUBSYSTEM: ' reactionSubsystem{i} '</html:p>']);
        end
        if ischar(ecNumbers{i}) && length(ecNumbers{i})>0
            fprintf(fid,['<html:p>PROTEIN_CLASS: ' ecNumbers{i} '</html:p>']);
        end
        fprintf(fid,'</notes>');
    end
    
    if COBRAFormat==false
        if ischar(ecNumbers{i}) && length(ecNumbers{i})>0
            toprint=['<annotation><rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '...
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '...
            'xmlns:vCard="http://www.w3.org/2001/vcard-rdf/3.0#" xmlns:bqbiol="http://biomodels.net/biology-qualifiers/" '...
            'xmlns:bqmodel="http://biomodels.net/model-qualifiers/">'...
            '<rdf:Description rdf:about="#metaid_R_' reactionIDs{i} '">'...
            '<bqbiol:is>'...
            '<rdf:Bag>'...
            '<rdf:li rdf:resource="urn:miriam:ec-code:' ecNumbers{i} '"/>'... 
            '</rdf:Bag>'...
            '</bqbiol:is>'...
            '</rdf:Description>'...
            '</rdf:RDF>'...
            '</annotation>'];
            fprintf(fid,toprint);
        end
    end
    
    fprintf(fid,'<listOfReactants>');

    %The reactants have negative values in the stochiometric matrix
    compounds=stochiometricMatrix(:,i);
    reactants=find(compounds<0);
    products=find(compounds>0);
    for j=1:length(reactants)
        %if COBRAFormat==false
        %	tempmetname=metaboliteIDs{reactants(j)};
        %else
        %    tempmetname=[metaboliteIDs{reactants(j)} '_' compAbbrev{metComps{reactants(j)}}]; 
        %end
        toprint=['<speciesReference species="M_' metaboliteIDs{reactants(j)} '" stoichiometry="' char(num2str(-1*compounds(reactants(j)))) '"/>'];
        fprintf(fid,toprint);
    end  
    fprintf(fid,'</listOfReactants><listOfProducts>');
    for j=1:length(products)
        %if COBRAFormat==false
       % 	tempmetname=metaboliteIDs{products(j)};
       % else
       %     tempmetname=[metaboliteIDs{products(j)} '_' compAbbrev{metComps{products(j)}}]; 
       % end
        toprint=['<speciesReference species="M_' metaboliteIDs{products(j)} '" stoichiometry="' char(num2str(compounds(products(j)))) '"></speciesReference>'];
        fprintf(fid,toprint);
    end  
    fprintf(fid,'</listOfProducts>');
    
    if COBRAFormat==false
        %Print modifiers if available.
        if ischar(geneAssociations{i}) && length(geneAssociations{i})>0
            %Loop through the number of modifiers (isoenzymes, complexes...)
            [crap crap crap crap crap crap mods]=regexp(geneAssociations{i},';');
            fprintf(fid,'<listOfModifiers>');
            for j=1:numel(mods)
                %Is it a complex?
                if isempty(strfind(mods{j},':'))
                    %Find the correct gene
                    index=find(strcmp(mods{j},geneNames),1);
                    %Assumes that it is found since that check should have been made
                    %before
                    fprintf(fid,'<modifierSpeciesReference species="E_%s"/>',num2str(index(1)));
                else
                    index=find(strcmp(mods{j},geneComplexes),1);
                    %Assumes that it is found since that check should have been made
                    %before
                    fprintf(fid,'<modifierSpeciesReference species="Cx_%s"/>',num2str(index(1)));
                end
            end
            fprintf(fid,'</listOfModifiers>');
        end
    end
    
    %Print constraints
    if isnan(upperBounds(i))
        upper=defaultUpper;
    else
        upper=upperBounds(i);
    end
    
    if isnan(lowerBounds(i))
       %Check for reversibility
       if reversibility(i)==1
           lower=defaultLower;
       else
           lower=0;
       end
    else
       lower=lowerBounds(i); 
    end
    
    %Print objectives
    if isnan(objectives(i))
        objective=0;
    else
        objective=objectives(i);
    end
    
    if COBRAFormat==false
        fprintf(fid,'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>FLUX_VALUE</ci></math><listOfParameters>');
        fprintf(fid,['<parameter id="LB_R_' reactionIDs{i} '" name="LOWER_BOUND" value="' sprintf('%15.8f',lower) '" units="mmol_per_gDW_per_hr"/><parameter id="UB_R_' reactionIDs{i} '" name="UPPER_BOUND" value="' sprintf('%15.8f',upper) '" units="mmol_per_gDW_per_hr"/><parameter id="OBJ_R_' reactionIDs{i} '" name="OBJECTIVE_COEFFICIENT" value="' sprintf('%15.8f',objective) '" units="dimensionless"/><parameter id="FLUX_VALUE" value="0.00000000" units="mmol_per_gDW_per_hr"/>']);
    else
        fprintf(fid,'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><apply><ci> LOWER_BOUND </ci><ci> UPPER_BOUND </ci><ci> OBJECTIVE_COEFFICIENT </ci></apply></math><listOfParameters>');
        fprintf(fid,['<parameter id="LOWER_BOUND" value="' sprintf('%15.8f',lower) '"/><parameter id="UPPER_BOUND" value="' sprintf('%15.8f',upper) '"/><parameter id="OBJECTIVE_COEFFICIENT" value="' sprintf('%15.8f',objective) '"/>']);    
    end
    fprintf(fid,'</listOfParameters></kineticLaw>');
    fprintf(fid,'</reaction>\n');
end

%Add reactions for the creation of complexes
if COBRAFormat==false
    for i=1:length(geneComplexes)
        fprintf(fid,['<reaction metaid="metaid_R_' int2str(length(reactionIDs)+i) '" id="R_' int2str(length(reactionIDs)+i) '" name="' strrep(geneComplexes{i},':',', ') '" reversible="false" sboTerm="SBO:0000176">']);
        fprintf(fid,'<listOfReactants>');

        for j=1:length(complexGenes{i,1})
            %Assumes that each gene name can be found
            toprint=['<speciesReference species="E_' int2str(find(strcmp(complexGenes{i,1}{j},geneNames))) '" stoichiometry="1"/>'];
            fprintf(fid,toprint);
        end

        fprintf(fid,'</listOfReactants><listOfProducts>');
        fprintf(fid,['<speciesReference species="Cx_' int2str(i) '" stoichiometry="1"></speciesReference>']);
        fprintf(fid,'</listOfProducts>');
        fprintf(fid,'<kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"><ci>FLUX_VALUE</ci></math><listOfParameters>');
        fprintf(fid,['<parameter id="LB_R_' int2str(length(reactionIDs)+i) '" name="LOWER_BOUND" value="0.00000000" units="mmol_per_gDW_per_hr"/><parameter id="UB_R_' int2str(length(reactionIDs)+i) '" name="UPPER_BOUND" value="' sprintf('%15.8f',defaultUpper) '" units="mmol_per_gDW_per_hr"/><parameter id="OBJ_R_' int2str(length(reactionIDs)+i) '" name="OBJECTIVE_COEFFICIENT" value="0.00000000" units="dimensionless"/><parameter id="FLUX_VALUE" value="0.00000000" units="mmol_per_gDW_per_hr"/>']);
        fprintf(fid,'</listOfParameters></kineticLaw>');
        fprintf(fid,'</reaction>\n');
    end
end
fprintf(fid,'</listOfReactions>');

%Write outro
outro='</model></sbml>';
fprintf(fid,outro);

fclose(fid);

%Replace the target file with the temporary file
delete(outputFile);
movefile(tempFile,outputFile,'f');

errorFlag=0;
end