function saveToExcel(model,filename)
% exportToExcelFormat
%   Exports a model structure to the Excel format or to a set of
%   tab-delimited text files.
%
%   model       a model structure
%   filename    filename of the Excel file. This could also be only a path,
%               in which case the model is exported to a set of tab-delimited
%               text files instead
%
%   The resulting Excel file can be used with SBMLFromExcel to generate a 
%   SBML file and to be used for modeling. If only text files should be
%   generated then they will will be saved in the path of filename
%   under the names excelRxns.txt, excelMets.txt, excelGenes.txt,
%   excelModel.txt, and excelComps.txt.
%
%   NOTE: Reactions and genes have no compatment in the model structure,
%         but must be assigned to a compartment to adhere to SBML
%         standards. Therefore the first compartment in model.comps is
%         used. This is purely annotation and has no effect on the
%         functionality of the model.
%
%   Usage: exportToExcelFormat(model,filename)
%
%   Rasmus Agren, 2010-12-16
%   
%   Modified by Ibrahim El-Semman 27 Dec 2012
%   Print the model As standard model

%The model must have an id and description
if ~isfield(model,'id') || ~isfield(model,'description') 
    throw(MException('','The model must have "id" and "description" fields'));
end
if ~any(model.id)
    throw(MException('','The model must have an "id" field'));
end
if ~any(model.description)
    fprintf('WARNING: The "description" field is empty. Uses the id as description.\n');
    model.description=model.id;
end

%Print a note
fprintf(['NOTE: Reactions and genes have no compartment in the model structure, but must be assigned to a '...
    'compartment to adhere to SBML standards. Therefore the first compartment in model.comps is '...
    'used. This is just for annotation and has no effect on the functionality of the model.\n']);
[filePath A B]=fileparts(filename);

if ~any(filePath)
   filePath=pwd; 
end
%If a filename and not only a path was submitted then store the filename
if any(A) && any(B)
    excelFile=[filePath '\' A B];
else
    excelFile=[];
end

if filePath(end)~='\'
    filePath=[filePath '\'];
end

%If the folder doesn't exist then create it
if ~exist(filePath,'dir')
    mkdir(filePath);
end

%Remove all leading and trailing white-spaces from all names. This should
%be ok, but I'm not 100% sure
model.rxns=strtrim(model.rxns);
model.mets=strtrim(model.mets);
model.metName=strtrim(model.metNames);
model.rxnNames=strtrim(model.rxnNames);
model.compNames=strtrim(model.compNames);

%The rest of the fields are checked further down because I don't want to
%make too many checks for structure names here

model.equations=constructEquations(model,model.rxns,true);

%Open for printing the rxn sheet
rxnFile=fopen(fullfile(filePath,'excelRxns.txt'),'wt');

%Print header
fprintf(rxnFile,'#\tRXNID\tNAME\tEQUATION\tEC-NUMBER\tGENE ASSOCIATION\tLOWER BOUND\tUPPER BOUND\tOBJECTIVE\tCOMPARTMENT\tSUBSYSTEM\tSBO TERM\tREPLACEMENT ID\n');

%Check if it should print EC-codes
if isfield(model,'eccodes');
    model.eccodes=strtrim(model.eccodes);
    printEC=true;
else
    printEC=false;
end

%Check if it should print genes
if isfield(model,'grRules');
    model.grRules=strtrim(model.grRules);    
    printRules=true;
    
    %Also do some parsing here
    rules=model.grRules;
    rules=strrep(rules,'(','');
    rules=strrep(rules,')','');
    rules=strrep(rules,' and ',':');
    rules=strrep(rules,' or ',';');
else
    printRules=false;
end

%Check if it should print subsystems
if isfield(model,'subSystems');
    model.subSystems=strtrim(model.subSystems);
    printSubSystems=true;
else
    printSubSystems=false;
end

%Loop through the reactions
for i=1:numel(model.rxns)
   fprintf(rxnFile,['\t' model.rxns{i} '\t' model.rxnNames{i} '\t' model.equations{i} '\t']);
   
   if printEC==true
        fprintf(rxnFile,[model.eccodes{i} '\t']);
   else
        fprintf(rxnFile,'\t');
   end
   
   if printRules==true
        fprintf(rxnFile,[rules{i} '\t']);
   else
        fprintf(rxnFile,'\t');
   end
   
   %Print bounds and objectives
   
   if model.lb(i)==model.ub(i)
       fprintf(rxnFile,[num2str(model.lb(i)) '\t' num2str(model.ub(i)) '\t']);
   else
       fprintf(rxnFile,'\t\t');
   end
   
   if model.c(i)~=0
        fprintf(rxnFile,[num2str(model.c(i)) '\t' ]);
   else
        fprintf(rxnFile,'\t');
   end
   
   fprintf(rxnFile,[model.comps{1} '\t']);
   
   if printSubSystems==true
        fprintf(rxnFile,[model.subSystems{i} '\t']);
   else
        fprintf(rxnFile,'\t');
   end
   
   %Print SBO-terms and Replacement IDs. This is not implemented yet
   fprintf(rxnFile,'\t\t');
   
   fprintf(rxnFile,'\n');
end

fclose(rxnFile);

%Open for printing the metabolites sheet
metFile=fopen(fullfile(filePath,'excelMets.txt'),'wt');

%Check if it should print unconstrained info
if isfield(model,'unconstrained')
    printUnconstrained=true;
else
    printUnconstrained=false;
end

%Check if it should print miriam info
if isfield(model,'metMiriams')
    printMiriam=true;
else
    printMiriam=false;
end

%Check if it should print formula info
if isfield(model,'metFormulas')
    model.metFormulas=strtrim(model.metFormulas);
    printFormulas=true;
else
    printFormulas=false;
end

%Check if it should print InChi info
if isfield(model,'inchis')
    model.inchis=strtrim(model.inchis);
    printInchis=true;
else
    printInchis=false;
end

%Print header
fprintf(metFile,'#\tMETID\tMETNAME\tUNCONSTRAINED\tMIRIAM\tCOMPOSITION\tInChI\tCOMPARTMENT\tREPLACEMENT ID\n');

%Loop through the metabolites
for i=1:numel(model.mets)
   fprintf(metFile,['\t' model.metNames{i} '[' model.compNames{find(ismember(model.comps,model.metComps(i)))} ']\t' model.metNames{i} '\t']);
   
   if printUnconstrained==true
       if model.unconstrained(i) %strcmp(model.metComps(i),model.comps(numel(model.comps)))
            fprintf(metFile,'true\t');
       else
            fprintf(metFile,'false\t');
       end
   else
       fprintf(metFile,'\t');
   end
   
   if printMiriam==true
       if ~isempty(model.metMiriams{i})
            fprintf(metFile,[strtrim(model.metMiriams{i}) '\t']);
       else
            fprintf(metFile,'\t');
       end
   else
       fprintf(metFile,'\t');
   end
   
   if printFormulas==true
       %Print all fomulas if there is no InChi fiels
       if printInchis==false
           fprintf(metFile,[model.metFormulas{i} '\t']);
       else
           %Check if there is an available InChi. If so don't print
           %anything
           if isempty(model.inchis{i})
               fprintf(metFile,[model.metFormulas{i} '\t']);
           else
               fprintf(metFile,'\t');
           end
       end
   else
       fprintf(metFile,'\t');
   end
   
   if printInchis==true
       fprintf(metFile,[model.inchis{i} '\t']);
   else
       fprintf(metFile,'\t');
   end
   
   fprintf(metFile,[model.metComps{i} '\t']);
   
   %There can be no replacement IDs in the structure, but it has to be
   %something to give working met IDs.
   fprintf(metFile,['m' int2str(i) '\t']);
   
   fprintf(metFile,'\n');
end

fclose(metFile);

if isfield(model,'genes')
    %Open for printing the genes sheet
    geneFile=fopen(fullfile(filePath,'excelGenes.txt'),'wt');

    %Check if it should print miriam structures
    if isfield(model,'geneMiriams')
        printMiriams=true;
    else
        printMiriams=false;
    end

    %Check if it should print short gene names
    if isfield(model,'geneShortNames');
        model.geneShortNames=strtrim(model.geneShortNames);
        printShortNames=true;
    else
        printShortNames=false;
    end

    %Print header
    fprintf(geneFile,'#\tGENE NAME\tGENE ID 1\tGENE ID 2\tSHORT NAME\tCOMPARTMENT\tKEGG MAPS\n');

    %Loop through the genes
    for i=1:numel(model.genes)
        fprintf(geneFile,['\t' model.genes{i} '\t']);

        if printMiriams==true
            %This is a little tricky. Should print those IDs that are not for
            %kegg maps.
            if ~isempty(model.geneMiriams{i})
                nonKegg=1:numel(model.geneMiriams{i}.name);
                kegg=strmatch('kegg.pathway',strtrim(model.geneMiriams{i}.name),'exact');
                nonKegg(kegg)=[];

                if ~isempty(nonKegg)
                    %The number here must be 1 or 2 as the converter is written
                    %now
                    fprintf(geneFile,[strtrim(model.geneMiriams{i}.name{nonKegg(1)}) ':' strtrim(model.geneMiriams{i}.value{nonKegg(1)}) '\t']);

                    if numel(nonKegg)>1
                        fprintf(geneFile,[strtrim(model.geneMiriams{i}.name{nonKegg(2)}) ':' strtrim(model.geneMiriams{i}.value{nonKegg(2)}) '\t']);
                    else
                        fprintf(geneFile,'\t');
                    end
                else
                    %Only kegg maps
                    fprintf(geneFile,'\t\t');
                end
            else
                fprintf(geneFile,'\t\t');
            end
        else
            fprintf(geneFile,'\t\t');
        end

        if printShortNames==true
            fprintf(geneFile,[model.geneShortNames{i} '\t']);
        else
            fprintf(geneFile,'\t');
        end

        %Since genes have no compartment in the model structure
        fprintf(geneFile,[model.comps{1} '\t']);

        if printMiriams==true
            %This is a little tricky. Should print those IDs that are for
            %kegg maps.
            if ~isempty(model.geneMiriams{i})
                kegg=strmatch('kegg.pathway',strtrim(model.geneMiriams{i}.name),'exact');
                if ~isempty(kegg)
                    for j=1:numel(kegg)
                       if j<numel(kegg)
                           pad=':';
                       else
                           pad='';
                       end
                       fprintf(geneFile,[strtrim(model.geneMiriams{i}.value{kegg(j)}) pad]);
                    end
                else
                    %No kegg maps
                    fprintf(geneFile,'\t');
                end
            else
                fprintf(geneFile,'\t');
            end
        else
            fprintf(geneFile,'\t');
        end

        fprintf(geneFile,'\n');
    end
    fclose(geneFile);
end

if isfield(model,'id')
    %Open for printing the model sheet
    modelFile=fopen(fullfile(filePath,'excelModel.txt'),'wt');

    %Print header
    fprintf(geneFile,'#\tMODELID\tMODELNAME\tDEFAULT LOWER\tDEFAULT UPPER\tCONTACT GIVEN NAME\tCONTACT FAMILY NAME\tCONTACT EMAIL\tORGANIZATION\tTAXONOMY\tNOTES\n');
    
    %Print model ID and name. It is assumed that the default lower/upper
    %bound correspond to min/max of the bounds
    fprintf(geneFile,['\t' model.id '\t' model.description '\t' num2str(-1000) '\t' num2str(1000) '\tRasmus\tAgren\trasmus.agren@chalmers.se\tChalmers University of Technology\t9606\t\n']);
    fclose(modelFile);
end

if isfield(model,'comps')
    %Open for printing the model sheet
    compsFile=fopen(fullfile(filePath,'excelComps.txt'),'wt');

    %Print header
    fprintf(compsFile,'#\tCOMPABBREV\tCOMPNAME\tINSIDE\tGO TERM\n');
    
    for i=1:numel(model.comps)
       fprintf(compsFile,['\t' model.comps{i} '\t' model.compNames{i} '\t' model.compOutside{i} '\t\n']);
    end
    fclose(compsFile);
end

%Now it has generated the text files. If a full file name was submitted
%then should those files be merged into one Excel file and then deleted
if any(excelFile)
    foundError=false;
    textFiles={'excelRxns.txt' 'excelMets.txt' 'excelGenes.txt' 'excelComps.txt' 'excelModel.txt'};
    formatStrings={'%s%s%s%s%s%s%n%n%n%s%s%s%s' '%s%s%s%s%s%s%s%s%s' '%s%s%s%s%s%s%s' '%s%s%s%s%s' '%s%s%s%n%n%s%s%s%s%s%s'};
    sheets={'RXNS' 'METS' 'GENES' 'COMPS' 'MODEL'};
    captions={{'#' 'RXNID' 'NAME' 'EQUATION' 'EC-NUMBER' 'GENE ASSOCIATION' 'LOWER BOUND' 'UPPER BOUND' 'OBJECTIVE' 'COMPARTMENT' 'SUBSYSTEM' 'SBO TERM' 'REPLACEMENT ID'};...
            {'#' 'METID' 'METNAME' 'UNCONSTRAINED' 'MIRIAM' 'COMPOSITION' 'InChI' 'COMPARTMENT' 'REPLACEMENT ID'};...
            {'#' 'GENE NAME' 'GENE ID 1' 'GENE ID 2' 'SHORT NAME' 'COMPARTMENT' 'KEGG MAPS'};...
            {'#' 'COMPABBREV' 'COMPNAME' 'INSIDE' 'GO TERM'};...
            {'#' 'MODELID' 'MODELNAME' 'DEFAULT LOWER' 'DEFAULT UPPER' 'CONTACT GIVEN NAME' 'CONTACT FAMILY NAME' 'CONTACT EMAIL' 'ORGANIZATION' 'TAXONOMY' 'NOTES'}};
    for i=1:numel(textFiles)    
        fid=fopen(fullfile(filePath,textFiles{i}),'r');
        C = textscan(fid,formatStrings{i},'Delimiter','\t','Headerlines',1,'Whitespace','');
        fclose(fid);

        %Since xlswrite requires a cell array I construct it here. This might
        %not be the fastest or best way to do this.
        cellArray=[];
        for j=1:numel(C)
            %This is because numbers are read as vectors
            if iscell(C{j})
                %It could be that the cell array contains empty values. Those
                %should be ''
                cellArray=[cellArray C{j}];
            else
                cellArray=[cellArray num2cell(C{j})];
            end
        end
        %I do like this rather than to read the captions from the file because
        %it makes it easier with the doubles
        cellArray=[captions{i};cellArray];
        errorFlag=xlswrite(excelFile,cellArray,sheets{i});
        
        if errorFlag==0
           fprintf('There was an error in writing the Excel file. Keeping the text files in the specified directory');
           foundError=true;
        end
    end
    
    %Delete the text files if the Excel sheet was generated sucessfully
    if foundError==false
        for i=1:numel(textFiles)
           delete(fullfile(filePath,textFiles{i}));
        end
    end
end
end
