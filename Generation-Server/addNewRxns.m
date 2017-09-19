function newModel=addNewRxns(model,rxnsToAdd,eqnType,compartment,allowNewMets)
% addRxns
%   Adds reactions to a model
%
%   model            a model structure
%   rxnsToAdd        the reaction structure can have the following fields:
%                    rxns       cell array with unique strings that 
%                               identifies each reaction
%                    equations  cell array with equation strings. Decimal 
%                               coefficients are expressed as ?1.2?. 
%                               Reversibility is indicated by ?<=>? or ?=>?
%                    rxnNames   cell array with the names of each reaction
%                               (opt, default '')
%                    lb         vector with the lower bounds (opt, default
%                               -inf for reversible reactions and 0 for
%                               irreversible)
%                    ub         vector with the upper bounds (opt, default
%                               inf)
%                    c          vector with the objective function
%                               coefficients (opt, default 0)
%                    eccodes    cell array with the EC-numbers for each
%                               reactions. Delimit several EC-numbers with
%                               ";" (opt, default '')
%                    subSystems cell array with the subsystems for each
%                               reaction (opt, default '')
%                    grRules    cell array with the gene-reaction
%                               relationship for each reaction. For example
%                               "(A and B) or (C)" means that the reaction 
%                               could be catalyzed by a complex between 
%                               A & B or by C on its own. All the genes 
%                               have to be present in model.genes. Add 
%                               genes with addGenes before calling this
%                               function if needed (opt, default '')
%   eqnType          double describing how the equation string should be
%                    interpreted
%                    1 - The metabolites are matched to model.mets. New
%                        metabolites (if allowed) are added to
%                        "compartment"
%                    2 - The metabolites are matched to model.metNames and
%                        all metabolites are assigned to "compartment". Any
%                        new metabolites that are added will be assigned
%                        IDs "m1", "m2"... If IDs on the same form are 
%                        already used in the model then the first available 
%                        integers will be used  
%                    3 - The metabolites are written as 
%                        "metNames[compNames]". Only compartments in
%                        model.compNames are allowed. Any
%                        new metabolites that are added will be assigned
%                        IDs "m1", "m2"... If IDs on the same form are 
%                        already used in the model then the first available 
%                        integers will be used  
%   compartment      a string with the compartment the metabolites should
%                    be placed in when using eqnType=2. Must match 
%                    model.compNames (opt when eqnType=1 or eqnType=3) 
%   allowNewMets     true if the function is allowed to add new
%                    metabolites. It is highly recommended to first add
%                    any new metabolites with addMets rather than
%                    automatically through this function. addMets supports
%                    more annotation of metabolites, allows for the use of
%                    exchange metabolites, and using it reduces the risk
%                    of parsing errors (opt, default false)
%                     
%   newModel         an updated model structure
%
%   NOTE: This function does not make extensive checks about formatting of
%   gene-reaction rules. 
%
%   NOTE: When adding metabolites to a compartment where it previously
%   doesn't exist, the function will copy any available information from
%   the metabolite in another compartment.
%
%   Usage: newModel=addRxns(model,rxnsToAdd,eqnType,compartment,allowNewMets)
%
%   Rasmus Agren, 2012-04-02
%

if nargin<4
    compartment=[];
end
if nargin<5
    allowNewMets=false;
end

newModel=model;

%If no reactions should be added
if isempty(rxnsToAdd)
    return;
end

%Check the input
if ~isnumeric(eqnType)
    throw(MException('','eqnType must be numeric'));
else
    if ~ismember(eqnType,[1 2 3])
        throw(MException('','eqnType must be 1, 2, or 3'));
    end
end

if eqnType==2 || (eqnType==1 && allowNewMets==true)
    if ~ischar(compartment)
        throw(MException('','compartment must be a string'));
    end
    if ~ismember(compartment,model.compNames)
        throw(MException('','compartment must match one of the compartments in model.compNames'));
    end
end

if ~isfield(rxnsToAdd,'rxns')
    throw(MException('','rxns is a required field in rxnsToAdd'));
else
    %To fit with some later printing
    rxnsToAdd.rxns=rxnsToAdd.rxns(:);
end
if ~isfield(rxnsToAdd,'equations')
    throw(MException('','equations is a required field in rxnsToAdd'));
end

if ~iscellstr(rxnsToAdd.rxns) && ~ischar(rxnsToAdd.rxns)
    %It could also be a string, but it's not encouraged
    throw(MException('','rxnsToAdd.rxns must be a cell array of strings'));
else
    rxnsToAdd.rxns=cellstr(rxnsToAdd.rxns);
end
if ~iscellstr(rxnsToAdd.equations) && ~ischar(rxnsToAdd.equations)
    %It could also be a string, but it's not encouraged
    throw(MException('','rxnsToAdd.equations must be a cell array of strings'));
else
    rxnsToAdd.equations=cellstr(rxnsToAdd.equations);
end

%Check some formatting
illegalCells=regexp(rxnsToAdd.rxns,'[^a-z_A-Z0-9]', 'once');
if ~isempty(cell2mat(illegalCells))
    errorText='Illegal character(s) in reaction IDs:\n';
    for i=1:length(illegalCells)
        if ~isempty(illegalCells{i})
            errorText=[errorText rxnsToAdd.rxns{i} '\n'];
        end
    end
    throw(MException('',errorText));
end
if isfield(rxnsToAdd,'rxnNames')
    %If the mets are metNames
    illegalCells=regexp(rxnsToAdd.rxnNames,'["%<>\\]', 'once');
    if ~isempty(cell2mat(illegalCells))
        errorText='Illegal character(s) in reaction names:\n';
        for i=1:length(illegalCells)
            if ~isempty(illegalCells{i})
                errorText=[errorText rxnsToAdd.rxnNames{i} '\n'];
            end
        end
        throw(MException('',errorText));
    end
end

nRxns=numel(rxnsToAdd.rxns);
nOldRxns=numel(model.rxns);
filler=cell(nRxns,1);
filler(:)={''};
largeFiller=cell(nOldRxns,1);
largeFiller(:)={''};

%***Add everything to the model except for the equations.
if numel(rxnsToAdd.equations)~=nRxns
   throw(MException('','rxnsToAdd.equations must have the same number of elements as rxnsToAdd.rxns'));
end

%Parse the equations. This is done at this early stage since I need the
%reversibility info
[S mets badRxns reversible]=constructS(rxnsToAdd.equations);
if any(badRxns)
	fprintf('WARNING: The following equations have one or more metabolites both as substrate and product. Only the net equations will be added\n');
    I=find(badRxns);
    for i=1:numel(I)
        fprintf(['\t' rxnsToAdd.rxns{I(i)} '\n']);
    end
end
if newModel.first==1
    newModel.rev=reversible;
    newModel.rxns=rxnsToAdd.rxns(:);
    
else
    newModel.rev=[newModel.rev;reversible];
    newModel.rxns=[newModel.rxns;rxnsToAdd.rxns(:)];
    
    
end

if isfield(rxnsToAdd,'rxnNames')
   if numel(rxnsToAdd.rxnNames)~=nRxns
       throw(MException('','rxnsToAdd.rxnNames must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'rxnNames')
       newModel.rxnNames=largeFiller;
   end
   newModel.rxnNames=[newModel.rxnNames;rxnsToAdd.rxnNames(:)];
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'rxnNames')
       newModel.rxnNames=[newModel.rxnNames;filler];
   end
end

if isfield(rxnsToAdd,'lb')
   if numel(rxnsToAdd.lb)~=nRxns
       throw(MException('','rxnsToAdd.lb must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'lb')
       newModel.lb=zeros(nOldRxns,1);
       newModel.lb(newModel.rev~=0)=-inf;
   end
   if newModel.first==1
       newModel.lb= rxnsToAdd.lb(:);
   else
       newModel.lb=[newModel.lb;rxnsToAdd.lb(:)];
   end
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'lb')
       I=zeros(nRxns,1);
       I(reversible~=0)=-inf;
       if newModel.first==1
          newModel.lb=I;
       else
          newModel.lb=[newModel.lb;I];
       end
   end
end

if isfield(rxnsToAdd,'ub')
   if numel(rxnsToAdd.ub)~=nRxns
       throw(MException('','rxnsToAdd.ub must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'ub')
       newModel.ub=inf(nOldRxns,1);
   end
   if newModel.first==1
       newModel.ub=rxnsToAdd.ub(:);
   else
       newModel.ub=[newModel.ub;rxnsToAdd.ub(:)];
   end
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'ub')
       newModel.ub=[newModel.ub;inf(nRxns,1)];
   end
end

if isfield(rxnsToAdd,'c')
   if numel(rxnsToAdd.c)~=nRxns
       throw(MException('','rxnsToAdd.c must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'c')
       newModel.c=zeros(nOldRxns,1);
   end
   if newModel.first==1
       newModel.c=rxnsToAdd.c(:);
   else
       newModel.c=[newModel.c;rxnsToAdd.c(:)];
   end
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'c')
       newModel.c=[newModel.c;zeros(nRxns,1)];
   end
end

if isfield(rxnsToAdd,'eccodes')
   if numel(rxnsToAdd.eccodes)~=nRxns
       throw(MException('','rxnsToAdd.eccodes must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'eccodes')
       newModel.eccodes=largeFiller;
   end
   if newModel.first==1
       newModel.eccodes=rxnsToAdd.eccodes(:);
   else
       newModel.eccodes=[newModel.eccodes;rxnsToAdd.eccodes(:)];
   end
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'eccodes')
       
       newModel.eccodes=[newModel.eccodes;filler];
   end
end

if isfield(rxnsToAdd,'subSystems')
   if numel(rxnsToAdd.subSystems)~=nRxns
       throw(MException('','rxnsToAdd.subSystems must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'subSystems')
       newModel.subSystems=largeFiller;
   end
   newModel.subSystems=[newModel.subSystems;rxnsToAdd.subSystems(:)];
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'subSystems')
       newModel.subSystems=[newModel.subSystems;filler];
   end
end

if isfield(rxnsToAdd,'grRules')
   if numel(rxnsToAdd.grRules)~=nRxns
       throw(MException('','rxnsToAdd.grRules must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'grRules')
       newModel.grRules=largeFiller;
   end
   newModel.grRules=[newModel.grRules;rxnsToAdd.grRules(:)];
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'grRules')
       newModel.grRules=[newModel.grRules;filler];
   end
end

if isfield(rxnsToAdd,'rxnFrom')
   if numel(rxnsToAdd.rxnFrom)~=nRxns
       throw(MException('','rxnsToAdd.rxnFrom must have the same number of elements as rxnsToAdd.rxns'));
   end
   %Fill with standard if it doesn't exist
   if ~isfield(newModel,'rxnFrom')
       newModel.rxnFrom=largeFiller;
   end
   newModel.rxnFrom=[newModel.rxnFrom;rxnsToAdd.rxnFrom(:)];
else
    %Fill with standard if it doesn't exist
   if isfield(newModel,'rxnFrom')
       newModel.rxnFrom=[newModel.rxnFrom;filler];
   end
end

%Check that ids contain no weird characters. This is only done if the
%equations are with metabolite ids
if eqnType==1
    illegalCells=regexp(mets,'[^a-z_A-Z0-9]', 'once');
    if ~isempty(cell2mat(illegalCells))
        errorText='Illegal character(s) in metabolite IDs:\n';
        for i=1:length(illegalCells)
            if ~isempty(illegalCells{i})
                errorText=[errorText mets{i} '\n'];
            end
        end
        throw(MException('',errorText));
    end
else
    %If the mets are metNames
    illegalCells=regexp(mets,'["%<>\\]', 'once');
    if ~isempty(cell2mat(illegalCells))
        errorText='Illegal character(s) in metabolite names:\n';
        for i=1:length(illegalCells)
            if ~isempty(illegalCells{i})
                errorText=[errorText mets{i} '\n'];
            end
        end
        throw(MException('',errorText));
    end
end

%***Start parsing the equations and adding the info to the S matrix
%The mets are matched to model.mets
if eqnType==1
    [I J]=ismember(mets,model.mets);
    if ~all(I)
        if allowNewMets==true
            %Add the new mets
            metsToAdd.mets=mets(~I);
            metsToAdd.metNames=metsToAdd.mets;
            metsToAdd.compartments=cell(numel(metsToAdd.mets),1);
            metsToAdd.compartments(:)={compartment};
            newModel=addMets(newModel,metsToAdd);
        else
            throw(MException('','One or more equations contain metabolites that are not in model.mets. Set allowNewMets to true to allow this function to add metabolites or use addMets to add them before calling this function'));
        end
    end
    %Calculate the indexes of the metabolites and add the info
    metIndexes=J;
    metIndexes(~I)=numel(newModel.mets)-sum(~I)+1:numel(newModel.mets);
end

%Do some stuff that is the same for eqnType=2 and eqnType=3
if eqnType==2 || eqnType==3
    %First find the first available indexes on
    %the form "m1", "m2"..
    maxCurrent=ceil(max(cellfun(@getInteger,model.mets)));
    
    if length(maxCurrent)==0
          maxCurrent=numel(model.mets);
    end
    %For later..
    [crap I]=ismember(model.metComps,model.comps);
    t2=strcat(model.metNames,'???',model.compNames(I));
end

%The mets are matched to model.metNames and assigned to "compartment"
if eqnType==2
    %%Check that the metabolite names aren't present in the same compartment.
    %Not the neatest way maybe..
    t1=strcat(mets,'???',compartment);
    [I J]=ismember(t1,t2);
 
    if ~all(I)
        if allowNewMets==true
            %Add the new mets
            metsToAdd.metNames=mets(~I);
            
            %Generate the ids
            m=maxCurrent+1:maxCurrent+numel(metsToAdd.metNames);
            metsToAdd.mets=strcat({'m'},num2str(m(:)));
            
            metsToAdd.compartments=cell(numel(metsToAdd.mets),1);
            metsToAdd.compartments(:)={compartment};
            newModel=addMets(newModel,metsToAdd);
        else
            throw(MException('','One or more equations contain metabolites that are not in model.metNames. Set allowNewMets to true to allow this function to add metabolites or use addMets to add them before calling this function'));
        end
    end
    
    %Calculate the indexes of the metabolites
    metIndexes=J;
    metIndexes(~I)=numel(newModel.mets)-sum(~I)+1:numel(newModel.mets);
end

%The equations are on the form metNames[compName]
if eqnType==3
    %Parse the metabolite names
    metNames=cell(numel(mets),1);
    compartments=metNames;
    for i=1:numel(mets)
        starts=max(strfind(mets{i},'['));
        ends=max(strfind(mets{i},']'));
        
        %Check that the formatting is correct
        if isempty(starts) || isempty(ends) || ends<numel(mets{i})
        	throw(MException('',['The metabolite ' mets{i} ' is not correctly formatted for eqnType=3']));
        end
        
        %Check that the compartment is correct
        compartments{i}=mets{i}(starts+1:ends-1);
        I=ismember(compartments(i),newModel.compNames);
        if ~I
            throw(MException('',['The metabolite ' mets{i} ' has a compartment that is not in model.compNames']));
        end
        metNames{i}=mets{i}(1:starts-1);
    end
    
    %Check if the metabolite exists already
    t1=strcat(metNames,'???',compartments);
    [I J]=ismember(t1,t2);
 
    if ~all(I)
        if allowNewMets==true
            %Add the new mets
            metsToAdd.metNames=metNames(~I);
            
            %Generate the ids
            m=maxCurrent+1:maxCurrent+numel(metsToAdd.metNames);
            metsToAdd.mets=strcat({'m'},num2str(m(:)));
            
            metsToAdd.compartments=compartments(~I);
            newModel=addMets(newModel,metsToAdd);
            if newModel.first==1
             %   newModel.S=0;
            end
        else
            throw(MException('','One or more equations contain metabolites that are not in model.metNames. Set allowNewMets to true to allow this function to add metabolites or use addMets to add them before calling this function'));
        end
    end
    
    %Calculate the indexes of the metabolites
    metIndexes=J;
    metIndexes(~I)=numel(newModel.mets)-sum(~I)+1:numel(newModel.mets);
end

%Add the info to the stoichiometric matrix and to the rxnGeneMat. I do this
%in a loop, but maybe not necessary
if newModel.first==0
   newModel.S=[newModel.S sparse(size(newModel.S,1),nRxns)];
end
for i=1:nRxns
    newModel.S(metIndexes,nOldRxns+i)=S(:,i);
    
    %Parse the grRules and add to rxnGeneMat
    if isfield(newModel,'grRules')
       rule=newModel.grRules{nOldRxns+i};
       rule=strrep(rule,'(','');
       rule=strrep(rule,')','');
       rule=strrep(rule,' or ',' ');
       rule=strrep(rule,' and ',' ');
       genes=regexp(rule,':','split');
       [I J]=ismember(genes,newModel.genes);
       if ~all(I) && any(rule)
            throw(MException('',['Not all genes for reaction ' rxnsToAdd.rxns{i} ' were found in model.genes. If needed, add genes with addGenes before calling this function']));
       end
       if any(rule)
            newModel.rxnGeneMat(nOldRxns+i,J)=1;
       else
            %If there are no genes for the reaction, the rxnGeneMat should
            %still be resized
            % We do not need it in reconstruction
            %Ibrahim
            %newModel.rxnGeneMat=[newModel.rxnGeneMat;sparse(1,numel(newModel.genes))];
          
       end

    end
end
if newModel.first==1
    newModel.first=0;
end
end

function I=getInteger(s)
    %Checks if a string is on the form "m1" and if so returns the value of
    %the integer
    I=0;
    if strcmpi(s(1),'m')
        t=str2double(s(2:end));
        if ~isnan(t) && ~isempty(t)
            I=t;
        end
    end
end