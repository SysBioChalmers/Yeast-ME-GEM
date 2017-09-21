function [rxnID rxnNames eq modification]=tRNAreactionsMito(base,modification_list,Reaction,Position,antiCodon)
eq={};
index=modification_list;
for k=1:numel(modification_list);
    i=index(k);
    modification_index=find(ismember(Reaction(1:58,4),Position(i,4)));
    modification_position=cell2mat(Reaction(modification_index,3));
    J=find(modification_position==cell2mat(Position(i,3)));
    %reaction index
    r=modification_index(J);
    
    Substrate=regexp(cell2mat(Reaction(r,5)),' \+ ','split');
    Products=regexp(cell2mat(Reaction(r,6)),' \+ ','split');
    
    compartment=' [mitochondrion]';
    
    if (strcmp(base,''))
        base=strtrim(sprintf('tRNA_%s_%s_%s',cell2mat(Position(i,1)),cell2mat(Position(i,2))),antiCodon);
        
    end
    
    modification=strtrim(sprintf('%s_%s_%d',base,cell2mat(Position(i,4)),cell2mat(Position(i,3))));
    
    rxnID {k,1} =sprintf('tRNA_%s_%s_%s_%s_%d_mitochondrion',cell2mat(Position(i,1)),cell2mat(Position(i,2)),antiCodon,cell2mat(Position(i,4)),cell2mat(Position(i,3)));
    rxnNames {k,1} = sprintf('tRNA_%s_%s_%s_%s_%d_mitochondrion',cell2mat(Position(i,1)),cell2mat(Position(i,2)),antiCodon,cell2mat(Position(i,4)),cell2mat(Position(i,3)));
    eq{k,1}=makeReactionMito(base,modification,Substrate, Products,compartment);
    base=modification;
    
end