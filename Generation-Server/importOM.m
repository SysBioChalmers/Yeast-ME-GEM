function ve_tom = importOM(model)

k=0;
ve_tom='';

for i=1:numel(model.rxns)
    name=regexp(cell2mat(model.rxns(i)),'importing_mitochondrion_IMS');
    if numel(name)>0
        if strcmp(ve_tom,'')
            ve_tom = sprintf('X%d',i);
        else
            ve_tom = sprintf('%s + X%d%c',ve_tom,i,sep);
        end
    end
    
    k=k+1;
    if mod(k,300)==0
        sep=char(10);
    else
        sep='';
    end
    
end
