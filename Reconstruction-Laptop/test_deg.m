degradation_constrain='';
degradation_constrain1='';
k=1;
for r=1:numel(model.rxns)
    name=regexp(cell2mat(model.rxns(r)),'_subunit_degradation');
    if numel(name)>0
        if k <= 800
        if strcmp(degradation_constrain,'')
            degradation_constrain = sprintf('X%d',r);
        else
            degradation_constrain = sprintf('%s + X%d',degradation_constrain,r);
        end
        else
            if strcmp(degradation_constrain1,'')
            degradation_constrain1 = sprintf('X%d',r);
        else
            degradation_constrain1 = sprintf('%s + X%d',degradation_constrain1,r);
        end
            
        end
        k=k+1;
    end
end    
constrain=sprintf('%s - XDEG1 =0\n XDEG1 + %s - =0\n',degradation_constrain,degradation_constrain1);

