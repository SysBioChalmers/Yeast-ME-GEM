k=1;
for i=12925:14209
    sol=testReaction(MyModel,cell2mat(MyModel.rxns(i)));
    if sol.f==0
        fprintf('%d %s\n',i,cell2mat(MyModel.rxns(i)));
        k=k+1;
    end
end