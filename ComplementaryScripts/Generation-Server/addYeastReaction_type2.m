function model=addYeastReaction_check_type2(model,eq,rxnID,rxnName,lb,ub,c,rule)
          enz={''};
          r.equations={eq};
          r.rxns=rxnID;
          r.rxnNames=rxnName;
          r.lb=lb;
          r.ub= ub;
          r.c=c;
          r.grRules=rule;
          r.subSystems={''};
          r.eccodes=enz;
          
          g=regexp(cell2mat(rule),':','split');
          for jj=1:numel(g)
              ng=find(ismember(model.genes,g(jj)));
              if numel(ng)==0
                  model.genes(numel(model.genes)+1)=g(jj);
              end
          end
    
         model=addNewRxns(model,r,2,'1',true);
          
              
          
