function model=addYeastReaction(model,eq,rxnID,rxnName,lb,ub,c,rule,gene)
          
          enz={''};
          r.equations={eq};
          r.rxns={rxnID};
          r.rxnNames={rxnName};
          r.lb=lb;
          r.ub= ub;
          r.c=c;
          r.grRules=rule;
          r.subSystems={''};
          r.eccodes=gene;
          model=addNewRxns(model,r,3,'1',true);
