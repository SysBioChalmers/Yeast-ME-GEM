function model=complex_medium(model,index)
%index is medium type
% 3 YPD
% 5 SC
%4 SD
%6 Rich_AA

[a b metabolite]=xlsread('TableS1.xlsx','Medium');
for i=13:numel(metabolite(:,1))
    eq=sprintf(' => %s',cell2mat(metabolite(i,2)));
    rxnID=sprintf('%s_import',cell2mat(metabolite(i,1)));
    rxnNames=rxnID;
    model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,cell2mat(metabolite(i,index)),0,{''});
end