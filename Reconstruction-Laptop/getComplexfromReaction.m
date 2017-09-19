function [newModel, reaction_list complex_list compartment substrate_list product_list  lb ub type_list] = getComplexfromReaction(model,newModel,saveFile)
reaction_list={};
complex_list={};
compartment={};
k=1;
for i=1:numel(model.grRules)

rule=model.grRules{i};
[comp S P lb ub C]=getCompartment(model,i);

rule1=strrep(rule,' AND ',':');
rule1=strrep(rule1,' OR ',';');
rule1=strrep(rule1,'(','');
rule1=strrep(rule1,')','');
complex=regexp(rule1,';','split');


for j=1:numel(complex)
    c=cell2mat(strrep(complex(j),':','_'));
    Sr=regexprep(S,'>','_GT_');
    Pr=regexprep(P,'>','_GT_');
    eq=sprintf('%s => %s',Sr,Pr);
%     if strcmp(c,'')==0
% %     if numel(complex)>=1
% %         S = sprintf('%s + %s_complex_forward[%s]',S,c,comp);
% %         P = sprintf('%s + %s_complex_forward_out[%s]',P,c,comp);
% %     end
%     end
    reaction_list(k,1)={sprintf('%s_forward',cell2mat(model.rxns(i)))};
    complex_list(k,1)=complex(j);
    compartment(k,1)={comp};
    substrate_list(k,1)={S};
    product_list(k,1)={P};
    ub_list(k,1)=ub;
    lb_list(k,1)=lb;
    type_list(k,1)={'forward'};   
    newModel=addYeastReaction(newModel,eq,reaction_list(k,1),{''},0,ub,C,complex_list(k,1));
    k=k+1;
    
    if lb<0
        if strcmp(complex,'')==0
%         if numel(complex)>=1
%             Sr = sprintf('%s + %s_complex_reverse_out[%s]',Sr,c,comp);
%             Pr = sprintf('%s + %s_complex_reverse[%s]',Pr,c,comp);
%         end
        end
    reaction_list(k,1)={sprintf('%s_reverse',cell2mat(model.rxns(i)))};
    complex_list(k,1)=complex(j);
    compartment(k,1)={comp};
    substrate_list(k,1)={Pr};
    product_list(k,1)={Sr};
    ub_list(k,1)=ub;
    if lb==-1000
      lb_list(k,1)=0;
    else
      lb_list(k,1)=lb;
    end
    type_list(k,1)={'reverse'};

    eq=sprintf('%s => %s',Pr,Sr);
    newModel=addYeastReaction(newModel,eq,reaction_list(k,1),{''},0,ub,C,complex_list(k,1));

    k=k+1;

    
    end
end
end

if(saveFile)
fptr=fopen('output_reaction.txt','w');
n=k-1;
for k=1:n
    fprintf(fptr,'%s\t%s\t%s\t%s\t%s\t%f\t%f\t%s\n',cell2mat(reaction_list(k,1)), cell2mat(substrate_list(k,1)), cell2mat(product_list(k,1)),cell2mat(complex_list(k,1)),cell2mat(type_list(k,1)), lb_list(k),ub_list(k),cell2mat(compartment(k,1)));
end
fclose(fptr)
end
