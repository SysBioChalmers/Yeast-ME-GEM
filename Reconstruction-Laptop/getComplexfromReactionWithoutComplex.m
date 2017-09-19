function [reaction_list complex_list compartment substrate_list product_list  lb ub type_list] = getComplexfromReactionWithoutComplex(model,saveFile)
reaction_list={};
complex_list={};
compartment={};
k=1;
for i=1:numel(model.grRules)
    i
rule=model.grRules{i};
[comp S P lb ub]=getCompartment(model,i);

rule1=strrep(rule,' AND ',':');
rule1=strrep(rule1,' OR ',';');
rule1=strrep(rule1,'(','');
rule1=strrep(rule1,')','');
complex=regexp(rule1,';','split');


for j=1:numel(complex)
    c=cell2mat(strrep(complex(j),':','_'));
    Sr=S;
    Pr=P;
    if strcmp(c,'')==0
    if numel(complex)>=1
        S = sprintf('%s + %s_complex_forward[%s]',S,c,comp);
        P = sprintf('%s + %s_complex_forward_out[%s]',P,c,comp);
    end
    end
    reaction_list(k,1)={sprintf('%s_forward',cell2mat(model.rxns(i)))};
    complex_list(k,1)=complex(j);
    compartment(k,1)={comp};
    substrate_list(k,1)={S};
    product_list(k,1)={P};
    ub_list(k,1)=ub;
    lb_list(k,1)=lb;
    type_list(k,1)={'forward'};    
    k=k+1;
    
    if lb<0
        if strcmp(complex,'')==0
        if numel(complex)>=1
            Sr = sprintf('%s + %s_complex_reverse_out[%s]',Sr,c,comp);
            Pr = sprintf('%s + %s_complex_reverse[%s]',Pr,c,comp);
        end
        end
    reaction_list(k,1)={sprintf('%s_reverse',cell2mat(model.rxns(i)))};
    complex_list(k,1)=complex(j);
    compartment(k,1)={comp};
    substrate_list(k,1)={Pr};
    product_list(k,1)={Sr};
    ub_list(k,1)=ub;
    lb_list(k,1)=0;
    type_list(k,1)={'reverse'};
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
