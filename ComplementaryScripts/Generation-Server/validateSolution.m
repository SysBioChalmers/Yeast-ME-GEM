function validateSolution(fileName,model)
%fileName='maltose_batch_NDI_left_1_left_2_right_3_left_4_left_5_right_6_left_7_right_8.lp';
fileNameResult=sprintf('%s.out', fileName);
fptr=fopen(fileName,'r');
%model=MyModel;
sol=readSoplex_results(fileNameResult,model);

%read LP
line = fgets(fptr);
line = fgets(fptr);
line = fgets(fptr);
line = fgets(fptr);
k=1;
ke=1;
msg={''};
 while ischar(line) 
    % repeat until find the complete line
%     fprintf('%d\n',k);
%     if (k==14242)
%         y=0;
%     end
    k=k+1;
    
    operator_index=regexp(line,'=');
    myCondition='';
    while numel(operator_index)==0
        if strcmp(myCondition,'')==1
            myCondition=line(1:length(line)-1);
        else
            myCondition=sprintf('%s %s', myCondition,line(1:length(line)-1));
        end
         operator_index=regexp(line,'=');
        if numel(operator_index)==0
            line = fgets(fptr);
        end

        

    end
    if strcmp(myCondition,'')==0
        line=myCondition;
    end
    operator='';
    %condition name
    condition_index=regexp(line,':');
    condition_name=line(1:condition_index-1);
    
 
    %condition operator
    operator_index=regexp(line,' =');
    if numel(operator_index)>0
        operator ='eq';
    else
        operator_index=regexp(line,'>=');
        if numel(operator_index)>0
        operator ='gt';
        else
        operator_index=regexp(line,'<=');
        if numel(operator_index)>0
        operator ='lt';
        end
        end
    end
    
    condition_eq=line(condition_index+1:operator_index-1);
    condition=condition_eq;
    %subsutite in condition
    var=regexp(condition,'X\d*','match');
    for i=1:numel(var)
        vi=str2num(strrep(cell2mat(var(i)),'X',''));
        if i==1 && numel(regexp(condition_eq(1:2),' X'))==1
           var_val=sprintf('%.15g',sol.X(vi));
        else
            var_postion=regexp(condition_eq,cell2mat(var(i)));
            if i<=numel(var) && numel(regexp(condition_eq(var_postion-2:var_postion),'+ X'))==1
                       var_val=sprintf('%.15g',sol.X(vi));
            elseif i<=numel(var) && numel(regexp(condition_eq(var_postion-2:var_postion),'- X'))==1
                       var_val=sprintf('%.15g',sol.X(vi));
             else
            var_val=sprintf('* %.15g',sol.X(vi));
            end
        end
        if i==numel(var)
            condition=strrep(condition,sprintf(' %s',cell2mat(var(i))),sprintf(' %s',var_val));
        else
            condition=strrep(condition,sprintf(' %s ',cell2mat(var(i))),sprintf(' %s ',var_val));
        end
    end
    %remove chr(10)
    %condition=strrep(condition,'\n','');
    %fprintf('%s\n',condition_name);
    
       if numel(find(ismember({'CtotalProtein'},condition_name)))==1
          eps_c=eval(condition);
       else
           try
             eps_c=eval(condition); 
           catch
               msg{ke}=sprintf('Error in computing %s',condition_name);
               ke=ke+1;
           end
       end
    %fprintf('%s;%s::%s::%.15f;%s\n',condition_name,condition_eq,condition,condition_value,operator); 
  
        
    if (abs(eps_c)>1.0e-10)
      fprintf('%s::%.15g\n',condition_name,eps_c);
    end
    
    line = fgets(fptr);
    if  numel(regexp(line,'Bounds'))>0
        break;
    end
end