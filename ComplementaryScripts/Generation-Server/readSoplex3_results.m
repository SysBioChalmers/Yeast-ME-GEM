function sol=readSoplex3_results(fileName,model)
fptr=fopen(fileName,'r');

line = fgets(fptr);
k=1;
status='NO objeective';
X=zeros(numel(model.rxns)+1,1);%number of variables in LP
while ischar(line)
  
    var=regexp(line,'X\d*\t\d*.\d*','match');
    if numel(var)>=1
        index=str2num(cell2mat(strrep(regexp(line,'X\d*','match'),'X','')));
        value=str2num(cell2mat(regexp(line,'\d*\.\d*','match')));
        if index<=length(X)
               X(index,1)=value;
        else
           % msg=sprintf('%s has error in solution reading',fileName);
           % error(msg);
           
        end
    end
    
    %the infeasible solution
    s=regexp(line,'SoPlex status       : problem is solved \[infeasible\]','match');
    if numel(s)>=1
       status='infeasible';
    end
    
    s=regexp(line,'SoPlex status       : problem is solved \[optimal\]');
    if numel(s)>=1
       status='optimal';
    end
    
    s=regexp(line,'SoPlex status       : error \[unspecified\]','match');
    if numel(s)>=1
       status='error';
    end
    
    
    
    objective=regexp(line,'Objective value     : [+-]?\d*\.\d*[eE][-+]\d*','match');
    if numel(objective)>=1
        obj=str2num(cell2mat(strrep(objective,'Objective value     : ','')));
    end
    line = fgets(fptr);
    k=k+1;
end
if strcmp(status,'optimal')==0
    sol.f=obj;
    sol.X=0;
    
    if strcmp(status,'error')==1
      sol.status='optimal';
      sol.status1='Unspecific';
    else
        sol.status='infeasible';
    end

else
    
sol.f=obj;
sol.X=X;
sol.status=status;
end
fclose(fptr);
