function sol=readSoplex_results(fileName,model)
fptr=fopen(fileName,'r');

line = fgets(fptr);
k=1;
obj='NO objeective';
X=zeros(numel(model.rxns)+5,1);%number of variables in LP %+1
while ischar(line)
    if k==133
        y=0;
    end
    var=regexp(line,'X\d*\t\d*\t\d*\.\d*','match');
    if numel(var)>=1
        index=str2num(cell2mat(strrep(regexp(line,'X\d*','match'),'X','')));
        value=str2num(cell2mat(regexp(line,'\d*\.\d*[eE][-+]\d*','match')));
        if index<=length(X)
               X(index,1)=value;
        else
           % msg=sprintf('%s has error in solution reading',fileName);
           % error(msg);
           
        end
    end
    s=regexp(line,'status=\w*([_][\w*])?','match');
    if numel(s)>=1
        status=strrep(cell2mat(s),'status=','');;
    end
    
    objective=regexp(line,'IEXAMP29 solution value is: [+-]?\d*\.\d*[eE][-+]\d*','match');
    if numel(objective)>=1
        obj=str2num(cell2mat(strrep(objective,'IEXAMP29 solution value is: ','')));
    end
    line = fgets(fptr);
    k=k+1;
end
if strcmp(obj,'NO objeective')==1
    sol.f=obj;
    sol.X=0;
    sol.status='No Solution';

else
    
sol.f=obj;
sol.X=X;
sol.status=status;
end
fclose(fptr);
