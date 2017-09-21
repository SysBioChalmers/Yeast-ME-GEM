function sol=readSoplex_results(fileName)
fptr=fopen(fileName,'r');

line = fgets(fptr);
k=1;
while ischar(line)
    if k==133
        y=0;
    end
    var=regexp(line,'X\d*\t\d*\t\d*\.\d*','match');
    if numel(var)>=1
        index=str2num(cell2mat(strrep(regexp(line,'X\d*','match'),'X','')));
        value=str2num(cell2mat(regexp(line,'\d*\.\d*[eE][-+]\d*','match')));
        X(index,1)=value;
    end
    s=regexp(line,'status=\w*([_][\w*])?','match');
    if numel(s)>=1
        status=strrep(cell2mat(s),'status=','');;
    end
    
    objective=regexp(line,'Solution value is: [+-]?\d*\.\d*[eE][-+]\d*','match');
    if numel(objective)>=1
        obj=str2num(cell2mat(strrep(objective,'Solution value is: ','')));
    end
    line = fgets(fptr);
    k=k+1;
end

sol.f=obj;
sol.X=X;
sol.status=status;

fclose(fptr);
