function  rr=isLPvisible(filename)

fid=fopen(filename,'r');
tline = fgets(fid);
while ischar(tline)
    s1=strfind(tline,'Problem status  : PRIMAL_AND_DUAL_FEASIBLE');
    if numel(s1)>0
        tline = fgets(fid);
        s2=strfind(tline,'Solution status : OPTIMAL');
        if numel(s2)>0
            tline = fgets(fid);
            rr=tline;
        end
      
    
    end
    tline = fgets(fid);
end
fclose(fid);