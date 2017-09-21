seq= 'ATG';
for i=1:numel(Codon)
    nc = regexp(cell2mat(Codon(i)),';','split');
    for j=1:numel(nc)
        nac = regexp(cell2mat(nc(j)),' ','split');
        
        for k=1:str2num(cell2mat(nac(1)))
            seq= sprintf('%s%s',seq,cell2mat(nac(2)));
        end
    end
end

