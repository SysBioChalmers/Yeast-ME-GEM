function eq=makeReactionMito(tRNA, tRNA_modification, Substrate, Products,compartment)

if (strcmp(Substrate(1),'non'))
    eq=sprintf('%s%s => %s%s',tRNA,compartment,tRNA_modification, compartment);
else
    s=sprintf('%s%s + %s%s',tRNA,compartment,strtrim(cell2mat(Substrate(1))),compartment);
       
    for j=2:numel(Substrate)
        s=sprintf('%s + %s%s' ,s,strtrim(cell2mat(Substrate(j))),compartment);
    end
    

    p=sprintf('%s%s + %s%s',tRNA_modification,compartment,  strtrim(cell2mat(Products(1))),compartment);
    
    for j=2:numel(Products)
        p=sprintf('%s + %s%s' ,p,strtrim(cell2mat(Products(j))),compartment);
    end
    
    eq=sprintf('%s => %s',s,p);
    
end

eq=strrep(eq,'[n]',' [nucleus]');
eq=strrep(eq,'[c]',' [cytoplasm]');
eq=strrep(eq,'+ H [','+ H+ [');
eq=strrep(eq,'2 H [','2 H+ [');

eq=strrep(eq,'SAM [','S-adenosyl-L-methionine [');
eq=strrep(eq,'SAH [','S-adenosyl-L-homocysteine [');
eq= strrep(eq,'dimethylallyl diphosphate (C00235)','prenyl diphosphate');
eq= strrep(eq,'NH4','ammonium');

eq=strrep(eq,' + non [mitochondrion]','');

       
    