seq='M';
for i=1:20
    for j=1:A(i)
       seq=sprintf('%s%s',seq,AAs{i});
    end
end
seq
aa2nt(seq)