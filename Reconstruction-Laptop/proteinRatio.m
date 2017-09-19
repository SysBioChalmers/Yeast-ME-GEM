function r=proteinRatio(gr)
if gr <=0.28
    r=51.88*gr + 36.407;
else
    r=-48.527*gr+ 64.311;
end
r=r/100;
