function step_average = randomScanning(a,p,d)

pos=[];
step_average=0;

for t=1:10000
    pos(t)=0;
    step=0;
  
while ~(pos (t) == a)
    r1=rand();
    r2=rand();
    if r2 < d
        pos(t)=0;
    end
    if r1 < p
        pos(t)=pos(t) + 1;%floor(step_length*rand());
        
    else
        pos(t)=pos(t) - 1;%floor(step_length*rand());
    end
    step=step+1;
end
step_average=step_average+step;
end
step_average= floor(step_average/10000);
