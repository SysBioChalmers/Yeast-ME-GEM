clear all
hold on
d=[0.002,0.001 0.0008 0.0004 0];
for k=1:numel(d)
t=0:20:1200;
for i=1:numel(t)
    step(i)=randomScanning(t(i),0.9,d(k));
end
plot(t,step,'o');
end
