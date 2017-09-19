close all;
clear ratio;
t=33:50;
topt=33;
tm=50;
for i=1:length(t)
    ratio(i)=(t(i)-topt)/(2*(tm-topt));
end
plot(t,ratio)