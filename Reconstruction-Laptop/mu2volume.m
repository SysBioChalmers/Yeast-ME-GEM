function [v t]=mu2volume(mu)

for i=1:numel(mu)
t(i)=log(2)/mu(i)*60;
a=23.4366;
b=18.4648;
c=0.0342305;
d=2.47842;
v(i)=a + b* exp(-exp(-c*t(i) + d));
end
