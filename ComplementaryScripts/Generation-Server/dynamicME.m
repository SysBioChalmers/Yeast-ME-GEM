function [t y]=dynamicME()
%Simulate simple model for glucose and growth rate mu
%glucose uptake rate

%growth
time_interval=[0 8];
Variables = ['mu','glc'];
intialValues = [0.013,11.11];


[t,y] = ode45(@DFBAproblem,time_interval,intialValues);



function dydt = DFBAproblem(t,y)
%
global MyModel complex_list compartment;
atpCost.gam=33;
atpCost.ngam=1;
T=33;%T=T1(k);
overExpressed.syn='mCherry1_complex';
 overExpressed.dilultion='mCherry1_complex_dilution';
 overExpressed.degradation='mCherry1_degradation';
 copies=0;
vmax_glc=10; % ~13
km_glucose=0.8;

mu= y(1);
glc= y(2);

%convert concenteration to flux bound
v_glc= vmax_glc*(glc/(glc+km_glucose));
[sol Vglc]= growthSearchDynamic(MyModel, 'minimum',atpCost,v_glc,T, copies, overExpressed, complex_list,compartment,0, .5,sprintf('Dynamic_ME.lp'),'Soplex-1.6',1);
mu=sol.mu;
%fprintf('%f\t%f\t%f\n',t,sol.mu,Vglc);
dydt = [mu*y(1); -Vglc*y(1)];

