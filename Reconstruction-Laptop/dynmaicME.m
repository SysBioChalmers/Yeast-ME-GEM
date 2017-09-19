function dynamicME()
%Simulate simple model for glucose and growth rate mu
%glucose uptake rate

%growth 
time_interval=[0 5];
Variables = ['X','glc'];
intialValues = [0.05,11.1012];

%options=odeset('MaxStep',10);

[t,y] = ode15s(@DFBAproblem,time_interval,intialValues);
subplot(1,2,1)
plot(t,y(:,1),'-o') 
subplot(1,2,2)
plot(t,y(:,2),'-o') 



function dydt = DFBAproblem(t,y)
%
global newModel1;
vmax_glc=10;
km_glucose=2;

mu= y(1);
glc= y(2);

%convert concenteration to flux bound
v_glc= vmax_glc*(glc/(glc+km_glucose));
%[mu Vglc]= growthSearch(MyModel, v_glc,complex_list,compartment);
newModel1.ub(find(ismember(newModel1.rxns,'r_1714_reverse')))=v_glc;
newModel1.lb(find(ismember(newModel1.rxns,'r_1714_reverse')))=0;
sol=solveLP(newModel1);
mu=-sol.f;
Vglc=sol.x(find(ismember(newModel1.rxns,'r_1714_reverse')));
fprintf('%f\t%f\t%f\n',t,mu,Vglc);
dydt = [mu*y(1); -Vglc*y(1)];

