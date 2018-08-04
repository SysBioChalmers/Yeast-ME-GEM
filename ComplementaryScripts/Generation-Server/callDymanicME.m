% call Dynmaic ME
clear all
fprintf('Calling Dynmaic ME\n');
global MyModel complex_list compartment;
load NewYeastnextGeneration.mat

[t y]=dynamicME()
save DynmaicME_results