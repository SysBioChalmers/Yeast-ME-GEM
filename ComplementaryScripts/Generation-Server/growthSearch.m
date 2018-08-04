function [sol1 glc]= growthSearch(model, uptake,T,complex_list,compartment,gr_min, gr_max)
%This function searches for growth rate (mu) for uptake rate (uptake). 
%It is a binary search
%Ibrahim Elsemman



key='myModel';
glc=0;
sigma=1;
while (abs(gr_max - gr_min)>=1e-4)
    gr=(gr_max + gr_min)/2;
    
    system('rm myModel.lp --force');
    
    sol=modelGenerationNew(model,gr,uptake,T, sprintf('MaxMu_T%d.lp',T),complex_list,compartment);
    
    if strcmp(sol.status, 'OPTIMAL')==1
        gr_min=gr;
        glc=sol.X(find(ismember(model.rxns,'r_1714_reverse')));
        sol1=sol;
        sol1.mu=gr_min;
    else
        gr_max=gr;
    end
    fprintf('growth rate <%f,%f>\n',gr_min,gr_max);
    
    
end
