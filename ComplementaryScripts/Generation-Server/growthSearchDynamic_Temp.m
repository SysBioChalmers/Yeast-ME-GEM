function [sol1 glc]= growthSearchDynamic_Temp(model, medium,atpCost,uptake,T, copies, overExpressed, complex_list,compartment,gr_min, gr_max,fileName,SoplexVersion,active_translation)
%This function searches for growth rate (mu) for uptake rate (uptake). 
%It is a binary search
%Ibrahim Elsemman


sol1.mu=0;
key='myModel';
glc=0;
sigma=1;
k=1;
while (abs(gr_max - gr_min)>=1e-3)
  gr=(gr_max + gr_min)/2;
    
  sol=modelGenerationGECKO_Temp(model,medium,gr,atpCost,'min_glucose',uptake,T,copies, overExpressed, fileName,complex_list,compartment,SoplexVersion,1,active_translation);
  
  if (strcmp(sol.status, 'OPTIMAL')&& strcmp(SoplexVersion,'Soplex-1.6')) || (strcmp(sol.status(1:7), 'optimal')&& strcmp(SoplexVersion,'Soplex-3.0'))
      %%check for status1
        gr_min=gr;
        sol1=sol;
        sol1.mu=gr_min;
        %system(sprintf('rm %s --force',fileName));
        fileName=strrep(fileName,'.lp',sprintf('_left_%d.lp',k));
      %  pause(30);
        k=k+1;
    else
       gr_max=gr;
       %system(sprintf('rm %s --force',fileName));
       fileName=strrep(fileName,'.lp',sprintf('_right_%d.lp',k));
        pause(30);
        k=k+1;
    end
    fprintf('growth rate <%f,%f>\n',gr_min,gr_max);
    
    
end
