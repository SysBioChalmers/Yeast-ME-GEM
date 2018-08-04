function [sol1 glc PredictedProteinMass newModel]= growthSearchDynamic(model, medium,atpCost,uptake,T, copies, overExpressed, complex_list,compartment,gr_min, gr_max,fileName,SoplexVersion,active_translation,ExpressionFileName)
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
    
  [sol newModel complexFileName molecule_number]=modelGenerationGECKO(model,medium,gr,atpCost,'min_carbon',uptake,T,copies, overExpressed, fileName,complex_list,compartment,SoplexVersion,1,active_translation);
  
  if (strcmp(sol.status, 'OPTIMAL')&& strcmp(SoplexVersion,'Soplex-1.6')) || (strcmp(sol.status(1:7), 'optimal')&& strcmp(SoplexVersion,'Soplex-3.0'))
      %%check for status1
        gr_min=gr;
        sol1=sol;
        sol1.mu=gr_min;
        sol1.fileName=fileName;
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
[geneList PredictedExpression ribosomeRatio] = modelProteinAllocation(model,sprintf('%s.out', sol1.fileName),sol1.mu,0.045,'Soplex-1.6');
 PredictedProteinMass = Estimate_abundance_metabolic(model,sol1.mu, 'complexName.txt',sprintf('%s.out', sol1.fileName), 'Soplex-1.6');
 if strcmp(ExpressionFileName,'')==0
 Sfactor= Estimate_abundance(model,sol1.mu, 'complexName.txt',sprintf('%s.out', sol1.fileName), ExpressionFileName,molecule_number,'Soplex-1.6',1);
 end

