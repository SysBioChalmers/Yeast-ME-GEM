function [AA weight]= aaComposition(model,solFile,mu)
sol=readSoplex_results(solFile);
flux=sol.X;
[a b Proteins]=xlsread('TableS1.xlsx','gene_seq');
[a b Proteins_mw]=xlsread('TableS1.xlsx','MW');
% get amino acids
weight=0;
AA=zeros(20,1);
for r=11272:numel(model.rxns)
     %fprintf('%d\n',r)
     name=regexp(cell2mat(model.rxns(r)),'_translation');
     if numel(name)>=1 
         %fprintf('%d\n',r);
         gene=regexp(cell2mat(model.rxns(r)),'\w*_translation','match');
         gene=strrep(gene,'_translation','');
          

         gene_deg={sprintf('%s_subunit_degradation',cell2mat(gene))};
         rxn_deg=find(ismember(model.rxns,gene_deg));
         rxn_syn=find(ismember(model.rxns,model.rxns(r)));
         MW_index=find(ismember(Proteins_mw(:,1),gene));
         
         index=find(ismember(Proteins(:,1),gene));
         seq=cell2mat(Proteins(index(1),2));
         
         con=sum(flux(rxn_syn))/(mu )*1000000;%(mu + sum(flux(rxn_deg)))*1000000; % convert mmol to fmol
         if numel(MW_index)>0
                      weight=weight+sum(flux(rxn_syn))/(mu )*cell2mat(Proteins_mw(MW_index,2))/1000;

         else
             fprintf('%s\n',cell2mat(gene));
         end
         AA=AA+aaConentration(con,seq);
         %get sequence.
         
     end
end
AA=AA/1000000;