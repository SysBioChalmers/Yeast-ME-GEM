function misFoldingYFPHsp10444(model,fptr,tempFactor,mu,kcat,kdeg,sigmaTemp,RefoldingRatio)
[a b Proteins_length]=xlsread('TableS1.xlsx','Length');

      
kcat = 1*sigmaTemp;
kdeg = kdeg;
c11=kcat / (mu + kdeg);
c22=mu / (mu + kdeg);
c33=kdeg / (mu + kdeg);

kcat = 20*sigmaTemp;
kdeg = kdeg;
c111=kcat / (mu + kdeg);
c222=mu / (mu + kdeg);
c333=kdeg / (mu + kdeg);

kcat = 10*sigmaTemp;
kdeg = kdeg;
c1111=kcat / (mu + kdeg);
c2222=mu / (mu + kdeg);
c3333=kdeg / (mu + kdeg);



X=zeros(numel(model.rxns)+1,1);
%This function print the misfolding constraints
k=1;
ve='';
ve_mitochondrion='';
ve_mitochondrion_refolding='';

for r=1:numel(model.rxns)
   
    name=regexp(cell2mat(model.rxns(r)),'YFP_misfolding_cytoplasm');
    if numel(name)>0
        %extract gene name
        rxnName=cell2mat(model.rxns(r));
        dashIndex=regexp(rxnName,'_');
        geneName=rxnName(1:dashIndex(1)-1);
        %find 
        lengthIndex =find(ismember(Proteins_length(:,1),geneName));
        L=cell2mat(Proteins_length(lengthIndex,2));
        
        kf=exp(16.15 - 1.28*sqrt(L(1)))*3600; %per hour
        kf=kf*sigmaTemp; %temp effect
        
        c1 = mu /   (kf +kdeg+mu);
        c2 = kdeg / (kf +kdeg + mu);
        c3 = kf /  (kf + mu +kdeg) *RefoldingRatio;
        
        rxn_misfold=strrep(cell2mat(model.rxns(r)),'_folding','_misfolding');
        syn=find(ismember(model.rxns,rxn_misfold));
        
        rxn_misfold=strrep(cell2mat(model.rxns(r)),'_folding','_degradation_misfolding');
        deg=find(ismember(model.rxns,rxn_misfold));
        
        rxn_misfold=strrep(cell2mat(model.rxns(r)),'_folding','_dilution_misfolding');
        dil=find(ismember(model.rxns,rxn_misfold));
        
        rxn_misfold=strrep(cell2mat(model.rxns(r)),'_folding','_refolding');
        refold=find(ismember(model.rxns,rxn_misfold));
        
        if numel(syn)==0
            warning('no misfolding for this reaction %s\n',cell2mat(model.rxns(r)));
        else
            if X(syn(1))==0
                  fprintf(fptr,'%s',sprintf('MISFOLD%d: X%d - %.15f X%d =0\n',k, syn(1),tempFactor,r));
                 % fprintf(fptr,'%s',sprintf('MISFOLDdeg%d: X%d - %.15f X%d + X%d >=0\n',k,deg(1),c2,syn(1),refold(1)));
                  fprintf(fptr,'%s',sprintf('MISFOLDdil%d: X%d - %.15f X%d =0\n',k,dil(1),c1,syn(1)));
                  fprintf(fptr,'%s',sprintf('MISFOLDrefold%d: X%d - %.15f X%d =0\n',k,refold(1),c3,syn(1)));

                if mod(k,300)==0
                    sep=char(10);
                else
                    sep='';
                end
                
                name_compartment=regexp(cell2mat(model.rxns(r)),'_folding_mitochondrion');
                if numel(name_compartment)>0
                    if strcmp(ve_mitochondrion,'')
                        ve_mitochondrion =sprintf('X%d',syn(1));
                    else
                        ve_mitochondrion =sprintf('%s + X%d%c',ve_mitochondrion, refold,sep);
                    end
                    %refolding HSP60_HSP10
                    if strcmp(ve_mitochondrion_refolding,'')
                        ve_mitochondrion_refolding =sprintf('X%d',syn(1));
                    else
                        ve_mitochondrion_refolding =sprintf('%s + X%d%c',ve_mitochondrion_refolding, refold,sep);
                    end
                    
                else
                    if strcmp(ve,'')
                        ve =sprintf('X%d',syn(1));
                    else
                        ve =sprintf('%s + X%d%c',ve, syn(1),sep);
                    end
                end
                
                X(syn(1))=1;
                k=k+1;
            end
        end
        
    end
end

%coupling HSP104
rxnID=sprintf('Hsp104_HSP70_HSP40_biosynthesis');
syn=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Hsp104_HSP70_HSP40_degradation');
deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Hsp104_HSP70_HSP40_dilution');
dil=find(ismember(model.rxns,rxnID));

 fprintf(fptr,'PQC%d: %s - %.15f X%d <= 0\n',1,ve,c11,syn );
 fprintf(fptr,'PQC%d: X%d - %.15f X%d = 0\n',2,deg,c33,syn );
 fprintf(fptr,'PQC%d: X%d - %.15f X%d = 0\n',3,dil,c22,syn );

% %coupling mtHSP78
% rxnID=sprintf('mtHSP78_mtHSP70_biosynthesis');
% syn=find(ismember(model.rxns,rxnID));
% 
% rxnID=sprintf('mtHSP78_mtHSP70_degradation');
% deg=find(ismember(model.rxns,rxnID));
% 
% rxnID=sprintf('mtHSP78_mtHSP70_dilution');
% dil=find(ismember(model.rxns,rxnID));
% 
%  fprintf(fptr,'PQC%d: %s - %.15f X%d <= 0\n',1,ve_mitochondrion,c111,syn );
%  fprintf(fptr,'PQC%d: X%d - %.15f X%d = 0\n',2,deg,c333,syn );
%  fprintf(fptr,'PQC%d: X%d - %.15f X%d = 0\n',3,dil,c222,syn );
% 
%  %Coupling HSP60_HSP10
%  rxnID=sprintf('mtHSP10_mtHSP60_biosynthesis');
% syn=find(ismember(model.rxns,rxnID));
% 
% rxnID=sprintf('mtHSP10_mtHSP60_degradation');
% deg=find(ismember(model.rxns,rxnID));
% 
% rxnID=sprintf('mtHSP10_mtHSP60_dilution');
% dil=find(ismember(model.rxns,rxnID));
% 
% fprintf(fptr,'PQCHSP60%d: %s - %.15f X%d <= 0\n',1,ve_mitochondrion_refolding,c1111,syn );
% fprintf(fptr,'PQCHSP60%d: X%d - %.15f X%d = 0\n',2,deg,c3333,syn );
% fprintf(fptr,'PQCHSP60%d: X%d - %.15f X%d = 0\n',3,dil,c2222,syn );
%  