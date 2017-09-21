function sol=modelGeneration(MyModel,media,mu1,atpCost,objective,glucoseUptake,T,copies, overExpressed, fileName,complex_list,compartment)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');
[a b Proteins_length]=xlsread('TableS1.xlsx','Length');
[a b subunit]=xlsread('TableS1.xlsx','Subunit');
[a b Protein_kcat]=xlsread('TableS1.xlsx','kcat');
[a b Protein_info]=xlsread('TableS1.xlsx','Annotation');
[a b Protein_plasma]=xlsread('TableS1.xlsx','Plasma_main');


k=0;
key=strrep(fileName,'.lp','');
command =sprintf('/zhome/88/9/107870/soplex-1.6.0/bin/soplex  -s0 -x  -q -C %s > %s.out %s',fileName,fileName);

sigma_plasma = 1;
    
c_mim=100;%31 is the best

% Tm=42;% 
% Tf=0;
% To=33;
% 
% TH=(2*Tm - T -To)/(2*(Tm-To));
% TL=(2*Tf - T -To)/(2*(Tf-To));

% Tm=43;% 
% Tf=7;
% Topt=33;
Tm=45;% 
Tf=10;
Topt=32;


if T>=Topt
    UnFoldRatio=(T-Topt)/(2*(Tm-Topt));
else
    UnFoldRatio=(T-Topt)/(2*(Tf-Topt));
end

Dt=1/(273.15 + T) - 1/(273.15 + Topt);

Df=exp(-1000*Dt);
sigmaTemp=Df;


mu=mu1; 
molecule_number=800; 
mitochondrionVolume=1000; %Ratio of mitochondrion to cell volume is equal 1-2%
cyto_volume =1000;
plasma_area=1000;
plasma_budjet=40;
Vmax_NH4 = 15;

sigma= 1;%3.059*mu + 0.1475;%the line is fitted from fermentation data
upr=0.1;

if sigma > 1
  sigma= 1;   
end
   

%the sector partition
Metabolic_Sector=0*proteinRatio(mu);
Ribosome_Sector=0*proteinRatio(mu);
Engergy_Sector=0*Metabolic_Sector;
Maintaince_Sector=upr*proteinRatio(mu);

kdeg=0.045;
kcat = 110*60*60;

% if T<=To
%     kcat=kcat*TL;
% else
%     kcat=kcat*TL;
% end
K=sigma*sigmaTemp*kcat/(mu + kdeg);
K1=sigma*kcat/(mu);
c1= mu / (kdeg + mu);
c2= kdeg / (kdeg + mu);
c1_unmodel=c1;
c2_unmodel=c2;

gam=atpCost.gam;
ngam=atpCost.ngam;

weight_constant = 1/(mu + kdeg);

weight_sigma =sigma*60*60/(mu + kdeg);

f=56517.2/1000.0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

m_rr=1.9e6;
m_aa=109;
f_rRNA=0.8;
kt=(1/0.1774);
r0=0.1445;

c_ribosome= m_rr/(m_aa * f_rRNA);

kcat_rate = c_ribosome * kt * (mu + kdeg)/(mu+ r0 *kt);

kcat_ribo= kcat_rate * sigmaTemp;

kdeg_ribo= kdeg;

c_ribo = kcat_ribo / (kdeg_ribo + mu);
c_ribo_deg = kdeg_ribo / (kdeg_ribo + mu);
c_ribo_dil = mu / (kdeg_ribo + mu);
c_ribo_constant =1347851.15/(1000* (kdeg_ribo + mu));

%translation factors
kcat_ribo_scanning = 10*60*60*sigmaTemp;
kcat_ribo_eEF1A =6*60*60*sigmaTemp;
kcat_ribo_eEF1B =6*60*60*sigmaTemp;
kcat_ribo_eEF2 =10*60*60*sigmaTemp;
kcat_ribo_eEF3 =940*60*sigmaTemp;
kcat_ribo_eRF=0.17*60*60*sigmaTemp;


kcat_ribo_mito=15*60*60;
kdeg_ribo_mito= kdeg;
c_ribo_mito = kcat_ribo_mito / (kdeg_ribo_mito + mu);
c_ribo_deg_mito = kdeg_ribo_mito / (kdeg_ribo_mito + mu);
c_ribo_dil_mito = mu / (kdeg_ribo_mito + mu);

%Assembly Factors
kcat_assembly=2000*60*sigmaTemp;
kdeg_asembly= kdeg;
c_assembly = kcat_assembly / (kdeg_asembly + mu);
c_assembly_deg = kdeg_asembly / (kdeg_asembly + mu);
c_assembly_dil = mu / (kdeg_asembly + mu);

%Proteasome
kcat_Proteasome=5*60*sigmaTemp;
kdeg_Proteasome= kdeg;
c_Proteasome = kcat_Proteasome / (kdeg_Proteasome + mu);
c_Proteasome_deg = kdeg_Proteasome / (kdeg_Proteasome + mu);
c_Proteasome_dil = mu / (kdeg_Proteasome + mu);

kcat_AARS = 3.5*60*60*sigmaTemp;
sigma_AARS=1;
K_AARS=sigma_AARS*kcat_AARS/(mu + kdeg);


%fileName=sprintf('Problem_test.lp');
model=MyModel;
model= fixReaction(model);

model = Inactivate_Glucose(model);

model.c=zeros(numel(model.c),1);
%eq=sprintf('lipid [cytoplasm] + 0.963 (1-_GT_3)-beta-D-glucan [cytoplasm] + 0.994 mannan [cytoplasm] + 0.667 glycogen [cytoplasm] + 0.085 trehalose [cytoplasm] + 0.051 AMP [cytoplasm] + 0.00359 dAMP [cytoplasm] + 0.051 GMP [cytoplasm] + 0.05 CMP [cytoplasm] + 0.067 UMP [cytoplasm] + 0.00243 dCMP [cytoplasm] + 0.00243 dGMP [cytoplasm] + 0.00359 dTMP [cytoplasm] => biomass_test [cytoplasm]');
%eq=sprintf('lipid [cytoplasm] + 1.14 (1-_GT_3)-beta-D-glucan [cytoplasm] + 0.821 mannan [cytoplasm] + 0.519 glycogen [cytoplasm] + 0.0234 trehalose [cytoplasm] + 0.051 AMP [cytoplasm] + 0.00359 dAMP [cytoplasm] + 0.051 GMP [cytoplasm] + 0.05 CMP [cytoplasm] + 0.067 UMP [cytoplasm] + 0.00243 dCMP [cytoplasm] + 0.00243 dGMP [cytoplasm] + 0.00359 dTMP [cytoplasm] => biomass_test [cytoplasm]');
eq=biomass_CEN_PK(mu);
rxnID='biomass_without_protein';
rxnNames=sprintf('biomass test');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,1000,0,{''});

eq=sprintf('%s => ','biomass_test [cytoplasm]');
rxnID='test_biomass';
rxnNames=sprintf('test biomass');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},mu,mu,0,{''});
%set biomass
model.c=zeros(numel(model.c),1);
model.c(find(ismember(model.rxns,'Ribosome_biosynthesis')))=0;

%gene deletion
% model=deleteGene(model,cell2mat(Strains(mutant)));
% model=deleteGene(model,'YDL085W');
% model=deleteGene(model,'YMR145C');
% model=deleteGene(model,'YOL059W');
% model.ub(find(ismember(model.rxns,'r_2115_forward')))=0;
% model.lb(find(ismember(model.rxns,'r_2115_forward')))=0;


% Fix the growth Media
model=Glucose_Growth(model,glucoseUptake,0);

if strcmp(media,'YPD')==1
    model=complex_medium(model,3);%YPD media
elseif strcmp(media,'SC')==1
    model=complex_medium(model,5);%SC media
elseif strcmp(media,'SD')==1
    model=complex_medium(model,4);%SD media
else
    %the growth media is minimum
end
%control NH4 uptake rate
model = closeAmmoniumOnly(model,Vmax_NH4,0);
%model = closeO2(model);
%model = Ethanol_Growth(model,1000,-1);
%model = Glycerol_Growth(model,1000,-1);


eq=sprintf('%f H2O [cytoplasm] + %f ATP [cytoplasm] => %f phosphate [cytoplasm] + %f ADP [cytoplasm] + %f H+ [cytoplasm]', gam,gam,gam,gam,gam);
rxnID='GAM';
rxnNames=sprintf('GAM reaction');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},mu,mu,0,{''});

eq=sprintf('H2O [cytoplasm] + ATP [cytoplasm] => phosphate [cytoplasm] + ADP [cytoplasm] + H+ [cytoplasm]');
rxnID='NGAM';
rxnNames=sprintf('NGAM reaction');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},ngam,ngam,0,{''});


%write to LP file 

fptr=fopen(fileName,'w');
fprintf(fptr,'Maximize\n');
if strcmp(objective,'min_glucose')==1
    index_obj=find(ismember(model.rxns,'r_1714_reverse'));
    fprintf(fptr,'obj: - X%d\n',index_obj);
elseif strcmp(objective,'max_unmodelled_protein')==1
    index_obj=find(ismember(model.rxns,'Unmodeled_Protein_biosynthesis'));
    fprintf(fptr,'obj: X%d\n',index_obj);
else
    fprintf(fptr,'obj: 0.000\n');
end



fprintf(fptr,'Subject To\n');
for i=1:numel(model.mets)
    J=find(full(model.S(i,:)));
    
    for k=1:numel(J)
        s=full(model.S(i,J(k)));
        if mod(k,200)==0
            sep=char(10);
        else
            sep='';
        end
        if k==1
           eq=sprintf('%.15f X%d',s,J(k));
        else
           if s>0
               eq=sprintf('%s + %.15f X%d%c',eq,s,J(k),sep);
           else
               eq=sprintf('%s %.15f X%d%c',eq,s,J(k),sep);
           end
        end
    end
    
    fprintf(fptr,'C%d: %s = 0\n',i,eq);

end


k1=0;

mitoProteins={'Q0045','Q0080','Q0085','Q0105','Q0130','Q0250','Q0275'};
report_reaction=zeros(numel(complex_list),1);

content='';
Engery_content='';
for i=1:numel(complex_list)
    
    if (1)% strcmp(compartment(i),'mitochondrion')==0
        peptide=regexp(complex_list{i},':','split');
       
        if numel(find(ismember(mitoProteins,peptide)))==0
            
            if strcmp(complex_list{i},'')==0 && report_reaction(i)==0
                %compute the molecular weight for the complex
                MW=0;
                
                    %complex is dimer 
                    if numel(peptide)==1
                        ns=1;
                        I=find(ismember(subunit(:,1),peptide));
                        if numel(I)==1
                            nsubunit=strrep(cell2mat(subunit(I,2)),'A','');
                            nsubunit=sprintf('%s ',nsubunit);
                            ns=str2num(nsubunit);
                        end
                    elseif (numel(peptide)>1)  %complex is a hetrodimer.
                        I=find(ismember(subunit(:,1),complex_list{i}));
                        ns=ones(numel(peptide),1);
                        if numel(I)==1
                            subunits=regexp(cell2mat(subunit(I,2)),',','split');
                            for jj=1:numel(subunits)
                               ns(jj)=str2num(subunits{jj});
                             end
                        end
                    end   
                for m1=1:numel(peptide)                              
                    peptide_index=find(ismember(Proteins_MW(:,1),strrep(peptide(m1),'-','')));
                    MW = MW + ns(m1) * cell2mat(Proteins_MW(peptide_index,2));
                end
                
                k1=k1+1;
                c=cell2mat(strrep(complex_list(i),':','_'));
                rxnID=sprintf('%s_complex_formation_%s',c,strrep(compartment{i},' ',''));
                rxnID=strrep(rxnID,'-','');
                
                s=find(ismember(model.rxns,rxnID));
                %add the condition collecting collect all catalyzed reactions by this complex
                rxns_by_this_complex=find(ismember(model.grRules,complex_list{i}));
                
                rxns=sprintf('X%d',rxns_by_this_complex(1));
                report_reaction(rxns_by_this_complex(1)) =1;
                for r=2:numel(rxns_by_this_complex)
                   rxns=sprintf('%s + X%d',rxns, rxns_by_this_complex(r)); 
                   report_reaction(rxns_by_this_complex(r)) =1;
                end
                %print kcat condition
                   index_kcat = find(ismember(Protein_plasma(:,1),complex_list{i}));
                   if numel(index_kcat)==1
                       if cell2mat(Protein_kcat(index_kcat,2)) <= 10
                           KK=weight_sigma*cell2mat(Protein_kcat(index_kcat,2));
                       else
                           if strcmp(cell2mat(Protein_kcat(index_kcat,3)),'G')==1
                               KK=sigmaTemp*3600*cell2mat(Protein_plasma(index_kcat,4))*cell2mat(Protein_plasma(index_kcat,3))*sigma_plasma;
                           else
                               KK=sigmaTemp*3600*cell2mat(Protein_plasma(index_kcat,4))*cell2mat(Protein_plasma(index_kcat,3))*sigma_plasma;
                           end
                       end
                       fprintf(fptr,'CM%d: %s - %.15f X%d <= 0\n',k1,rxns,KK,s); %K for average kcat or KK for different kcats
                   else
                       fprintf(fptr,'CM%d: %s - %.15f X%d <= 0\n',k1,rxns,K,s);
                   end
                
                if mod(i,300)==0
                    sep=char(10);
                else
                    sep='';
                end
        
                
                if strcmp(content,'')==1
                   content = sprintf('%.15f X%d',weight_constant*MW/1000.0,s);
                else
                   content = sprintf('%s + %.15f X%d%c',content, weight_constant*MW/1000.0,s,sep);
                end
                
              
                rxnID=sprintf('%s_complex_degradation_%s',c,strrep(compartment{i},' ',''));
                rxnID=strrep(rxnID,'-','');
                
                d=find(ismember(model.rxns,rxnID));
                fprintf(fptr,'CM%d: X%d - %.15f X%d = 0\n',k1,d(1),c2,s );
                
                rxnID=sprintf('%s_complex_dilution_%s',c,strrep(compartment{i},' ',''));
                rxnID=strrep(rxnID,'-','');
                
              %    k1=k1+1;
                d=find(ismember(model.rxns,rxnID));
                fprintf(fptr,'CM%d: X%d - %.15f X%d = 0\n',k1,d(1),c1,s);
                
            end
        else % mito 
            
            I=find(ismember(subunit(:,1),complex_list{i}));
            ns=ones(numel(peptide),1);
            if numel(I)==1
                subunits=regexp(cell2mat(subunit(I,2)),',','split');
                for jj=1:numel(subunits)
                    ns(jj)=str2num(subunits{jj});
                end
            end
            
            MW=0;
            for m1=1:numel(peptide)
                peptide_index=find(ismember(Proteins_MW(:,1),strrep(peptide(m1),'-','')));
                MW = MW + ns(m1) * cell2mat(Proteins_MW(peptide_index,2));
            end
                
                
                k1=k1+1;
                c=cell2mat(strrep(complex_list(i),':','_'));
                rxnID=sprintf('%s_complex_formation_mitochondrion',c);
                rxnID=strrep(rxnID,'-','');
                
                s=find(ismember(model.rxns,rxnID));
                
                index_kcat = find(ismember(Protein_kcat(:,1),complex_list{i}));
                if numel(index_kcat)==1
                    KK=weight_sigma*cell2mat(Protein_kcat(index_kcat,2))/sigma;
                    fprintf(fptr,'CM%d: X%d - %.15f X%d <= 0\n',k1,i,K,s); % K for average kcat
                else
                    fprintf(fptr,'CM%d: X%d - %.15f X%d <= 0\n',k1,i,K,s);
                    
                end
                
               %degradation
               rxnID=sprintf('%s_complex_degradation_mitochondrion',c);
               rxnID=strrep(rxnID,'-','');
               d=find(ismember(model.rxns,rxnID));
               fprintf(fptr,'CMEMDEG%d: X%d - %.15f X%d <= 0\n',k1,d,c2,s);
               
               %dilution
               rxnID=sprintf('%s_complex_dilution_mitochondrion',c);
               rxnID=strrep(rxnID,'-','');
               dil=find(ismember(model.rxns,rxnID));
               fprintf(fptr,'CMMEMDIL%d: X%d - %.15f X%d <= 0\n',k1,dil,c1,s);
               
               
               content = sprintf('%s + %.15f X%d%c',content, (1/mu)*MW/1000.0,s,sep);
               
               if strcmp(Engery_content,'')==1
                   Engery_content = sprintf('%.15f X%d%c', (1/mu)*MW/1000.0,s,sep);         
               else
                   Engery_content = sprintf('%s + %.15f X%d%c',Engery_content, (1/mu)*MW/1000.0,s,sep);
               end
                
        end
        
    end
end
%fprintf(fptr,'MetabolicContent: %s >= %f \n',content,0.30*proteinRatio(mu));
%Glycolysis Subsector
%fprintf(fptr,'MetabolicContent: %s <= %f \n',GlycolysisCondition(model, 0.1,0.04),Glycolysis_Subsector); 
%fprintf(fptr,'MetabolicEngery: %s <= %f \n',Engery_content,Engergy_Sector);


% adding loadig amino acids constrains

k1=0;
for i=1:numel(model.rxns)
    name=regexp(cell2mat(model.rxns(i)),'_loading','match');
    if numel(name)>0 && report_reaction(i)==0
        if strcmp(model.grRules(i),'')==0
            peptide=regexp(cell2mat(model.grRules(i)),':','split');
            fprintf('%s\n',cell2mat(model.rxns(i)));
            if strcmp(peptide,'')==0
                
                k1=k1+1;
                c=cell2mat(strrep(model.grRules(i),':','_'));
                
                %check if AARS is cytoplasm or Mito
                name=regexp(cell2mat(model.rxns(i)),'_loading_mitochondrion','match');
                if numel(name)>0
                    rxnID=sprintf('%s_complex_formation_mitochondrion',c);
                    comp='mitochondrion';
                    
                else
                rxnID=sprintf('%s_complex_formation_cytoplasm',c);
                comp='cytoplasm';
                end
                rxnID=strrep(rxnID,'-','');
                
                s=find(ismember(model.rxns,rxnID));
                
               %add the condition collecting collect all catalyzed reactions by this complex
                rxns_by_this_complex=find(ismember(model.grRules,model.grRules{i}));
                
                rxns=sprintf('X%d',rxns_by_this_complex(1));
                report_reaction(rxns_by_this_complex(1)) =1;
                for r=2:numel(rxns_by_this_complex)
                   rxns=sprintf('%s + X%d',rxns, rxns_by_this_complex(r)); 
                   report_reaction(rxns_by_this_complex(r)) =1;
                end

                
                fprintf(fptr,'CAARS%d: %s - %.15f X%d <= 0\n',k1,rxns,K_AARS,s);
                
                % add the constrains about the  charged tRNA is less than or
                % equal the uncharged tRNA
                %         fprintf(fptr,'CAARS%d: X%d - X%d = 0\n',k1,i,i-1);
                
                
                rxnID=sprintf('%s_complex_degradation_%s',c,comp);
                rxnID=strrep(rxnID,'-','');
                
                k1=1+1;
                d=find(ismember(model.rxns,rxnID));
                fprintf(fptr,'CAARS%d: X%d - %.15f X%d = 0\n',k1,d(1),c2,s );
                
                rxnID=sprintf('%s_complex_dilution_%s',c,comp);
                rxnID=strrep(rxnID,'-','');
                
                k1=k1+1;
                d=find(ismember(model.rxns,rxnID));
                fprintf(fptr,'CAARS%d: X%d - %.15f X%d = 0\n',k1,d(1),c1,s);
                
            end
        end
    end
end

% add ribosome capacity constrain
rxnID=sprintf('Ribosome_biosynthesis');
syn_ribo=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Ribosome_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Ribosome_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

%add set of proteins biosyenthesis reactions


ve='';
ve_scanning='';
ve_eEF1A='';
ve_eEF1B='';
ve_eEF2='';
ve_eEF3='';

ve_total='';
ve_total_molecules='';
k=0;
r=1;
condition1='';
condition2='';
condition3='';

for i=1:numel(model.rxns)
    if i==6897
      y=0;
    end
    name=regexp(cell2mat(model.rxns(i)),'_translation');
    name1=regexp(cell2mat(model.rxns(i)),'_mitochondrion');
    if numel(name)>0 && numel(name1)==0
         
         if numel(name1)>0
             error('mito protein in cytosoloic ribosome');
         end
        
        k=k+1;
        %protein Length
        
        gene=strrep(model.rxns(i),'_translation','');
        gene_index=find(ismember(Proteins_length(:,1),gene));
        L= cell2mat(Proteins_length(gene_index(1),2));
        L_scanning= cell2mat(Proteins_length(gene_index(1),3));
        %some genes have no 5'UTR. We assume that the number of steps equal
        %1
        if L_scanning ==0
            L_scanning=1;
        end
        if L==0
            error(sprintf('Length of gene%d: %s is zero',r,cell2mat(gene)));
        end
        %print the translation constraint
        fprintf(fptr,'CTranslation%d: X%d - %.15f R%d =0\n',r,i,kcat_ribo/L,r);
        fprintf(fptr,'CTranslationScanning%d: X%d - %.15f IF%d =0\n',r,i,kcat_ribo_scanning/L_scanning,r);
        
        % elongation
        fprintf(fptr,'CTranslationElongation1A%d: X%d - %.15f eEF1A%d =0\n',r,i,kcat_ribo_eEF1A/L,r);
        fprintf(fptr,'CTranslationElongation1B%d: X%d - %.15f eEF1B%d =0\n',r,i,kcat_ribo_eEF1B/L,r);
        fprintf(fptr,'CTranslationElongation2%d: X%d - %.15f eEF2%d =0\n',r,i,kcat_ribo_eEF2/L,r);
        fprintf(fptr,'CTranslationElongation3%d: X%d - %.15f eEF3%d =0\n',r,i,kcat_ribo_eEF3/L,r);
        %termination
        fprintf(fptr,'CTranslationTermination%d: X%d - %.15f eRF%d =0\n',r,i,kcat_ribo_eRF,r);
         
        %misfolding effect
        
        gene=strrep(model.rxns(i),'_translation','');
        gene_index=find(ismember(Proteins_MW(:,1),gene));
        MW=cell2mat(Proteins_MW(gene_index,2))/1000/(mu + kdeg);
        
        %mass fraction
        gene=strrep(model.rxns(i),'_translation','');
        gene_index=find(ismember(Proteins_MW(:,1),gene));
        MW=cell2mat(Proteins_MW(gene_index,2))/1000/(mu + kdeg);
       
        if mod(k,300)==0
            sep=char(10);
        else
            sep='';
        end
        
        if strcmp(ve,'')
                ve = sprintf('R%d',r);
                ve_scanning = sprintf('IF%d',r);
                ve_eEF1A = sprintf('eEF1A%d',r);
                ve_eEF1B = sprintf('eEF1B%d',r);
                ve_eEF2 = sprintf('eEF2%d',r);
                ve_eEF3 = sprintf('eEF3%d',r);
                ve_eRF = sprintf('eRF%d',r);
                %total mass
                ve_total=sprintf('%f X%d',MW, i);
                ve_total_molecules=sprintf('X%d', i);
                
        else
               %ribosome
               ve=sprintf('%s + R%d%c',ve, r,sep);
               %scanning
               ve_scanning=sprintf('%s + IF%d%c',ve_scanning, r,sep);
               %elongation
               ve_eEF1A=sprintf('%s + eEF1A%d%c',ve_eEF1A, r,sep);
               ve_eEF1B=sprintf('%s + eEF1B%d%c',ve_eEF1B, r,sep);
               ve_eEF2=sprintf('%s  + eEF2%d%c',ve_eEF2, r,sep);
               ve_eEF3=sprintf('%s  + eEF3%d%c',ve_eEF3, r,sep);
               ve_eRF=sprintf('%s  + eRF%d%c',ve_eRF, r,sep);
               %total mass
               ve_total=sprintf('%s + %f X%d%c',ve_total,MW,i,sep);
               ve_total_molecules=sprintf('%s + X%d%s', ve_total_molecules,i,sep);
        end
        r=r+1;
    end
end
 

fprintf(fptr,'CR%d: %s - %f X%d <= 0\n',1,ve,1/(mu+kdeg), syn_ribo );
fprintf(fptr,'CR%d: X%d - %.15f X%d = 0\n',2,syn_deg,c_ribo_deg,syn_ribo );
fprintf(fptr,'CR%d: X%d - %.15f X%d = 0\n',3,syn_dil,c_ribo_dil,syn_ribo );

%Translation factors
%scanning
[syn_IF dil_IF deg_IF]=getETF(MyModel,{'Intialtion_Factors'});
fprintf(fptr,'CRS%d: %s - %f X%d <= 0\n',1,ve_scanning,1/(mu+kdeg), syn_IF );
fprintf(fptr,'CRS%d: X%d - %.15f X%d = 0\n',2,deg_IF,c_ribo_deg,syn_IF );
fprintf(fptr,'CRS%d: X%d - %.15f X%d = 0\n',3,dil_IF,c_ribo_dil,syn_IF );

%eEF1A
[syn_IF dil_IF deg_IF]=getETF(MyModel,{'YPR080W'});
fprintf(fptr,'CREF1A%d: %s - %f X%d <= 0\n',1,ve_eEF1A,1/(mu+kdeg), syn_IF );
fprintf(fptr,'CREF1A%d: X%d - %.15f X%d = 0\n',2,deg_IF,c_ribo_deg,syn_IF );
fprintf(fptr,'CREF1A%d: X%d - %.15f X%d = 0\n',3,dil_IF,c_ribo_dil,syn_IF );
% % %EF1B
[syn_IF dil_IF deg_IF]=getETF(MyModel,{'eEF1B'});
fprintf(fptr,'CREF1B%d: %s - %f X%d <= 0\n',1,ve_eEF1B,1/(mu+kdeg), syn_IF );
fprintf(fptr,'CREF1B%d: X%d - %.15f X%d = 0\n',2,deg_IF,c_ribo_deg,syn_IF );
fprintf(fptr,'CREF1B%d: X%d - %.15f X%d = 0\n',3,dil_IF,c_ribo_dil,syn_IF );

  % %eEF2
  [syn_IF dil_IF deg_IF]=getETF(MyModel,{'YDR385W'});
  [syn_IF1 dil_IF1 deg_IF1]=getETF(MyModel,{'YOR133W'});
%  
  fprintf(fptr,'CReEF2%d: %s - %f X%d - %f X%d <= 0\n',1,ve_eEF2,1/(mu+kdeg), syn_IF, 1/(mu+kdeg), syn_IF1 );
%  
  fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',2,deg_IF,c_ribo_deg,syn_IF );
  fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',3,dil_IF,c_ribo_dil,syn_IF );
%  
  fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',4,deg_IF1,c_ribo_deg,syn_IF1 );
  fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',5,dil_IF1,c_ribo_dil,syn_IF1);

 %eEF3
 [syn_IF dil_IF deg_IF]=getETF(MyModel,{'YLR249W'});
 [syn_IF1 dil_IF1 deg_IF1]=getETF(MyModel,{'YNL014W'});
 
 fprintf(fptr,'CReEF2%d: %s - %f X%d - %f X%d <= 0\n',1,ve_eEF3,1/(mu+kdeg), syn_IF, 1/(mu+kdeg), syn_IF1 );
 
 fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',2,deg_IF,c_ribo_deg,syn_IF );
 fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',3,dil_IF,c_ribo_dil,syn_IF );
 
 fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',4,deg_IF1,c_ribo_deg,syn_IF1 );
 fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',5,dil_IF1,c_ribo_dil,syn_IF1);
 
 %eRF
 [syn_IF dil_IF deg_IF]=getETF(MyModel,{'Relase_Factors'});
 fprintf(fptr,'CReEF2%d: %s - %f X%d  <= 0\n',1,ve_eRF,1/(mu+kdeg), syn_IF );
 fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',2,deg_IF,c_ribo_deg,syn_IF );
 fprintf(fptr,'CReEF2%d: X%d - %.15f X%d = 0\n',3,dil_IF,c_ribo_dil,syn_IF );
 
 
% add mitocondiral ribosome capacity constrain
rxnID=sprintf('MitoRibosome_biosynthesis');
syn_ribo=find(ismember(model.rxns,rxnID));

rxnID=sprintf('MitoRibosome_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('MitoRibosome_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

%add set of proteins biosyenthesis reactions
ve_mito='';
for i=1:numel(model.rxns)
    name1=regexp(cell2mat(model.rxns(i)),'_translation_mitochondrion');
    if numel(name1)>0 
        
    if strcmp(ve_mito,'')
        ve_mito=sprintf('X%d',i);
    else
        ve_mito=sprintf('%s + X%d%c',ve_mito,i,sep);
    end
      gene=strrep(model.rxns(i),'_translation_mitochondrion','');
      
      gene_index=find(ismember(Proteins_MW(:,1),gene));
      MW=cell2mat(Proteins_MW(gene_index,2))/1000/(mu + kdeg);
      
      %total mass
      ve_total=sprintf('%s + %f X%d',ve_total,MW,i);
      ve_total_molecules=sprintf('%s + X%d', ve_total_molecules,i);
    end
end

fprintf(fptr,'CMR%d: %s - %.15f X%d <= 0\n',1,ve_mito,c_ribo_mito,syn_ribo );
fprintf(fptr,'CMR%d: X%d - %.15f X%d = 0\n',2,syn_deg,c_ribo_deg_mito,syn_ribo );
fprintf(fptr,'CMR%d: X%d - %.15f X%d = 0\n',3,syn_dil,c_ribo_dil_mito,syn_ribo );


%protein content constrains
fprintf(fptr,'CtotalProtein: %s = %f\n',ve_total,proteinRatio(mu)/(1-UnFoldRatio));
total_number=molecule_number*1.0e6/(13*6.023e8)*(mu+kdeg);
fprintf(fptr,'CveTotal: %s <= %.15f\n',ve_total_molecules, total_number);


%fprintf(fptr,'RibosomeSector: %.15f X%d >= %f\n',c_ribo_constant,syn_ribo,Ribosome_Sector);
%mRNA constrain
%fprintf(fptr,'CmRNA%d: %s <= 9e-04 \n',1,ve );

 
% %relase factors
%  genes={'Relase_Factors'};
%  kdeg=kdeg; %Petri et al
%  kcat= 0.17*60*60;% per hour from  
%  c1=kcat / (mu + kdeg);
%  c2=mu / (mu + kdeg);
%  c3=kdeg / (mu + kdeg);
%  printeETF(model,fptr,genes, c1,c2,c3,ve);
% 


%Chaperone constrain
%Chaperone constrain
%misFoldingConstraint(MyModel,fptr,UnFoldRatio);
misFoldingConstraintHsp1044(model,fptr,UnFoldRatio,mu,5,kdeg,sigmaTemp);

% SSB
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'SSB');

% RAC
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'RAC');

% SSA
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'SSA');

% SSE
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'SSE');

% HSP40
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'HSP40');

% HSP90
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'HSP90');

% % PFD
 kcat = 60*sigmaTemp;
 kdeg = kdeg;
 c1=kcat / (mu + kdeg);
 c2=mu / (mu + kdeg);
 c3=kdeg / (mu + kdeg);
% 
 addChaperoneConstrain(fptr,model,c1,c2,c3,'PFD');

% CCT
kcat = 60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'CCT');



%add content constrains
%fprintf(fptr,'CCMP1: %s + = 0.25\n',content );
%mass_content=sprintf('%s + %.15f X%d \n',content,1347851.15/(1000.0*(mu+kdeg_ribo)),syn_ribo);

%add the boundary

%add intermembrane
fprintf(fptr,'CMIM: %s <= %f\n',membraneCondition(model, mu,kdeg),c_mim);
fprintf(fptr,'PMR:   %s <= %f\n',plasmaCondition(model, mu,kdeg),plasma_area);
%fprintf(fptr,'PMR1   %s <= %f\n', outMembraneConstraints(),plasma_budjet);
%NH4 Vmax
fprintf(fptr,'PMR2: X4172 <= %f\n',Vmax_NH4);
fprintf(fptr,'PMR3: X4631 <= 15\n');
fprintf(fptr,'PMR2: X4681 <= 15\n');

fprintf(fptr,'MitochondrionVolume: %s <= %f\n',mitoVolumeCondition(model, mu,kdeg),mitochondrionVolume);


%add Cytosolic Constrain
%sc=cytosolCondition(model, mu,kdeg);


%soplex needs 8109 characters in the line. We devide the condition into
%two lines.
   %[constrain1 constrain2 constrain3 constrain4]=cytosolConditionAll(model, mu,kdeg);
   [constrain1 constrain2 constrain3 ]=cytosolCondition(model, mu,kdeg);
   fprintf(fptr,'CCV1: %s  - XXCV =  0  \n',constrain1);
   fprintf(fptr,'CCV2: %s  + XXCV - XXCV2  =  0  \n',constrain2);
   fprintf(fptr,'CCV4: %s  + XXCV2 <= %f\n',constrain3,cyto_volume);

    %[constrain1 constrain2 constrain3 constrain4]=cytosolConditionAll(model, mu,kdeg);
    %printf(fptr,'CCV1: %s  - XXCV =  0  \n',constrain1);
    %fprintf(fptr,'CCV2: %s  + XXCV - XXCV2  =  0  \n',constrain2);
    %fprintf(fptr,'CCV3: %s  + XXCV2 - XXCV3  =  0  \n',constrain3);
    %printf(fptr,'CCV4: %s  + XXCV3 <= %f\n',constrain4,cyto_volume);


%add the assembly factors constrains
rxnID=sprintf('Assembly_Factors_biosynthesis');
syn_ribo=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Assembly_Factors_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Assembly_Factors_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Ribosome_biosynthesis');
syn_ribo_cyto=find(ismember(model.rxns,rxnID));


fprintf(fptr,'CAF%d: X%d - %.15f X%d <= 0\n',1,syn_ribo_cyto,c_assembly,syn_ribo );
fprintf(fptr,'CAF%d: X%d - %.15f X%d = 0\n',2,syn_deg,c_assembly_deg,syn_ribo );
fprintf(fptr,'CAF%d: X%d - %.15f X%d = 0\n',3,syn_dil,c_assembly_dil,syn_ribo );

% %add the proteasome constrain
rxnID=sprintf('proteasome_biosynthesis');
syn_ribo=find(ismember(model.rxns,rxnID));

rxnID=sprintf('proteasome_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('proteasome_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

%get all degradation reaction
rxnID=sprintf('Ribosome_biosynthesis');
syn_ribo_cyto=find(ismember(model.rxns,rxnID));
% 

degradation_constrain='';
degradation_constrain1='';
k=1;
for r=1:numel(model.rxns)
    name=regexp(cell2mat(model.rxns(r)),'_subunit_degradation');
    if numel(name)>0
        if k <= 800
            if strcmp(degradation_constrain,'')
                degradation_constrain = sprintf('X%d',r);
            else
                degradation_constrain = sprintf('%s + X%d',degradation_constrain,r);
            end
        else
            if strcmp(degradation_constrain1,'')
                degradation_constrain1 = sprintf('X%d',r);
            else
                degradation_constrain1 = sprintf('%s + X%d',degradation_constrain1,r);
            end
            
        end
        k=k+1;
    end
end

fprintf(fptr,'CPD1: %s - X%d =0\nCPD2: X%d + %s - %.15f X%d =0\n',degradation_constrain,numel(model.rxns)+1,numel(model.rxns)+1,degradation_constrain1,c_Proteasome, syn_ribo);
fprintf(fptr,'CPD3: X%d - %.15f X%d = 0\n',syn_deg,c_Proteasome_deg,syn_ribo );
fprintf(fptr,'CPD4: X%d - %.15f X%d = 0\n',syn_dil,c_Proteasome_dil,syn_ribo );

%add unmodeled proteins
rxnID=sprintf('Unmodeled_Protein_biosynthesis');
syn=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Unmodeled_Protein_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('Unmodeled_Protein_dilution');
syn_dil=find(ismember(model.rxns,rxnID));


fprintf(fptr,'CUMP%d: X%d - %.15f X%d = 0\n',1,syn_deg,c2_unmodel,syn );
fprintf(fptr,'CUMP%d: X%d - %.15f X%d = 0\n',2,syn_dil,c1_unmodel,syn );
f=48260*weight_constant/1000.0 ;
fprintf(fptr,'MaintainceSector: %.15f X%d  >= %f\n',f,syn,Maintaince_Sector);

%protein import in mitochondrion
kcat = 60*60;
kdeg=kdeg;

c1=kcat/(mu + kdeg);
c2= mu/(mu+kdeg);
c3=kdeg/(mu+kdeg);

addTIM23Constraint(fptr,model,c1,c2,c3);
addTIM22Constraint(fptr,model,c1,c2,c3);
addTOMConstraint(fptr,model,c1,c2,c3);

%add glycine importer into mitochondrion

rxnID=sprintf('YDL119C_complex_formation');
syn=find(ismember(model.rxns,rxnID));

rxnID=sprintf('YDL119C_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('YDL119C_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

fprintf(fptr,'CNP%d: X%d - %.15f X%d = 0\n',1,syn_deg,c1,syn );
fprintf(fptr,'CNP%d: X%d - %.15f X%d = 0\n',2,syn_dil,c2,syn );


fprintf(fptr,'CNP%d: X4371 - %.15f X%d <= 0\n',3,K,syn );

% %overexpression constraint, such as YFP or GFP proteins
%     e=copies/(13*6.023e8);
%     MMu=mu*e;
%     KKdeg= kdeg*e;
%     Syn=(mu+kdeg)*e;
%     
%      fprintf(fptr,'Overexpreesed1: X%d = %.15f\n',find(ismember(model.rxns, overExpressed.syn)),Syn);
%      fprintf(fptr,'Overexpreesed2: X%d = %.15f\n',find(ismember(model.rxns, overExpressed.dilultion)),MMu);
%      fprintf(fptr,'Overexpreesed3: X%d = %.15f\n',find(ismember(model.rxns, overExpressed.degradation)),KKdeg);

%     model.ub(find(ismember(model.rxns,overExpressed.dilultion)))=MMu;
%     model.lb(find(ismember(model.rxns,overExpressed.dilultion)))=MMu;
%     model.ub(find(ismember(model.rxns,overExpressed.degradation)))=KKdeg;
%     model.lb(find(ismember(model.rxns,overExpressed.degradation)))=KKdeg;
% 

%
%Print the bounds for each flux.
fprintf(fptr,'Bounds\n');
one_copy_per_cell=(mu+kdeg)*(1/(13*6.023e8));
N_cpc=20000000;
for i=1:numel(model.rxns)
    name=regexp(cell2mat(model.rxns(i)),'_complex_formation','match');
    if numel(name)>0
         fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,N_cpc*one_copy_per_cell);
    else
        
        if model.ub(i)>=100
            fprintf(fptr,'%f <= X%d <= +infinity\n',model.lb(i),i);
        else
            fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,model.ub(i));
        end
    end
    
end
fprintf(fptr,'End\n');
fclose(fptr);
%



fprintf('finished\n');

system(command,'-echo');
system(sprintf('./myResults1 %s',strrep(fileName,'.lp','')),'-echo');

fileNameOutput=sprintf('%s.out',fileName);
sol=readSoplex_results(fileNameOutput,model);