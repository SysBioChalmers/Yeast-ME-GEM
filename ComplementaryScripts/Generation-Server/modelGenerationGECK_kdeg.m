function [sol model complexFileName molecule_number]=modelGenerationGECKO_kdeg(MyModel,media,mu1,atpCost,objective,glucoseUptake,T,copies, overExpressed, fileName,complex_list,compartment,SoplexVersion,Iteration,active_translation)

[a b Proteins_MW]=xlsread('TableS1.xlsx','MW');
[a b Proteins_length]=xlsread('TableS1.xlsx','Length');
[a b subunit]=xlsread('TableS1.xlsx','Subunit');
[a b Protein_kcat]=xlsread('TableS1.xlsx','kcat');
[a b Protein_info]=xlsread('TableS1.xlsx','Annotation');
[a b Protein_plasma]=xlsread('TableS1.xlsx','Plasma_main');
[a b KcatGEKO]=xlsread('TableS1.xlsx','KcatGEKO');


k=0;
key=strrep(fileName,'.lp','');
%
%SoplexPath='/zhome/88/9/107870/solpex-solvers/soplex-3.0.1'
%SoplexPath='/zhome/88/9/107870/solpex-solvers/soplex-3.0.1/bin/soplex -i2147483647 -f0 -o0 --writebas=mybase.bas -c  -x --int:checkmode=2 -i2147483647 -q --solvemode=2 --readmode=1 --int:syncmode=1';
SoplexPath='/zhome/88/9/107870/soplex-2.2.0/bin/soplex -f1e-20  -o1e-20 --int:syncmode=1 --solvemode=2 --readmode=1 -x -q -c'; 
%Soplex 1.6
if strcmp(SoplexVersion,'Soplex-1.6')==1
  command =sprintf('/zhome/88/9/107870/soplex-solvers/soplex-1.6.0/bin/soplex  -s0 -x  -q -C %s > %s.out %s',fileName,fileName);
end
%soplex 3.0.1
if strcmp(SoplexVersion,'Soplex-3.0')==1
  command =sprintf('%s %s > %s.out %s',SoplexPath,fileName,fileName);
end

% Tm=42;% 
% Tf=0;
% To=33;
% 
% TH=(2*Tm - T -To)/(2*(Tm-To));
% TL=(2*Tf - T -To)/(2*(Tf-To));

% Tm=43;% 
% Tf=7;
% Topt=33;
Tm=41;% 
Tf=10;
Topt=33;
RefoldingRatio=0.0;

if T>=Topt
    UnFoldRatio=(T-Topt)/(2*(Tm-Topt));
else
    UnFoldRatio=0.0;
end


Dt=1/(273.15 + T) - 1/(273.15 + Topt);

Df=exp(-3200*Dt);
sigmaTemp=1;


%% Parameter setting
mu=mu1; 
sigma_plasma = 1;
c_mim=67; %67 is the best parameter

average_amino_acids=128;
average_protein_length=300;

molecule_number=floor(proteinRatio(mu)*13/(average_protein_length*average_amino_acids)*6.023e23/1e12); 
molecule_number_mu_max=200;%floor(proteinRatio(0.46)*13/(average_protein_length*average_amino_acids)*6.023e23/1e12); 

mitochondrionVolume=0.88;%Ratio of mitochondrion to cell volume is equal 1-2%
cyto_volume =4;
plasma_area=60;
plasma_budjet=10;
Vmax_NH4 = 7;

sigma=1;%3.7*mu;% + 0.1475;%the lobjine is fitted from fermentation data
upr=-1.6667*mu + 0.76;
upr=0.26; %best parameter is 0.26

    
if sigma>1
    sigma=1;
end

% if mu>=0.30
%     upr=0.30;
% end
%the sector partition
Metabolic_Sector=0*proteinRatio(mu);
Ribosome_Sector=0*proteinRatio(mu);
Engergy_Sector=0*Metabolic_Sector;
Maintaince_Sector=upr*proteinRatio(mu);
min_kcat=1;
kdeg=0.1;
kcat = 90*60*60;
kcat_for_missing=85;
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

weight_sigma =sigma/(mu + kdeg);
f=56517.2/1000.0;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                

m_rr=1.9e6;
m_aa=average_amino_acids;
f_rRNA=0.85;
kt=(1/0.22);
r0=0.14;

ribosome_active_ratio=active_translation;
if mu<=0.05
    ribosome_active_ratio= 4*mu + 0.75;

end

c_ribosome= m_rr/(m_aa * f_rRNA);

kcat_rate = c_ribosome * kt * (mu + kdeg)/(mu+ r0 *kt);

%kcat_rate= 20.657*mu + 1.4949;
kcat_ribo= kcat_rate * sigmaTemp;

kdeg_ribo= kdeg;

c_ribo = kcat_ribo / (kdeg_ribo + mu);
c_ribo_deg = kdeg_ribo / (kdeg_ribo + mu);
c_ribo_dil = mu / (kdeg_ribo + mu);
c_ribo_constant =1347851.15/(1000* (kdeg_ribo + mu));

%translation factors
sTF=0.8;
kcat_ribo_scanning = sTF*kcat_ribo*60*60*sigmaTemp;
kcat_ribo_eEF1A =sTF*6*60*60*sigmaTemp;
kcat_ribo_eEF1B =sTF*6*60*60*sigmaTemp;
kcat_ribo_eEF2 =sTF*kcat_ribo*60*60*sigmaTemp;
kcat_ribo_eEF3 =sTF*970*kcat_ribo*60*sigmaTemp;
kcat_ribo_eRF=sTF*0.17*60*60*sigmaTemp;

Mito_ribosome_active_ratio=1;
kcat_ribo_mito= (1/10)*20*60*60;
kdeg_ribo_mito= kdeg;
c_ribo_mito = kcat_ribo_mito / (kdeg_ribo_mito + mu);
c_ribo_deg_mito = kdeg_ribo_mito / (kdeg_ribo_mito + mu);
c_ribo_dil_mito = mu / (kdeg_ribo_mito + mu);

%Assembly Factors
kcat_assembly=(1/4)*2000*60*sigmaTemp;
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

kcat_AARS = 3*60*60*sigmaTemp;
sigma_AARS=1;
K_AARS=sigma_AARS*kcat_AARS/(mu + kdeg);


%fileName=sprintf('Problem_test.lp');
model=MyModel;
model= fixReaction(model);

model = Inactivate_Glucose(model);

%% Biomass Reaction
eq=biomass_CEN_PK(mu);
%fprintf('%s/n',eq);
rxnID='biomass_without_protein';
rxnNames=sprintf('biomass test');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},0,1000,0,{''});

eq=sprintf('%s => ','biomass_test [cytoplasm]');
rxnID='test_biomass';
rxnNames=sprintf('test biomass');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},mu,mu,0,{''});
%set biomass
if strcmp(objective,'min_carbon')==0
    model.c=zeros(numel(model.c),1);
end
model.c(find(ismember(model.rxns,'Ribosome_biosynthesis')))=0;

%%
%gene deletion
% model=deleteGene(model,cell2mat(Strains(mutant)));
% model=deleteGene(model,'YDL085W');
% model=deleteGene(model,'YMR145C');
% model=deleteGene(model,'YOL059W');
% model.ub(find(ismember(model.rxns,'r_2115_forward')))=0;
% model.lb(find(ismember(model.rxns,'r_2115_forward')))=0;


%%  growth Media

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
%model = closeAmmoniumOnly(model,Vmax_NH4,0);
%model = closeO2(model);
%model = Ethanol_Growth(model,1000,-1);
%model = Glycerol_Growth(model,1000,-1);

%% GAM and NGAM reactions
eq=sprintf('%f H2O [cytoplasm] + %f ATP [cytoplasm] => %f phosphate [cytoplasm] + %f ADP [cytoplasm] + %f H+ [cytoplasm]', gam,gam,gam,gam,gam);
rxnID='GAM';
rxnNames=sprintf('GAM reaction');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},mu,mu,0,{''});

eq=sprintf('H2O [cytoplasm] + ATP [cytoplasm] => phosphate [cytoplasm] + ADP [cytoplasm] + H+ [cytoplasm]');
rxnID='NGAM';
rxnNames=sprintf('NGAM reaction');
model=addYeastReaction(model,eq,{rxnID},{rxnNames},ngam,ngam,0,{''});


%% writing into LP file 
fptr=fopen(fileName,'w');
fprintf(fptr,'Maximize\n');
if strcmp(objective,'min_glucose')==1
    index_obj=find(ismember(model.rxns,'r_1714_reverse'));
    fprintf(fptr,'obj: - X%d\n',index_obj);
elseif strcmp(objective,'max_unmodelled_protein')==1
    index_obj=find(ismember(model.rxns,'Unmodeled_Protein_biosynthesis'));
    fprintf(fptr,'obj: X%d\n',index_obj);
elseif strcmp(objective,'min_carbon')==1
    index_obj=find(model.c==-1);
    fprintf(fptr,'obj: -X%d\n',index_obj);
else
    fprintf(fptr,'obj: 0.000\n');
end


fprintf(fptr,'Subject To\n');
%% SV=0
%This part write the constrain SV=0
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

%% V <= kcat [e]
%This part write the constrain 
k1=0;

mitoProteins={'Q0045','Q0080','Q0085','Q0105','Q0130','Q0250','Q0275'};
report_reaction=zeros(numel(model.grRules),1);

content='';
Engery_content='';
% estimate abundace
% fileName_abundance=fileAbundance.fileName;
% sheet_abundance=fileAbundance.sheet;
% Expression=Estimate_abundance(complex_list,fileName_abundance,sheet_abundance);
Expression_syn=zeros(numel(complex_list),1);
fileExpression=sprintf('S_%.2g.xlsx',0.1) % replace 0.1 with mu
%  if(sprintf('%s.out.sfactor',fileName),'file')==2
%  else
     
 if (Iteration==1)
     S_factor=ones(6000,1);
 else
    [a b c express_data express_model S_factor]=textread(strrep(sprintf('%s.out.sfactor',fileName),sprintf('k_%d',Iteration),sprintf('k_%d',Iteration-1)),'%s%s%s%f%f%f');

 end
 % print S factor for each complex
 complexFileName = strrep(fileExpression,'.xlsx','.complexName.txt');
 fptr_sfactor=fopen( complexFileName,'w');
 
ks=1;

for i=1:numel(complex_list) %check again for 1:n
    
    if (ks<=80000)% strcmp(compartment(i),'mitochondrion')==0
               
        peptide=regexp(complex_list{i},':','split');
        
        if numel(find(ismember(mitoProteins,peptide)))==0
            
            if strcmp(complex_list{i},'')==0 && report_reaction(i)==0
                fprintf(fptr_sfactor,'S%d\t%s\t%s\n',ks,complex_list{i},strrep(compartment{i},' ','_'));
                if ks<=numel(S_factor(:))
                    SS= S_factor(ks);
                    if SS==0
                        SS=sigma;
                    end
                else
                    SS=sigma;
                    
                end
              
              
              % SS=1;
              %scape this 
                ks=ks+1;
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
                %print the molecular mass of a complex in a file
                k1=k1+1;
                c=cell2mat(strrep(complex_list(i),':','_'));
                rxnID=sprintf('%s_complex_formation_%s',c,strrep(compartment{i},' ',''));
                rxnID=strrep(rxnID,'-','');
                
              s=find(ismember(model.rxns,rxnID));
              %add the condition collecting collect all catalyzed reactions by this complex
              rxns_by_this_complex=find(ismember(model.grRules,complex_list{i}));
              %find only the metabolic index
              metabolic_rxn_index  = find(rxns_by_this_complex<=numel(complex_list));
              rxns_by_this_complex= rxns_by_this_complex(metabolic_rxn_index);
              
              
              rxn_kcat=cell2mat(KcatGEKO(rxns_by_this_complex(1)+1,3));
              if rxn_kcat==0
                  rxn_kcat = kcat_for_missing;
              elseif rxn_kcat<=min_kcat
                  rxn_kcat=min_kcat;
                  
              end
               newKcat=3600*rxn_kcat; %convert per second to per hour
               rxns=sprintf('%.15f X%d', 1/newKcat, rxns_by_this_complex(1)); 
               
               %rxns=sprintf('X%d',rxns_by_this_complex(1));
               report_reaction(rxns_by_this_complex(1)) =1;
               sep1='';
               for r=2:numel(rxns_by_this_complex)
                   %%kcat for this reaction
                   rxn_kcat=cell2mat(KcatGEKO(rxns_by_this_complex(r)+1,3));
                   if rxn_kcat==0
                       rxn_kcat = kcat_for_missing;
                   elseif rxn_kcat<=min_kcat
                       rxn_kcat=min_kcat;
                       
                   end
                    newKcat=3600*rxn_kcat; %convert per second to per hour
                    rxns=sprintf('%s + %.15f X%d%s',rxns, 1/newKcat, rxns_by_this_complex(r),sep1); 
                    report_reaction(rxns_by_this_complex(r)) =1;
                    
                    if mod(r,150)==0
                        sep1=char(10);
                    else
                        sep1='';
                    end
              end
                %print kcat condition
                %print kcat condition
                glcouse_transporter_genes={'YDL245C' 'YDL247W' 'YDR342C' 'YDR343C' 'YDR345C' 'YDR536W' 'YEL069C' 'YFL011W' 'YHR092C' 'YHR094C' 'YHR096C' 'YJL214W' 'YJL219W' 'YJR158W' 'YJR160C' 'YLR081W' 'YMR011W' 'YNR072W' 'YOL156W'};
                saturation_transporter=1;
                if numel(find(ismember(complex_list(i),Protein_plasma(:,1))))>=1
                    
                   if numel(find(ismember(glcouse_transporter_genes, cell2mat(complex_list(i)))))>0
                       if mu<=0.27
                        saturation_transporter=1;%0.28*mu + 0.1355;%0.75/sigma;
                       else
                        saturation_transporter=1;% 2.3027*mu - 0.4047;%0.75/sigma;
                       end

                           
                   end
                    
                   
                    
                end
                
               fprintf(fptr,'CM%d: %s - %.15f X%d <= 0\n',k1,rxns,saturation_transporter*SS*weight_sigma,s); %K for average kcat or KK for different kcats, KK1 for GEKO
                   
                
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
                
                %store the sysnthesis reaction index
                Expression_syn(i)=s;
                
                
            end
        else % mito 
            fprintf(fptr_sfactor,'S%d\t%s\t%s\n',ks,complex_list{i},strrep(compartment{i},' ','_'));
            if ks<=numel(S_factor(:))
                    SS= S_factor(ks);
                    if SS==0
                        SS=sigma;
                    end
                else
                    SS=sigma;
                    
                end

               %SS=1;    
               if SS<=0.01
                 SS=0.01;
               end
              if SS>=1
                  SS=1;
              end
              
              
               ks=ks+1;
               
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
                
              rxn_kcat=cell2mat(KcatGEKO(i+1,3));
              if rxn_kcat==0
                  rxn_kcat = kcat_for_missing;
              elseif rxn_kcat<=1
                  rxn_kcat=1;
                  
              end
               newKcat=3600*rxn_kcat; %convert per second to per hour
                          
              fprintf(fptr,'CM%d: %.15f X%d - %.15f X%d <= 0\n',k1,1/newKcat, i,SS*weight_sigma,s);
                    
                
                
               %degradation
               rxnID=sprintf('%s_complex_degradation_mitochondrion',c);
               rxnID=strrep(rxnID,'-','');
               d=find(ismember(model.rxns,rxnID));
               fprintf(fptr,'CMEMDEG%d: X%d - %.15f X%d = 0\n',k1,d,c2,s);
               
               %dilution
               rxnID=sprintf('%s_complex_dilution_mitochondrion',c);
               rxnID=strrep(rxnID,'-','');
               dil=find(ismember(model.rxns,rxnID));
               fprintf(fptr,'CMMEMDIL%d: X%d - %.15f X%d = 0\n',k1,dil,c1,s);
               %store expression sythensis reaction
                Expression_syn(i)=s;
                
                %molecular weight
               content = sprintf('%s + %.15f X%d%c',content, (1/mu)*MW/1000.0,s,sep);
               
               if strcmp(Engery_content,'')==1
                   Engery_content = sprintf('%.15f X%d%c', (1/mu)*MW/1000.0,s,sep);         
               else
                   Engery_content = sprintf('%s + %.15f X%d%c',Engery_content, (1/mu)*MW/1000.0,s,sep);
               end
                
        end
        
    end
end

%close fptr_sfactor
fclose(fptr_sfactor);

%%

%tRNA loading
%
k1=1;
for i=1:numel(model.rxns)
    if k1<=100
    name=regexp(cell2mat(model.rxns(i)),'_loading');
    name1=regexp(cell2mat(model.rxns(i)),'_loading_mitochondrion');
    if numel(name1)==0
        comp='cytoplasm';
    else
        comp='mitochondrion';
    end
    if numel(name)>=1 && report_reaction(i)==0 && length(cell2mat(model.grRules(i)))>=1
        %find_reactions
        rxns_by_this_complex=find(ismember(model.grRules,model.grRules{i}));
        %find only the metabolic index
        metabolic_rxn_index  = find(rxns_by_this_complex>=numel(complex_list));
        rxns_by_this_complex= rxns_by_this_complex(metabolic_rxn_index);
       
        
        false=ones(numel(rxns_by_this_complex),1);
        for r=1:numel(rxns_by_this_complex)
            if strcmp(comp,'cytoplasm')==1
                name1=regexp(cell2mat(model.rxns(rxns_by_this_complex(r))),'_loading_mitochondrion');
                if numel(name1)>=1
                    false(r)=0; % remove this reaction
                end
            else %now the loading is in mitochondrion
                name1=regexp(cell2mat(model.rxns(rxns_by_this_complex(r))),'_loading_mitochondrion');
                if numel(name1)==0
                    false(r)=0; % remove this reaction
                end
            end
        end
        rxns_by_this_complex(false==0)=[];
        
        newKcat=kcat_AARS; %convert per second to per hour
        rxns=sprintf('%.15f X%d', 1/newKcat, rxns_by_this_complex(1));
        
        %rxns=sprintf('X%d',rxns_by_this_complex(1));
        report_reaction(rxns_by_this_complex(1)) =1;
        sep1='';
        for r=2:numel(rxns_by_this_complex)
            %%kcat for this reaction
            newKcat=kcat_AARS; %convert per second to per hour
            rxns=sprintf('%s + %.15f X%d%s',rxns, 1/newKcat, rxns_by_this_complex(r),sep1);
            report_reaction(rxns_by_this_complex(r)) =1;
            

        end
        
               
         c=cell2mat(strrep(MyModel.grRules(i),':','_'));
         rxnID=sprintf('%s_complex_formation_%s',c,comp);
         rxnID=strrep(rxnID,'-','');
         s=find(ismember(model.rxns,rxnID));
         fprintf(fptr,'CMLoading%d: %s - %.15f X%d <= 0\n',k1,rxns,weight_sigma,s); %K for average kcat or KK for different kcats, KK1 for GEKO

         

         rxnID=sprintf('%s_complex_degradation_%s',c,comp);
         rxnID=strrep(rxnID,'-','');
         d=find(ismember(model.rxns,rxnID));
         fprintf(fptr,'CMloading%d: X%d - %.15f X%d = 0\n',k1+1,d(1),c2,s);

         rxnID=sprintf('%s_complex_dilution_%s',c,comp);
         rxnID=strrep(rxnID,'-','');
         d=find(ismember(model.rxns,rxnID));
         fprintf(fptr,'CMloading%d: X%d - %.15f X%d = 0\n',k1+2,d(1),c1,s);
         
               k1=k1+1;
     
    end
    

    end
end


%% Ribosome capacity constrain
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
        fprintf(fptr,'CTranslation%d: X%d - %.15f R%d =0\n',r,i,kcat_ribo/(L+L_scanning + 4),r);
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
 

fprintf(fptr,'CR%d: %s - %.15f X%d <= 0\n',1,ve,ribosome_active_ratio/(mu+kdeg), syn_ribo );
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
 
 
%% Mitocondiral ribosome capacity constrain
% add mitocondiral ribosome capacity constrain
rxnID=sprintf('MitoRibosome_biosynthesis');
syn_ribo=find(ismember(model.rxns,rxnID));

rxnID=sprintf('MitoRibosome_degradation');
syn_deg=find(ismember(model.rxns,rxnID));

rxnID=sprintf('MitoRibosome_dilution');
syn_dil=find(ismember(model.rxns,rxnID));

%add set of proteins biosyenthesis reactions
ve_mito='';
r=0;
for i=1:numel(model.rxns)
    name1=regexp(cell2mat(model.rxns(i)),'_translation_mitochondrion');
     if numel(name1)>0 
        %index of mito protein r
        r=r+1;
        %protein Length
        gene=strrep(model.rxns(i),'_translation_mitochondrion','');
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
        fprintf(fptr,'CTranslation%d: X%d - %.15f MitoRibosome%d =0\n',r,i,kcat_ribo_mito/L,r);
        
        
    if strcmp(ve_mito,'')
        ve_mito=sprintf('MitoRibosome%d',r);

       
    else
        ve_mito=sprintf('%s + MitoRibosome%d',ve_mito,r);
    end
      gene=strrep(model.rxns(i),'_translation_mitochondrion','');
      gene_index=find(ismember(Proteins_MW(:,1),gene));
      MW=cell2mat(Proteins_MW(gene_index,2))/1000/(mu + kdeg);
      
      %total mass
      ve_total=sprintf('%s + %f X%d',ve_total,MW,i);
      ve_total_molecules=sprintf('%s + X%d', ve_total_molecules,i);
    end
end

fprintf(fptr,'CMR%d: %s - %.15f X%d <= 0\n',1,ve_mito,Mito_ribosome_active_ratio/(mu+kdeg), syn_ribo );
fprintf(fptr,'CMR%d: X%d - %.15f X%d = 0\n',2,syn_deg,c_ribo_deg_mito,syn_ribo );
fprintf(fptr,'CMR%d: X%d - %.15f X%d = 0\n',3,syn_dil,c_ribo_dil_mito,syn_ribo );

%% protein content constrains
%protein content constrains
fprintf(fptr,'CtotalProtein: %s = %f\n',ve_total,proteinRatio(mu)/(1-UnFoldRatio*(1-RefoldingRatio)));
total_number=molecule_number/(13*6.023e8)*(mu+kdeg);
fprintf(fptr,'CveTotal: %s <= %.15f\n',ve_total_molecules, total_number);

%%


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

%% Chaperone constraint
%Chaperone constrain
%misFoldingConstraint(MyModel,fptr,UnFoldRatio);
%%misFoldingYFPHsp10444(model,fptr,UnFoldRatio,mu,5,kdeg,sigmaTemp,RefoldingRatio);

% SSB
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'SSB');

% RAC
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'RAC');

% SSA
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'SSA');

% SSE
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'SSE');

% HSP40
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'HSP40');

% HSP90
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'HSP90');

% % PFD
 kcat = (1/1.5)*60*sigmaTemp;
 kdeg = kdeg;
 c1=kcat / (mu + kdeg);
 c2=mu / (mu + kdeg);
 c3=kdeg / (mu + kdeg);
% 
 addChaperoneConstrain(fptr,model,c1,c2,c3,'PFD');

% CCT
kcat = (1/1.5)*60*sigmaTemp;
kdeg = kdeg;
c1=kcat / (mu + kdeg);
c2=mu / (mu + kdeg);
c3=kdeg / (mu + kdeg);

addChaperoneConstrain(fptr,model,c1,c2,c3,'CCT');



%add content constrains
%fprintf(fptr,'CCMP1: %s + = 0.25\n',content );
%mass_content=sprintf('%s + %.15f X%d \n',content,1347851.15/(1000.0*(mu+kdeg_ribo)),syn_ribo);

%add the boundary

%% Space constraints
fprintf(fptr,'CMIM: %s <= %f\n',membraneCondition(model, mu,kdeg),c_mim);
fprintf(fptr,'PMR:   %s <= %f\n',plasmaCondition(model, mu,kdeg),plasma_area);
%fprintf(fptr,'PMR1   %s <= %f\n', outMembraneConstraints(),plasma_budjet);
%NH4 Vmax
% fprintf(fptr,'PMR2: X4172 <= %f\n',Vmax_NH4);
% fprintf(fptr,'PMR3: X4631 <= 15\n');
% fprintf(fptr,'PMR2: X4681 <= 15\n');

fprintf(fptr,'MitochondrionVolume: %s <= %f\n',mitoVolumeCondition(model, mu,kdeg),mitochondrionVolume);


%% Cytosolic Constrain
%sc=cytosolCondition(model, mu,kdeg);


%soplex needs 8109 characters in the line. We devide the condition into
%two lines.
   %[constrain1 constrain2 constrain3 constrain4]=cytosolConditionAll(model, mu,kdeg);
   [constrain1 constrain2 constrain3 ]=cytosolCondition(model, mu,kdeg);
   fprintf(fptr,'CCV1: %s  - XXCV =  0  \n',constrain1);
   fprintf(fptr,'CCV2: %s  + XXCV - XXCV2  =  0  \n',constrain2);
   fprintf(fptr,'CCV4: %s  + XXCV2 <= %f\n',constrain3,cyto_volume);



%% Ribosome Assembly constraints
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

%% proteasome constraint
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

%% un-modeled proteins
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

%% TOM constraints
%protein import in mitochondrion
kcat = (1/20)*60*60;
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

%% over-expression constraint, such as YFP or GFP proteins
%overexpression constraint, such as YFP or GFP proteins
    e=copies/(13*6.023e8);
    MMu=mu*e;
    KKdeg= kdeg*e;
    Syn=(mu+kdeg)*e;
    
     fprintf(fptr,'Overexpreesed1: X%d = %.15f\n',find(ismember(model.rxns, overExpressed.syn)),Syn);
     fprintf(fptr,'Overexpreesed2: X%d = %.15f\n',find(ismember(model.rxns, overExpressed.dilultion)),MMu);
     fprintf(fptr,'Overexpreesed3: X%d = %.15f\n',find(ismember(model.rxns, overExpressed.degradation)),KKdeg);

%     model.ub(find(ismember(model.rxns,overExpressed.dilultion)))=MMu;
%     model.lb(find(ismember(model.rxns,overExpressed.dilultion)))=MMu;
%     model.ub(find(ismember(model.rxns,overExpressed.degradation)))=KKdeg;
%     model.lb(find(ismember(model.rxns,overExpressed.degradation)))=KKdeg;
% 

%
%Print the bounds for each flux.
fprintf(fptr,'Bounds\n');
one_copy_per_cell=(mu+kdeg)*(1/(13*6.023e8));
N_cpc=100000000;
kk=1;
for i=1:numel(model.rxns)
           
            if model.ub(i)>=100
                fprintf(fptr,'%f <= X%d <= +infinity\n',model.lb(i),i);
            else
                fprintf(fptr,'%f <= X%d <= %f\n',model.lb(i),i,model.ub(i));
            end
      
end
fprintf(fptr,'End\n');
fclose(fptr);
pause(10);
%



%fprintf('finished\n');

system(command,'-echo');
%system(sprintf('./myResults1 %s',strrep(fileName,'.lp','')),'-echo');
%save('MyModelNew.mat', 'model');

fileNameOutput=sprintf('%s.out',fileName);

if strcmp(SoplexVersion,'Soplex-1.6')==1
   sol=readSoplex_results(fileNameOutput,model);
end
%soplex 3.0
if strcmp(SoplexVersion,'Soplex-3.0')==1
  sol=readSoplex3_results(fileNameOutput,model);
end
