clear SGD_names 
clear rxn_complex
clear kcat

[a b ProteinNames]=xlsread('kcatmap.xlsx','uniprot_map');
[a b ProteinUniprot]=xlsread('kcatmap.xlsx','uniprot');
[a b kcat]=xlsread('kcatmap.xlsx','kcat');

[a b reactions]=xlsread('YeastNewThermaFinish3.xlsx','RXNS');
SGD_names={''};
rxn_complex={''};
rxn_kact=[];
r=1;
for i=1:numel(ProteinUniprot(:,1))

    for j=1:20
        if strcmp(cell2mat(ProteinUniprot(i,j)),'[]')==0
            p=regexp(cell2mat(ProteinUniprot(i,j)),' ','split');
            yp={''};
            yp_complex='';
            for k=1:numel(p)
                %map p to yp
                index=find(ismember(ProteinNames(:,1),p(k)));
                if numel(index)==0
                    error(sprintf('%s is not mapped',cell2mat(p(k))));
                else
                    yp(k)=ProteinNames(index,2);
                    if k==1
                        yp_complex=sprintf('%s',cell2mat(yp(k)));
                    else
                        yp_complex=sprintf('%s:%s',yp_complex,cell2mat(yp(k)));
                    end
                  
                end
            end
            SGD_names(i,j)={yp_complex};
            rxn_complex(r,1)={yp_complex};
            rxn_kcat(r,1)=cell2mat(kcat(i,j))/3600;
            index=regexp(cell2mat(reactions(r+2,2)),'_reverse');
            if numel(index)>0
                r=r+1;
                if r==3630 || r==3631 ||  r== 3607 || r== 3606 r== 3625 || r== 3635 || r== 3637
                     r=r+1;
                     rxn_complex(r,1)={''};
                     rxn_kcat(r,1)=0;
                else
                rxn_complex(r,1)={yp_complex};
                rxn_kcat(r,1)=cell2mat(kcat(i,j))/3600;
                  end
            end
           r=r+1; 
           if r==3630 || r==3631 ||  r== 3607 || r== 3606 || r== 3625 || r== 3635 || r== 3637
               r=r+1;
           end
        else
         if j==1
            index=regexp(cell2mat(reactions(r+2,2)),'_reverse');
            if numel(index)>0
                r=r+1;
                if r==3630 || r==3631 ||  r== 3607 || r== 3606 r== 3625 || r== 3635 || r== 3636
                     r=r+1;
                end
            end
             r=r+1;
             if r==3630 || r==3631 ||  r== 3607 || r== 3606 r== 3625  || r== 3635 || r== 3636
                 r=r+1;
             end
         end
         
        if r==101
            y=0;
        end
        end
      
    end
    end
