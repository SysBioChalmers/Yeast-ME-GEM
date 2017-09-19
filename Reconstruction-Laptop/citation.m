clear PaperID
clear Name Sections
[ a b citations]=xlsread('Citations.xlsx','scopus');
names=citations(:,4);

Name={''};
Sections={''};
PaperID={''};
k=1;
for i=2:numel(names)
    n=regexp(cell2mat(names(i)),'\.\, ','split');
    for j=1:numel(n)
        nn=cell2mat(n(j));
        if numel(n)==j
            nn=strrep(nn,'.','');
        end
        Name{k,1}=nn;
        Sections{k,1}=cell2mat(citations(i,3));
        PaperID{k,1}=cell2mat(citations(i,19));
        k=k+1;
    end
    
end

%buid the Graph Matrix
[ a b Section]=xlsread('Citations.xlsx','Sections');
PI=Section(:,2);
Sections=Section(:,1);
GRAPH=zeros(17,17);
paper=unique(PaperID);
for i=1:numel(paper)
    % find the author
    index=find(ismember(PaperID,paper(i)));
    Author=Name(index);
    Index=[];
    k=1;
    for j=1:numel(Author)
        for p=1:numel(PI)
            aa=regexp(cell2mat(Author(j)),cell2mat(PI(p)));
            if numel(aa)>0
                Index(k)=p;
                k=k+1;
            end
            
        end
    end
    
    %count the papers
    Index=unique(Index);
    
    for j=1:numel(Index)
        
        GRAPH(Index(j),Index(j))=GRAPH(Index(j),Index(j)) + 1;
        for k=j+1:numel(Index)
            GRAPH(Index(j),Index(k))=GRAPH(Index(j),Index(k)) + 1;
            GRAPH(Index(k),Index(j))=GRAPH(Index(k),Index(j)) + 1;
        end
    end
        
    end
    
 
fid=fopen('Citation.xgmml','w');
fprintf(fid,'<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n');
fprintf(fid,'<graph label="Reporter Subnetwork" \n');
fprintf(fid,'xmlns:xlink="http://www.w3.org/1999/xlink" \n');
fprintf(fid,'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"\n');
fprintf(fid,'xmlns:cy="http://www.cytoscape.org"\n');
fprintf(fid,'xmlns="http://www.cs.rpi.edu/XGMML" \n');
fprintf(fid,'directed="1">\n');

 for i=1:numel(Sections)
     fprintf(fid,'<node label="%s" id="%s">\n',strrep(cell2mat(Section(i)),' ','_'),strrep(cell2mat(Section(i)),' ','_'));
     fprintf(fid,'<graphics type="ELLIPSE" h="%d" w="%d" fill="#936c90" width="0" outline="#000000" transparency="255"/>\n',GRAPH(i,i)+5,GRAPH(i,i)+5);
     fprintf(fid,'</node>\n');
     

 end
 
 
 k=1;
 for i=1:numel(Sections)
     for j=i+1:numel(Sections)
         if GRAPH(i,j)~=0
             fprintf(fid,'<edge label="%d" source="%s" target="%s">\n',k,strrep(cell2mat(Section(i)),' ','_'),strrep(cell2mat(Section(j)),' ','_'));
             fprintf(fid,'<graphics width="%d" fill="#1f1f1f" transparency="50"/>\n',GRAPH(i,j)+1);
             fprintf(fid,'</edge>\n');
             k=k+1;
         end
     end
 end

 fprintf(fid,'</graph>\n');
fclose(fid);


 

 
 
