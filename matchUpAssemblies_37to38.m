%This script demonstrates an effort to match SNP data from MouseRef37
%assembly to MouseRef38.  A work in progress.

fileRe=fopen('C:\FastLMM\Plink\remapped_Test_trimmed.txt');
file37=fopen('C:\FastLMM\Plink\masterNewGenos.txt');
%fileReference=fopen('C:\FastLMM\Plink\Reference.txt');

remapped=textscan(fileRe,'%u %u %u %u','delimiter','\t');
headers37=textscan(file37,repmat('%s',1,18),1,'delimiter','\t');
data37=textscan(file37,'%*s %c %c %u %u %u %c %c %c %c %c %c %c %c %c %c %c %c','HeaderLines',1,'delimiter','\t');
%reference=textscan(fileReference,'%u %u %u %u','delimiter','\t');
fclose(fileRe);
fclose(file37);
%fclose(fileReference);

tempRe=remapped;
remove=(remapped{4}==0);
tempRe{4}(remove)=[];
tempRe{3}(remove)=[];
tempRe{2}(remove)=[];
tempRe{1}(remove)=[];
temp37=data37;
byChromRe=cell(21,1);
byChrom37=cell(21,1);
for i=1:19
    byChromRe{i}=tempRe{3}(tempRe{1}==i);
    byChrom37{i}=data37{5}(data37{3}==i);
end
byChromRe{20}=tempRe{3}(tempRe{1}==23);
byChrom37{20}=data37{5}(data37{3}==23);
byChromRe{21}=tempRe{3}(tempRe{1}==24);
byChrom37{21}=data37{5}(data37{3}==24);
match=cell(21,1);
match_index37=cell(21,1);
match_indexRe=cell(21,1);
for i=1:21
    [match{i}, match_index37{i}, match_indexRe{i}]=intersect(byChrom37{i},byChromRe{i},'stable');
end
size=0;
for i=1:length(match)
    size=size+length(match{i});
end
%size
%del_index37=cell(21,1);
%for i=1:21
%   
%    del_index37{i}
%size=0;


%for i=1:21
    
%    take from 2 in Re
%    4-chrom37

data37Chrom=cell(21,1);
remapChrom=cell(21,1);
fields37=length(data37);
for i=1:19
    remapChrom{i}=tempRe{2}(tempRe{1}==i);
    data37Chrom{i}=cell(fields37,1);
    for j=1:fields37
        data37Chrom{i}{j}=data37{j}(data37{3}==i);
    end
end
remapChrom{20}=tempRe{2}(tempRe{1}==23);
data37Chrom{20}=cell(fields37,1);
for j=1:fields37
    data37Chrom{20}{j}=data37{j}(data37{3}==23);
end
data37Chrom{21}=cell(fields37,1);
remapChrom{21}=tempRe{2}(tempRe{1}==24);
for j=1:fields37
    data37Chrom{21}{j}=data37{j}(data37{3}==24);
end

notMatch=cell(21,1);
for i=1:21
    notMatch{i}=(1:length(data37Chrom{i}{1}))';
    notMatch{i}(match_index37{i})=[];
end
for i=1:21
    data37Chrom{i}{4}(match_index37{i})=remapChrom{i}(match_indexRe{i});
    for j=1:length(data37)
        data37Chrom{i}{j}(notMatch{i})=[];
    end
end
writeGeno=fopen('C:\FaSTLMM\Plink\translatedGenotypes.txt','w');
writeLines=fopen('C:\FaSTLMM\Plink\matchedLinesSNPs.txt','w');  
fprintf(writeGeno,repmat('%s\t',1,17), headers37{3}{1}, headers37{4}{1}, headers37{5}{1}, headers37{6}{1}, headers37{7}{1}, headers37{8}{1}, headers37{9}{1}, headers37{10}{1}, headers37{11}{1}, headers37{12}{1}, headers37{13}{1}, headers37{14}{1}, headers37{15}{1}, headers37{16}{1}, headers37{17}{1}, headers37{18}{1});
fprintf(writeLines,'Match Line\tIndex from build 37\tIndex from build 38\tOld GenoPos\tNew GenoPos\n');
for i=1:21
    %for j=1:length(match{i})        
    %    fprintf(writeLines,'%u\t%u\t%u\t%u\n',match{i}(j), match_index37{i}(j), match_indexRe{i}(j),byChrom37{4}(match_index37{i}(j)),(data37Chrom{i}{4}(j)));
   % end
    for j=1:length(data37Chrom{i}{1})
        fprintf(writeGeno,'\n%c\t%c\t%u\t%u\t%u\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%c\t%c', data37Chrom{i}{1}(j), data37Chrom{i}{2}(j), data37Chrom{i}{3}(j), data37Chrom{i}{4}(j), data37Chrom{i}{5}(j), data37Chrom{i}{6}(j), data37Chrom{i}{7}(j), data37Chrom{i}{8}(j), data37Chrom{i}{9}(j), data37Chrom{i}{10}(j), data37Chrom{i}{11}(j), data37Chrom{i}{12}(j), data37Chrom{i}{13}(j), data37Chrom{i}{14}(j), data37Chrom{i}{15}(j), data37Chrom{i}{16}(j), data37Chrom{i}{17}(j));
    end
end
fclose(writeLines);
fclose(writeGeno);
    