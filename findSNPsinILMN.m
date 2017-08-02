function findSNPsinILMN( root, chr, threshold, minPos, maxPos, cov)
% This function is used for extracting a set of SNPs that meet a q-value
% cutoff threshold within a certain chromosomal position range; useful for
% eQTL analysis, where a prior GWAS on expression data (e.g., ILMN/Illumina
% IDs) has resulted in results for a number of expression phenotypes - this
% function will scan across all these results and pull out SNPs that meet
% the cutoff within the specified chromosomal range.  Range is given in
% terms of base pairs.
% root = dir for merged results within "output" dir ('Output\' root
% '\Merged\' root '.' [your phenotypes] '.all.fastlmm.txt')
% chr = chromosome we're looking for the eQTL on
% threshold = q value cutoff
% minPos = minimum chromosomal position
% maxPos = maximal chromosomal position
% cov = whether covariates were used in FaST-LMM analyses of transcripts
%
% The function will output its results to 'Output\' root '\Merged\sigSNPs\sigSNPs_' [your phenotypes] '.txt'
% if there were SNPs that met the critera; otherwise, it will output files
% 'Output\' root '\Merged\sigSNPs\no_sigSNPs_' [your phenotypes] '.txt'.
if nargin < 6 
    error ('Not enough arguments (needs 6)');
end
format LONGG;
Qvals = dir(['Output\' root '\Merged\*LowestPsQs.txt']);
if size(Qvals,1)>1 
    error ('More than 1 LowestPsQs.txt found');
end
Qvals =Qvals(1).name;
Qvals = ['Output\' root '\merged\' Qvals];
fileID = fopen(Qvals);
lowestPsQs = textscan(fileID,'%s %*f %*f %f %*[^\n]','HeaderLines',1,'delimiter','\t');
ignore = lowestPsQs{2}(:)<=threshold;
lowestPsQs{1} = lowestPsQs{1}(ignore);
lowestPsQs{2} = lowestPsQs{2}(ignore);

fclose(fileID);

for i=1:size(lowestPsQs{1},1)
    if exist(['Output\' root '\Merged\' root '.' lowestPsQs{1}{i} '.all.fastlmm.txt'],'file')
        fileID=fopen(['Output\' root '\Merged\' root '.' lowestPsQs{1}{i} '.all.fastlmm.txt']);
        if ~cov
            columnFormat = repmat('%s',1,28);
        else
            columnFormat = repmat('%s',1,30);
        end
        dataHeaders = textscan(fileID,columnFormat,1,'Delimiter','\t');
        fclose(fileID);
        columnFormat = '%.0u %s  %*u %u %s %.10f %*f %.15f %.15f %.15f %.0u16 %*[^\n]';
        fileID=fopen(['Output\' root '\Merged\' root '.' lowestPsQs{1}{i} '.all.fastlmm.txt']);
        data = textscan(fileID,columnFormat,'HeaderLines',1,'Delimiter','\t');
        fclose(fileID);
        chrLogical = data{1}==chr;
        for j=1:size(data,2)
            data{j}=data{j}(chrLogical);
        end
       chrLogical = data{3} <= maxPos & data{3} >= minPos;
       for j=1:size(data,2)
            data{j}=data{j}(chrLogical);
       end
        chrLogical = data{7} <= threshold;
        for j=1:size(data,2)
            data{j}=data{j}(chrLogical);
        end
        outSize = size(data{1},1);
        if outSize > 0
            [data{5}, idx] = sort(data{5},1);
            for j=1:size(data,2)
                if j~=5
                    data{j}=data{j}(idx);
                end
            end
        end
        
        outSize = size(data{1},1);  
        if outSize > 0 
            if ~exist(['Output\' root '\Merged\sigSNPs\'], 'dir')
                mkdir(['Output\' root '\Merged\sigSNPs\']);
            end
            fileResults = fopen(['Output\' root '\Merged\sigSNPs\sigSNPs_' lowestPsQs{1}{i} '.txt'],'w');    
            fprintf(fileResults,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n',dataHeaders{1}{1}, dataHeaders{2}{1}, dataHeaders{4}{1}, dataHeaders{5}{1}, dataHeaders{6}{1}, dataHeaders{8}{1}, dataHeaders{9}{1}, dataHeaders{10}{1}, dataHeaders{11}{1});
            columnFormat = '%.0u\t %s\t %u\t %s\t %.10f\t %.15f\t %.15f\t %.15f\t %.0u \n';
            for j=1:outSize
                fprintf(fileResults,columnFormat,data{1}(j),data{2}{j},data{3}(j),data{4}{j},data{5}(j),data{6}(j),data{7}(j),data{8}(j),data{9}(j)); %
            end
        else
            fileResults = fopen(['Output\' root '\Merged\sigSNPs\no_sigSNPs_' lowestPsQs{1}{i} '.txt'],'w');
            fprintf(fileResults,'No SNPs meet criteria.\n');
        end
        fclose(fileResults);


    else
        error(['file Output\' root '\merged\' lowestPsQs{1}{i} '.all.fastlmm.txt does not exist!']); 

    end
end
end
