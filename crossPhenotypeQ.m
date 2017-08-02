function crossPhenotypeQ(filenames, output)
%can calculate Q values based upon multiple phenotypes; usually only works
%with a few phenotypes (liekly due to memory limitations)
%INPUTS:
%   filenames, a cell array of strings ocntaining the filenames
%   corresponding to the multiple phenotypes (note the '.all.fast.lmm.txt'
%   suffix, implying these are meant to be files that have been merged
%   across all chromosomes  Directory is hard coded as 'Output\Merged'.
%   output, a filename for the merged results, will carry the suffix
%   'mergedPs.createQvals.txt' Directory hard coded as 'Output\Merged\MergedQs'.
if nargin ~= 2
    display('Not enough arguments passed to function. 2 are required');
return
end
if iscell(filenames)
    format longg;
    currFile = ['Output\Merged\' filenames{1} '.all.fastlmm.txt'];
    if exist(currFile,'file')
        fileID=fopen(['Output\Merged\' filenames{1} '.all.fastlmm.txt']);
        columnFormat = [repmat('%s',1,7) '%*[^\n]'];
        dataHeaders = textscan(fileID,columnFormat,1,'Delimiter','\t');
        fclose(fileID);
    else
        disp(['File ' currFile ' not found.']);
        return
    end
    columnFormat = '%u8 %s %u8 %f %s %.15f %.15f';
    data = cell(length(filenames),1);
    for phenotype=1:length(filenames) 
        currFile = ['Output\Merged\' filenames{phenotype} '.all.fastlmm.txt'];
        if exist(currFile,'file')
            fileID=fopen(currFile);
            data{phenotype} = textscan(fileID,columnFormat,'HeaderLines',1,'Delimiter','\t');
            fclose(fileID);
        else
            disp(['File ' currFile ' not found.']);
            return   
        end
    end
    %
    for phenotype=1:length(filenames)
        for columnNum=1:size(dataHeaders,2)
            if ~(columnNum==6)
                data{phenotype}{columnNum} = [];
            end
        end
    end
    mergedPs = [];
    %
    for phenotype=1:length(filenames)
        mergedPs = [mergedPs data{phenotype}{6}(:)'];
    end

    %[fdr, q] = mafdr(mergedPs);
    q = mafdr(mergedPs);
    fileID=fopen(['Output\Merged\MergedQs\' output '.mergedPs.createQvals.txt'],'w');
    for i=1:size(dataHeaders,2)
        fprintf(fileID,'%s\t',dataHeaders{i}{1}); 
    end
    additionalCols = cell(2,1);
    additionalCols{1} = 'all_pFDR';
    additionalCols{2} = 'all_Q_values';
    fprbintf(fileID,'%s\t%s\n',additionalCols{1},additionalCols{2});
    columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.15f\t%.15f\t%.15f\t%.15f\n';
    count = 0;
    for i=1:size(data,1)
        phenoFileID=fopen(['Output\Merged\MergedQs\' filenames{i} '.mergedPs.createQvals.txt'],'w');
        for j=1:size(dataHeaders,2)
            fprintf(phenoFileID,'%s\t',dataHeaders{j}{1}); 
        end
        fprintf(phenoFileID,'%s\t%s\n',additionalCols{1},additionalCols{2});
        for j=1:size(data{i}{1},1)
            count = count+1;
            fprintf(phenoFileID,columnFormat,data{i}{1}(j), data{i}{2}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count));
            fprintf(fileID,columnFormat,data{i}{1}(j), data{i}{2}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count));
        end
        fclose(phenoFileID);
    end
    fclose(fileID);
    format short;
else
    disp('Input must be a cell array containing string(s).');
end
end
