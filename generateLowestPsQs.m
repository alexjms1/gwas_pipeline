function generateLowestPsQs(root, filenames,testQvals,renamePreviousAppendOverwrite)
% Helper function for generateFullSummary, which (re-)creates '.LowestPsQs.txt'
% or '.LowestPs.txt' files (determined by the boolean testQvals) which are
% summary files containing lowest Ps / Qs and RS information for a set of
% phenotypes specified by filenames.  Outputted as 'Output\' root 
% '\Merged\regen_' root '.LowestPs.txt'.
% INPUTS: 
%   root: root directory for your FaST-LMM results files, such that results
%       are found in 'Output\' root ' Merged'.
%   filenames: phenotype file names of the format 'Output\' root '\Merged\'
%       root '.' filename '.all.fastlmm.txt'
%   testQvals: a boolean, specifying whether to test Q values.  Can take
%       the values 0-4 with positive values greater than 1 for consistency
%       with the mergeResults function; however, all values greater than 1
%       are treated the same here (i.e., only the most commonly used measures,
%       FDR, Q, and BHFR created).
%   renamePreviousAppendOverwrite: 0 = overwrite, 1 = append, 2 = copy
%       previous files (i.e., saving them with prefix 'BACKUP_') and make a new one
if nargin ~= 4
    error('Not enough arguments passed to function. 4 are required');
end
if testQvals > 3
    error('Q value test integer not within range of 0 - 4');
end
format LONGG;
numFiles=length(filenames);
lowestPs = zeros(1,numFiles);
if testQvals
    lowestStQs = zeros(1,numFiles);
    lowestStFDRs = zeros(1,numFiles);
    lowestBHFDRs = zeros(1,numFiles);
end 
lowestP_chr = zeros(1,numFiles);
lowestP_pos = zeros(1,numFiles);
lowestP_rs = cell(1,numFiles);
data=cell(numFiles,1);

if numFiles
    if renamePreviousAppendOverwrite == 2
        writePriv = 'w';
        oldFiles=dir(['Output\' root ' Merged\*LowestPs*.txt']);
        for i=1:length(oldFiles)
            movefile(['Output\' root ' Merged\' oldFiles(i).name],['Output\' root ' Merged\BACKUP_ ' strrep(oldFiles(i).name,'LowestPs','')]);
            
        end
    elseif renamePreviousAppendOverwrite == 1
        writePriv = 'a';
    elseif renamePreviousAppendOverwrite == 0
        writePriv = 'w';
    else
        error('Bad write priviledge entry; 0 for overwrite, 1 for append, 2 for copy previous and write a new file');
    end
    if ~testQvals
        columnFormat = '%u %s %u %f %s %.10f %*[^\n]';
    else
        columnFormat = '%u %s %u %f %s %.10f %f %f %.15f %.15f %.15f %*[^\n]'; 
    end
    for phenotype=1:numFiles
            fileID = fopen(['Output\' root '\Merged\' root '.' filenames{phenotype} '.all.fastlmm.txt']);
            if fileID<1
                display(['Error opening ' filenames{phenotype}]);
            else
                data{phenotype} = textscan(fileID,columnFormat,'Delimiter','\t','Headerlines',1);
                fclose(fileID);
            end
    end
    for phenotype=1:numFiles
        [lowestPs(phenotype), index]= min(data{phenotype}{6});
        lowestP_chr(phenotype) = data{phenotype}{1}(index);
        lowestP_rs{phenotype} = data{phenotype}{2}{index};
        lowestP_pos(phenotype) = data{phenotype}{4}(index);
        if testQvals
           [fdr, q] = mafdr(pvals);
           lowestStQs(phenotype)=min(q);
           lowestStFDRs(phenotype)=min(fdr);
           fdrBH = mafdr(pvals,'BHFDR',true);
           lowestBHFDRs(phenotype)=min(fdrBH);  
        end            
    end

    
          
    if testQvals
        fileID = fopen(['Output\' root '\Merged\regen_' root '.LowestPsQs.txt'],writePriv);  
        fprintf(fileID,'Phen\tPval\tpFDR\tQval\tBHFDR\tbestChrPos\t\tbestRS\n');
        columnFormat = '%s\t%.10f\t%.15f\t%.15f\t%.15f\t%.0u\t%u\t%s\n';
        for cnt=1:numFiles
            fprintf(fileID,columnFormat,filenames{cnt}, filenames{cnt}, lowestPs(cnt), lowestStFDRs(cnt), lowestStQs(cnt), lowestBHFDRs(cnt), lowestP_chr(cnt), lowestP_pos(cnt), lowestP_rs{cnt});
        end
    else
        fileID = fopen(['Output\' root '\Merged\regen_' root '.LowestPs.txt'],writePriv);  
         fprintf(fileID,'Phen\tPval\tbestChrPos\t\tbestRS\n');
        columnFormat = '%s\t%.10f\t%.0u\t%f\t%s\n';
        for cnt=1:numFiles
            fprintf(fileID,columnFormat,filenames{cnt}, lowestPs(cnt), lowestP_chr(cnt), lowestP_pos(cnt), lowestP_rs{cnt});
        end
    end 
    fclose(fileID);

else
    display('No filenames given.');
end
format;
end
