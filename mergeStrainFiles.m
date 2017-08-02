function mergeStrainFiles(source,genome,append,remake,retain_local_copy,delete_geno_copies)
% This function is for use after extractStrainsGenotypes, when .ped files
% for each individual strain have been made, along with a single allStrains.map.
% This function parses your phenotypes file, determining which are in your sample,
% and assembles copies of genomes as appropriate to represent your sample.
% Inserts subject number into .ped file as appropriate for FaST-LMM, generating
% a copy for each subject named according to starin/subject combinations.
% Merges genomes into a single binary plink format genome filtered such that 
% minor allele frequency > 0.05 and SNP missingness < .1 using
% allStrains.map (SNP map identical across strains).  Places resultant
% files in C:\FastLMM\Cpp_MKL\.
% INPUTS:
%   source: your phenotypes file, as specified by
%       'C:\FastLMM\Cpp_MKL\Files\' source '.txt'.  See example referenced in readme.
%   genome: a name for the genome file to be outputted, under 'C:\FastLMM\CPP_MKL\'
%   append: a boolean, specifying whether to append to an already existing
%       binary plink file 'merged.bim'
%   remake: a boolean, specifying whether to force-remake the individual
%       genome files
%   retain_local_copy: a boolean, specifying whether to retain the merged
%       binary plink file in C:\FastLMM\Plink (merged.bed/.bim/.fam)
%   delete_geno_copies: a boolean, specifying whether to retain the
%       individual strain/subject combination .map/.ped files 
%       in 'C:\FastLMM\Plink\.  Does not delete allStrains.map.
if nargin ~= 6
    error('Not enough arguments passed to generateLowestPsQs. 6 are required: filename, output genome name, append (true/false), remake (true/false), retain local copy in plink folder (true/false), delete copies of genomes made during merge process (true/false).');
end
if append && ~(exist('C:\FastLMM\Plink\merged.bim','file'))
    display('No file to merge to (merged.bim).');
    return;
elseif ~exist(['C:\FastLMM\Cpp_MKL\Files\' source '.txt'],'file')
    display(['Input data ' source '.txt not found']);
    return;
end
if remake || ~append
	if exist('C:\FastLMM\Plink\merged.bim','file')
        delete('C:\FastLMM\Plink\merged.bim');
	end
    if exist('C:\FastLMM\Plink\merged.bed','file')
        delete('C:\FastLMM\Plink\merged.bed');
    end
	if exist('C:\FastLMM\Plink\merged.fam','file')
        delete('C:\FastLMM\Plink\merged.fam');
	end 
end
fileID = fopen(['C:\FastLMM\Cpp_MKL\Files\' source '.txt']);
sample = textscan(fileID,'%s %s %*[^\n]','delimiter','\t');
fclose(fileID);
numSubjects = size(sample{1},1);
includeStrain = zeros(numSubjects,1);
startStrain = 0;
for i=2:numSubjects
    if ~exist(['C:\FastLMM\Plink\' sample{1}{i} sample{2}{i} '.ped'],'file')
        if exist(['C:\FastLMM\Plink\' sample{1}{i} '.ped'],'file')
            rep = replaceinfile([sample{1}{i} ' ' sample{1}{i}],[sample{1}{i} ' ' sample{2}{i}],['C:\FastLMM\Plink\' sample{1}{i} '.ped'],['C:\FastLMM\Plink\' sample{1}{i} sample{2}{i} '.ped']);
            if rep
                fileID = fopen(['C:\FastLMM\CPP_MKL\' genome '.log'],'a');
                fprintf(fileID,'%s genome error: could not insert subject number into ped file. Check file permissions.',sample{1}{i}); 
                cont = menu([sample{1}{i} ' was not found in database.\nWould you like to continue?'],'Yes','No');
                if cont == 2
                    fprintf(fileID,'\nUser aborted.');
                    fclose(fileID);
                    return
                end
                fclose(fileID);
            end
            includeStrain(i)=1;           
        else
            fileID = fopen(['C:\FastLMM\CPP_MKL\' genome '.log'],'a');
            fprintf(fileID,'%s genome was not found (check the Plink folder).',sample{1}{i}); 
            cont = menu([sample{1}{i} ' was not found in database.\nWould you like to continue?'],'Yes','No');
            if cont == 2
                 fprintf(fileID,'\nUser aborted.');
                 fclose(fileID);
                 return
            end
             fclose(fileID);
        end
    elseif remake
        if ~exist(['C:\FastLMM\Plink\' sample{1}{i} '.ped'],'file')
            fileID = fopen(['C:\FastLMM\CPP_MKL\' genome '.log'],'a');
            fprintf(fileID,'%s genome was not found (check the Plink folder).',sample{1}{i}); 
            cont = menu([sample{1}{i} ' was not found in database.\nWould you like to continue?'],'Yes','No');
            if cont == 2
                 fprintf(fileID,'\nUser aborted.');
                 fclose(fileID);
                 return
            end
             fclose(fileID);
        else
            includeStrain(i)=1;
        end
    end
end

if sum(includeStrain) ~= numSubjects-1
	fileID = fopen(['C:\FastLMM\CPP_MKL\' genome '.log'],'a');
	fprintf(fileID,'Number of strains detected not equal to the number of strains in the phenotype file (%u of %u).\n', sum(includeStrain),numSubjects-1);
    cont = menu(['Not all strains were found - ' sum(includeStrain) ' of ' numSubjects-1 ' found..\nContinue?  Otherwise, check strain names and case sensitivity.'],'Yes','No');
    if cont == 2

        fprintf(fileID,'User aborted - check strain names and case sensitivity.');
        fclose(fileID);
        return
    end
    fclose(fileID);
end
for i=2:numSubjects
    if includeStrain(i)
        startStrain = i;
        if ~(exist(['C:\FastLMM\Plink\' sample{1}{i} sample{2}{i} '.map'],'file')) 
            copyfile('C:\FastLMM\Plink\allStrains.map',['C:\FastLMM\Plink\' sample{1}{i} sample{2}{i} '.map']);
        end 
        break
    end
end
if i == 0
    disp('No subjects to merge were found.');
    return
end
if startStrain
    fileID=fopen('C:\FastLMM\Plink\strainMerge.txt','w');
    for i=startStrain+1:numSubjects-1
        if includeStrain(i)
           fprintf(fileID, '%s%s.ped allStrains.map\n', sample{1}{i}, sample{2}{i});
        end
    end
    if includeStrain(numSubjects) &&  ~(numSubjects==startStrain)
         fprintf(fileID, '%s%s.ped allStrains.map', sample{1}{numSubjects}, sample{2}{numSubjects}); 
    end
    fclose(fileID);
    if ~append
        command = ['cd C:\FastLMM\Plink\ & plink --file ' sample{1}{startStrain} sample{2}{startStrain}];
        if ~(numSubjects==startStrain)
            command = [command ' --merge-list strainMerge.txt --make-bed --out merged'];
        else
            command = [command ' --make-bed --out merged'];
        end
        dos(command,'-echo');
    else
        command = ['cd C:\FastLMM\Plink\ & plink --file ' sample{1}{startStrain} sample{2}{startStrain}];
        if ~(numSubjects==startStrain)
            command = [command ' --merge-list strainMerge.txt --make-bed --out toBeMerged'];
        else
            command = [command ' --make-bed --out toBeMerged'];
        end
        dos(command,'-echo');
        command = 'cd C:\FastLMM\Plink\ & plink --bfile merged --bmerge toBeMerged.bed toBeMerged.bim toBeMerged.fam --make-bed --out mergeAppended';
        dos(command,'-echo');
        delete('toBeMerged.bed');
        delete('toBeMerged.bim');
        delete('toBeMerged.fam');
         if exist('C:\FastLMM\Plink\merged.bim','file')
             delete('C:\FastLMM\Plink\merged.bim');
         end
         if exist('C:\FastLMM\Plink\merged.bed','file')
             delete('C:\FastLMM\Plink\merged.bed');
         end
         if exist('C:\FastLMM\Plink\merged.fam','file')
             delete('C:\FastLMM\Plink\merged.fam');
         end
        dos('ren C:\FastLMM\Plink\mergeAppended.bed merged.bed','-echo')
        dos('ren C:\FastLMM\Plink\mergeAppended.bim merged.bim','-echo')
        dos('ren C:\FastLMM\Plink\mergeAppended.fam merged.fam','-echo')
    end
    command = ['cd C:\FastLMM\Plink\ & plink --bfile merged --maf 0.05 --geno 0.1 --make-bed --out ' genome];
    dos(command,'-echo'); % '--write-snplist'
    command = ['for %f in (C:\FastLMM\Plink\' genome '.bed, C:\FastLMM\Plink\' genome '.bim, C:\FastLMM\Plink\' genome '.fam, C:\FastLMM\Plink\' genome '.log, C:\FastLMM\Plink\' genome '.nosex) do move /Y %f C:\FastLMM\Cpp_MKL\'];
    dos(command,'-echo');
    if exist('C:\FastLMM\Plink\strainMerge.txt','file')   
        delete('C:\FastLMM\Plink\strainMerge.txt');
    else
        display('Warning ... no strainMerge.txt found for deletion');
    end
    if ~retain_local_copy
       if exist('C:\FastLMM\Plink\merged.bim','file')
           delete('C:\FastLMM\Plink\merged.bim');
        end
        if exist('C:\FastLMM\Plink\merged.bed','file')
           delete('C:\FastLMM\Plink\merged.bed');
        end
        if exist('C:\FastLMM\Plink\merged.fam','file')
            delete('C:\FastLMM\Plink\merged.fam');
        end
    end
else
    display('No strains to merge or all strains already have plink copies (if so use remake = 1).')
end
if delete_geno_copies
    for i=2:numSubjects
        if exist(['C:\FastLMM\Plink\' sample{1}{i} sample{2}{i} '.ped'],'file') && ~strcmpi(sample{1}{i},sample{2}{i}) && ~isempty(sample{2}{i})
            delete(['C:\FastLMM\Plink\' sample{1}{i} sample{2}{i} '.ped']);
        end
    end
    mapStart = ['C:\FastLMM\Plink\' sample{1}{startStrain} sample{2}{startStrain} '.map'];
    if exist(mapStart,'file') && ~strcmp(mapStart,'C:\FastLMM\Plink\allStrains.map')
        delete(mapStart);
    end
end
sprintf('\nMerge script finished.')
end         