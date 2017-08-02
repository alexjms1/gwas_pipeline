function extractStrainsGenotypes(source)
% Uses an input genome file containing all strains (here, transposed plink format) to output
% ped/map files, one for each strain.  Uses parallel processing.
% Note directories hard coded to 'C:\FastLMM\Plink'
if nargin ~= 1
    error('Incorrect number of arguments passed to function. 1 is required.');
end
% Creates allStrains.map/ped from the transposed plink format
dos(['cd C:\FastLMM\Plink\ & plink --tfile ' source ' --recode --out allStrains'],'-echo');
fileID = fopen(['C:\FastLMM\Plink\' source '.tfam']);
refSubjects = textscan(fileID,'%s %*[^\n]','delimiter','\t');
refSubjects = refSubjects{1};
fclose(fileID);
numRefStrains = length(refSubjects);
% Iterates through all strains, creating a file,
% strainIsolator_[strain].txt for each, used for isolating a particular
% strain for extraction from allStrains.map/ped; these files are not kept
parfor i=1:numRefStrains
    fileID=fopen(['C:\FastLMM\Plink\strainIsolator_' refSubjects{i} '.txt'],'w');
    fprintf(fileID, '%s\t%s', refSubjects{i}, refSubjects{i});
    fclose(fileID);
    command = ['cd C:\FastLMM\Plink\ & plink --file allStrains --keep strainIsolator_' refSubjects{i} '.txt --recode --out ' refSubjects{i}];
    dos(command,'-echo');
    delete(['C:\FastLMM\Plink\strainIsolator_' refSubjects{i} '.txt']);
end
mapFiles = dir('C:\FastLMM\Plink\*.map');
% Deletes all .map files for individual strains (leaving allStrains.map)
% since they are genotyped at the same SNP map, and are therefore redundant
for i=1:length(mapFiles)
    if ~strcmp(mapFiles(i).name,'allStrains.map')
        delete(['C:\FastLMM\Plink\' mapFiles(i).name]);
    end
end
        
end