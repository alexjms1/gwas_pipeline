function generatePed(file,columns)
% creates test .ped files for each strain specified in 'columns', from a 
% file 'file' which contains all strains (and converts strain names with a '/'
% to a strain name with a '.'). Also outputs .log files for each strain, 
% cotaining information about total number of SNPs and total number of
% missing SNPs ('N' -> '0' values).  May also be conducted using plink
% directly (see extractStrainsGenotypes.m).
numStrains=length(columns);
columnFormat=repmat('%s',1,numStrains);
fileID = fopen(file);
strains = textscan(fileID,columnFormat,1,'Delimiter','\t');
for i=1:numStrains
    strains{i}{1} = strrep(strains{i}{1},'/','.');
end
columnFormat=repmat('%c',1,numStrains);
genotypes=textscan(fileID,columnFormat,'HeaderLines',1,'Delimiter','\t');
fclose(fileID);
for i=1:numStrains
    missing=genotypes{i}=='N';
    genotypes{i}(missing)='0';
    numSNPs=length(genotypes{i});
    fileID=fopen(['C:\FastLMM\Plink\' strains{i}{1} '.test.ped'],'w');
    fprintf(fileID,'%s %s 0 0 0 -9 ',strains{i}{1},strains{i}{1});
    for j=1:numSNPs
        fprintf(fileID,'%c %c ',genotypes{i}(j), genotypes{i}(j));
    end
    fclose(fileID);
    fileID=fopen(['C:\FastLMM\Plink\' strains{i}{1} '.test.log'],'w');
    fprintf(fileID,'SNPs: %u , Missing: %u',numSNPs,missing);
    fclose(fileID);
end

end

