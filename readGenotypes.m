function readGenotypes(path,mapPath)
% This function is an attempt at converting imputed genotypes from MouseRef
% assembly 37 to MouseRef 38, outputting summary statistics regarding the
% success rate of this conversion.  A work in progress.
format LONGG;
if path(end) == '\'
    path = path(1:end-1);
end
if mapPath(end) == '\'
    mapPath = mapPath(1:end-1);
end
if ~exist(path,'dir')
    disp([path ' is invalid.\n'])
    return
end
if ~exist(mapPath,'dir')
    disp([mapPath ' is invalid.\n'])
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fID = fopen([mapPath '\allStrains.map']);
map = textscan(fID,'%u8%*s%*u%u32%*[^\n]','Delimiter','\t');
fclose(fID);
mapByChrom = cell(20,1);
for i=1:19
    mapByChrom{i} = map{2}(map{1}==i);
end
mapByChrom{20} = map{2}(map{1}==23);
numMapSNPs = length(map{1});
clear map;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
remapFiles = dir([path '\report_*.txt']);
filenamesOnlyRemapped=cell(length(remapFiles),1);
parfor i=1:length(remapFiles)
   filenamesOnlyRemapped{i}=remapFiles(i).name;
end
filenamesOnlyRemapped = sortn(filenamesOnlyRemapped);
numRemappedFiles = length(filenamesOnlyRemapped);
%remapped = cell(numRemappedFiles,1);
%fID = cell{numRemappedFiles,1};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

genotypeFiles = dir([path '\*.csv']);
filenamesOnlyGenotypes=cell(length(genotypeFiles),1);
for i=1:length(genotypeFiles)
   filenamesOnlyGenotypes{i}=genotypeFiles(i).name;
end
%numGenotypeFiles = length(filenamesOnlyGenotypes);
%chrGenotypes=cell(numGenotypeFiles,1);
fID=fopen([path '\' filenamesOnlyGenotypes{1}]);
numFieldsGenotypes = 23;
fieldNamesGenotypes = textscan(fID,repmat('%s',1,numFieldsGenotypes),1,'Delimiter',',');
fclose(fID);
numGenomes = (length(chrGenotypes)-5)/2; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numRemappedFiles
    fID = fopen([path '\' filenamesOnlyRemapped{i}]);    
  %  remapped{i} = textscan(fID,'%*s%*u%*u%*s%s%*u%*u%u32%*u%*c%*u%*u%u32%*u%c%*u%*s%*s','Headerlines',1,'Delimiter','\t');
  %  remapped{i}{
    remapped = textscan(fID,'%*s%*s%*s%*s%s%*s%*s%s%*s%*s%*s%*s%s%*s%s%*s%*s%*s','Headerlines',1,'Delimiter','\t');
    fclose(fID);
    fID=fopen([path '\' filenamesOnlyGenotypes{i}]);
    chrGenotypes = textscan(fID,'%u%u32%*u%*s%*s%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8','Headerlines',1,'Delimiter',',');
    fclose(fID);
    if length(remapped{1}) ~= length(chrGenotyes{1})
        disp(['non correspondance for chr ' i]);
        return
    end
    cleaning1 = strcmp('NULL',remapped{1});
    cleaning2 = strcmp('-',remapped{4});
    cleaning = cleaning1 | cleaning2;
    cleaning = ~cleaning;
    if sum(cleaning) < length(remapped{1})
        for j=1:4
            remapped{j}=remapped{j}(cleaning);
        end
        for j=1:numFieldsGenotypes
            chrGenotypes{j}=chrGenotypes{j}(cleaning);
        end
    end
    numImpSNPs = length(remapped{1});
    remapped{2} = str2double(remapped{2});
    remapped{3} = str2double(remapped{3});
    
    
end
remapped = vertcat(remapped{1:end});
remMerged = cell(1,4);
parfor i=1:4
	remMerged{i} = vertcat(remapped{:,i});
end
remapped=remMerged;
clear  files;

numImpSNPs = length(remapped{1});
remMerged = zeros(numImpSNPs,4);
chrCopy = remapped{1};
parfor i=1:numImpSNPs
    remMerged(i,1) = str2double(chrCopy{i}(4:end));
end
parfor i=2:4
    remMerged(:,i) = str2double(remapped{i});
end
remapped=remMerged;
clear remMerged chrCopy
%idx = find(data{1}>9,1,'first');
%remapped{1} = cell2mat(remapped{1}(4:end));
%remapped{1} = cell2mat(remapped{1});

% num = zeros(numImpSNPs,2);
% chrsTemp= remapped{1};
% parfor i=1:size(remapped{1},1)
%     num(i) = str2double(chrsTemp{i}); %use str2num if sorting needed
% end
% remapped{1} = num;
% clear num chrsTemp
%parfor(i=2:4)
%    remapped{i} = cell2mat(remapped{i});
%end
%remapped = sortrows(remapped,[1 2]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
genotypeFiles = dir([path '\*.csv']);
filenamesOnly=cell(length(genotypeFiles),1);
for i=1:length(genotypeFiles)
   filenamesOnly{i}=genotypeFiles(i).name;
end
numGenotypeFiles = length(filenamesOnly);
chrGenotypes=cell(numGenotypeFiles,1);
fID=fopen([path '\' filenamesOnly{1}]);
fieldNamesGenotypes = textscan(fID,repmat('%s',1,23),1,'Delimiter',',');
fclose(fID);
for i=1:numGenotypeFiles
    fID=fopen([path '\' filenamesOnly{i}]);
    chrGenotypes{i} = textscan(fID,'%u%u32%*u%*s%*s%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8%c%u8','Headerlines',1,'Delimiter',',');
    fclose(fID);
end
clear filenamesOnly genotypeFiles
chrGenotypes = vertcat(chrGenotypes{1:end});
remMerged=cell(1,numGenotypeFiles);
parfor i=1:numGenotypeFiles
	remMerged{i} = vertcat(chrGenotypes{:,i});
end
chrGenotypes=remMerged;
clear remMerged

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fID = fopen([mapPath '\allStrains.map']);
map = textscan(fID,'%u8%*s%*u%u32%*[^\n]','Delimiter','\t');
fclose(fID);
numMapSNPs = length(map{1});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%i = current SNP in map file being worked on
%k = current Chr in imputed set being worked with
%j = current SNP in imputed set being worked with


[chrGen_hitsGenA,chrGen_idxGenoB]=ismember([chrGenotypes{1} chrGenotypes{2}],[remapped{1} remapped{2}],'rows');
assemblyConvertFails = remapped{4}~='+';
sprintf('%u SNPs failed assembly converison.\n',sum(assemblyConvertFails))
sprintf('%u of %u imputed SNPs failed assembly converison.\n',sum(~chrGen_hitsGenA), length(chrGenotypes{1}))
for i=1:length(chrGenotypes)
    chrGenotypes{i}(~chrGen_hitsGenA)=[];
    chrGenotypes{i}(assemblyConvertFails)=[];
end

chrGen_idxGenoB=chrGen_idxGenoB(chrGen_idxGenoB~=0);
chrGen_idxGenoB=uint32(chrGen_idxGenoB);

for i=length(remapped)
    remapped{i}=remapped{i}(chrGen_idxGenoB);
    remapped{i}(assemblyConvertFails)=[];
end
   
clear assemblyConvertFails;
remapped(4)=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rem_hitsMap,rem_idxMap]=ismember([remapped{1} remapped{3}], [map{1} map{2}],'rows');
sprintf('%u of %u assembly converted imputed SNPs matched with strains map file.\n',sum(rem_hitsMap), length(chrGenotypes{1}))
for i=1:length(remapped)
    remapped{i}(~rem_hitsMap)=[];
end
for i=1:length(chrGenotypes)
   chrGenotypes{i}(~rem_hitsMap)=[];
end

rem_idxMapNoZeros=rem_idxMap(rem_idxMap~=0);
rem_idxMapNoZeros=uint32(rem_idxMapNoZeros);
rem_idxMap=uint32(rem_idxMap);
%if length(rem_idxMap~=numMapSNPs)
    
map{1}(rem_idxMap)
% if length(map{1}(rem_idxMapNoZeros)) ~= numMapSNPs
%     sprintf('%u of %u map file SNPs failed to convert.\n',numMapSNPs-length(map{1}(rem_idxMapNoZeros)),numMapSNPs);
%     %map = map{1}(~rem_idxMap);
%     
%     for i=1:
%         if 
%     sorted_rem_idxMap=sort(rem_idxMap);
%     counter = 0;
%     sorted_rem_idxMap(sorted_rem_idxMap==
%     idx = find(sorted_rem_idxMap~=
%     for i=1:length(sorted_rem_idxMap)
%      
%     numGenomes = (length(chrGenotypes)-5)/2;    
%     for n=1:numGenomes
%         genomes(n)(chrGenotypes{5+n*2}
%     end   
%     
    
    
numImpSNPs=length(remapped{1});
% for i=length(map)
%     map{i}=map{i}(rem_idxMap);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% misses = hits==0;
% sprintf('%u SNPs were not matched.\n',sum(misses));
numGenomes = (length(chrGenotypes)-5)/2; 
for n=1:numGenomes
    genomes{n}=chrGenotypes{5+n*2};
end      
fID = fopen([path '\log.txt']);
for i=1:numMapSNPs
    if chrGenotypes{6+n*2}(i) < 2 || rem_hitsMap(i)==0
        fprintf(fID,'Imputed Col %u\tRow %u\tChr %u\tPos %u32\tConf\t%u\nMap file row: %u\n',5+n*2,hits(i),chrGenotypes{1}(rem_idxMap(i)),chrGenotypes{2}(rem_idxMap(i)),chrGenotypes{6+n*2}(hits(i)),i);
    end
end
sprintf('%u of %u SNPs matched from a set of %u32.\n',length(hits)-misses,numMapSNPs,numImpSNPs)
strains=cell(numGenomes,1);
for i=1:numGenomes
    strains{i} = fieldNamesGenotypes{5+i*2};
end
parfor i=1:numGenomes
    missing=genomes{i}=='N';
    genomes{i}(missing)='0';
    lenGenome=length(genomes{i});
    fID=fopen(['C:\FastLMM\Plink\' strains{i} '.test.log'],'w');
    fprintf(fID,'SNPs: %u , Missing: %u\n\tReference of %u with %u found during assembly conversion.',lenGenome,missing,numMapSNPs,sum(hits))
    sprintf('%s: SNPs: %u , Missing: %u',strains{i},lenGenome,missing)
    fclose(fID);
    
	fID = fopen([strains{i} '.ped']);
	fprintf(fID,'%s %s 0 0 0 -9 ',strains{i}, strains{i});
	for j=1:length(genomes{i})
		fprintf(fID,'%c %c ',genomes{i}(j),genomes{i}(j));
	end
end
end	
% for 1:length(idxB)
%     if remapped{4}
%     end
% end
% hitsA=hitsA(remapped{4}(idxB)~='+')=0;
% 
% 
% for i=1:numMapSNPs
%     k=1; 
%     while(remapped{1}(k)==map{1}(i) && hits(i)==0)
%         for j=1:numImpSNPs
%             
%             if remapped{3}(j) == map{2}(i) && remapped{4}(j) == '+'
%                 m=1;
%                 while(chrGenotypes{1}(m)==map{1}(i))
%                     if chrGenotypes{2}(m) == remapped{2}(j)
%                          hits(i)=m;
%                          j=numImpSNPs+1;
%                          break
%                     end
%                     m=m+1;
%                 end
%             end
%         end
%         k=k+1;
%     end
% end
% %  [LIA,LOCB] = ismember(A,B,'rows') returns a vector containing true where the rows of A are
% %     also rows of B and false otherwise. also returns an array LOCB containing the
% %     lowest absolute index in B for each element in A which is a member of
%     B and 0 if there is no such index.
% 

%               
%     tableLength = length(chrGenotypes{i}{1});
%     while j<=tableLength && n<=numMapSNPs
%         if chrGenotypes{i}{2}(j)==remapped{2}(k)
%             chrom = str2double(remapped{1}{k});
%             if (chrGenotypes{i}{1}(j) < 9 && chrGenotypes{i}{1}(j)==chrom) ||  (chrGenotypes{i}{1}(j) >9 && chrGenotypes{i}{1}(j)==chrom)
%                 if remapped{4}(k) == '+'
%                     found = n;
%                     for z=1:numMapSNPs
%                         if map{1}(z) == remapped{1}(k)
%                             break
%                         end
%                     end
%                     while map{1}(z) == chrom && z<=numMapSNPs
%                         if  map{2}(z) == remapped{3}(k)                             
%                             for m=0:numGenomes
%                                 genomes{m}(n)=chrGenotypes{i}{5+m*2}(j);
%                                 hits(z)=k;
%                                 if chrGenotypes{i}{6+m*2}(j) < 2
%                                     fprintf(fID,'Col %u\tRow %u\tChr %u\tPos %u32\tConf\t%u\n',5+m*2,j,chrGenotypes{i}{1}(j),chrGenotypes{i}{2}(j),chrGenotypes{i}{6+m*2}(j));
%                                 end
%                             end
%                             n=n+1;
%                             break
%                         end
%                         z=z+1;
%                     end
%                 end
%                 if found == n
%                     fprintf(fID,'table %u\trow%u\tremappedChr %u\tremappedRow %u\tremappedBP %u\tzMAPindex %u\tkRemappedMergedIndex %u\n',i,j,remappedChrs(n),n,remappedMerged{2}(n),z,k); 
%                 end
%             end
%             k=k+1;
%         end
%         j=j+1;
%     end
%     j=1;
% end
% fclose(fID);
    
	
