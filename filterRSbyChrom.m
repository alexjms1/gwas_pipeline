function filterRSbyChrom(genome)
%Input is a plink binary genome file, function returns a series of tab delimited
%file - for each chromosome, 1 file with SNPs from that chromosome only & 1
%file for SNPs found on all other chromosomes (to prevent proximal
%contamination, and calculated relatedness random effect matrix using SNPs
%not in LD with the SNPs you're presently testing.  Mouse chromosome 23 is
%treated as chromosome 20 (skips 20-22).
%The function will create the directory 'C:\FastLMM\Cpp_MKL\Files\SNPs\' genome, 
%containing 2 sets of files, 'filtered.snp_ids.chr.' chr '.included.txt' containing SNP IDs from only a single
%chromosome, and 'filtered.snp_ids.chr.' chr '.excluded.txt', containing SNPS from all other chromosomes, for
%a total of 40 files (20 x 2).  Also creates a file genome
%'.ChromSNP.Count.txt', a count of SNPs per chromosome.
%INPUTS:
%   genome: binary plink genome file, found in hard-coded directory
%       'C:\FastLMM\Cpp_MKL\'
if nargin ~= 1
    error('Not enough arguments passed to function. 1 is required');
end

if ~exist(['C:\FastLMM\Cpp_MKL\' genome '.bim'],'file')
    display(['Input genome ' genome ' not found.'])
    return
end
if ~exist(['C:\FastLMM\Cpp_MKL\Files\SNPs\' genome ],'dir')
    mkdir(['C:\FastLMM\Cpp_MKL\Files\SNPs\' genome ]);
end
fileID = fopen(['C:\FastLMM\Cpp_MKL\' genome '.bim']);
SNPs=textscan(fileID,'%.0u %s %*[^\n]','delimiter','\t');
fclose(fileID);
filteredSNPs = cell(20,1);
for chrom = 1:19
    filteredSNPs{chrom}=SNPs{2}(SNPs{1}==chrom);
end
filteredSNPs{20}=SNPs{2}(SNPs{1}==23);

filename = ['C:\FastLMM\Cpp_MKL\Files\SNPs\' genome '\' genome '.ChromSNP.Count.txt'];
fileID = fopen(filename,'w');
fprintf(fileID,'Chromosome:\t SNP count\n');

for i=1:19
    fprintf(fileID,'Chromosome %u:\t%u\n',i,size(filteredSNPs{i},1));
end
fprintf(fileID,'Chromosome %u:\t%u\n',23,size(filteredSNPs{20},1));
fclose(fileID);
clear filename

parfor chrom=1:20
    filename = ['C:\FastLMM\Cpp_MKL\Files\SNPs\' genome '\' genome '.filtered.snp_ids.chr.' num2str(chrom) '.included.txt'];
    fileID=fopen(filename,'w');
    fprintf(fileID, '%s\n', filteredSNPs{chrom}{:});
    fclose(fileID);
end
for chrom=1:20
    filename = ['C:\FastLMM\Cpp_MKL\Files\SNPs\' genome '\' genome '.filtered.snp_ids.chr.' num2str(chrom) '.excluded.txt'];
    fileID=fopen(filename,'w');
    for i=1:20
       if chrom~=i          
           fprintf(fileID, '%s\n', filteredSNPs{i}{:});
       end
    end
    fclose(fileID);
end  
end
    
    
    