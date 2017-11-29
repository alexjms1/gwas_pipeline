function [phenoNames, phenoNumbers] = fastLMM(genome, root, phenotypeFile, phenotypesTesting, return_only_names, varargin)
% This function calls FaST-LMM on specified phenotypes using a specified
% genome, performing leave-one-chromosome-out linkage disequilibrium
% control
% INPUTS:
%   genome: a plink binary genome file, corresponding to your sample, in the directory 'C:\FastLMM\Cpp_MKL'
%   root: a name for the root directory of your result output,
%       corresponding to 'Output\' root.  Filenames outputted will be of the 
%       format 'Output\' root '\' root '.' phenotype '.' chromosome '.fastlmm.txt'. 
%   phenotypeFile: name of the file containing your strain/phenotype data,
%       found in 'Files\' phenotypeFile '.txt' 
%   phenotypesTesting: either a cell array of phenotype names, or an a
%       vector of coumn numbers, corresponding to the phenotypes you want
%       analyzed; if numbers, column 3 in your data file should be considered
%       the first phenotype
%   return_only_names: a boolean, if true, FaST-LMM is not run and the
%       function returns phenoNames and phenoNumbers; uusually for diagnostic purposes 
%   vargin: for covariate use, specifiy filename; if only 1 file, it is used for all analyses,
%       if more than one, they are matched to successive phenotypes, until number
%       of covariate files is expended, at which point it will use the last-most
%       file
% OUTPUTS:
%   phenoNames: column names corresponding to the phenotypes you selected
%       for analysis
%   phenoNumbers: phenotype column numbers corresponding to the phenotypes you
%   selected for analyses (column num - 2)
if nargin < 5
    error('Not enough arguments passed to function. 5 are required');
end
%varargin = 
phenos = importdata(['C:\FastLMM\CPP_MKL\Files\' phenotypeFile '.txt'],'\t');
if iscell(phenotypesTesting)
    phenoNumbers = zeros(size(phenotypesTesting,2));
    for i=1:size(phenotypesTesting,2)
        for j=3:size(phenos.textdata,2)
            if phenotypesTesting{i} == phenos.textdata{1,j}
                phenoNumbers(i) = j-2;
                j=size(phenos.textdata,2)+1;
            end
        end
        if phenoNumbers(i) == 0
            disp(['!!Invalid phenotype ' phenotypesTesting{i} '!!']);
            return;
        end
    end
    phenoNames = phenotypesTesting;
else
    phenoNames = cell(size(phenotypesTesting,2),1);
    for i=1:size(phenotypesTesting,2)
        if phenotypesTesting(i)+2 <= size(phenos.textdata,2)
            phenoNames{i} = phenos.textdata{1,phenotypesTesting(i)+2};
        else
            disp(['!!Invalid phenotype ' phenotypesTesting(i) '!!']);
            return;
        end
    end
    phenoNumbers = phenotypesTesting;
end
for i=1:size(phenotypesTesting,2)
    phenoNames{i} = [phenoNames{i} '_' num2str(phenoNumbers(i))];
end
if return_only_names
    return
end
if ~exist(['C:\FastLMM\CPP_MKL\Output\' root],'dir')
    mkdir(['C:\FastLMM\CPP_MKL\Output\' root]);
end
if ~exist(['C:\FastLMM\CPP_MKL\Output\' root '\Merged'],'dir')
    mkdir(['C:\FastLMM\CPP_MKL\Output\' root '\Merged']);
end
for i=1:size(phenotypesTesting,2)
    parfor chromosome=1:20
        command = ['cd C:\FastLMM\Cpp_MKL & fastlmmc -bfile ' genome ' -bfileSim ' genome ' -pheno Files\' phenotypeFile '.txt -mpheno '  num2str(phenoNumbers(i)) ' -setOutputPrecision 5 -verboseout -out Output\' root '\' root '.' phenoNames{i} '.Chrom' num2str(chromosome) '.fastlmm.txt'];
        if(size(varargin,2)>0)
            if size(varargin,2>1)
                if size(varargin,2)<=i
                    command = [command ' -covar Files\' varargin{i} '.txt']; 
                else
                    display('Warning: number of covariate files is greater than the number of phenotypes');
                end
            else
                command = [command ' -covar  Files\' varargin{1} '.txt']; 
            end
        end
        command = [command ' -v -extract Files\SNPs\' genome '\' genome '.filtered.snp_ids.chr.' num2str(chromosome) '.included.txt -extractSim Files\SNPs\' genome '\' genome '.filtered.snp_ids.chr.' num2str(chromosome) '.excluded.txt -maxThreads 1']; 
        dos(command,'-echo');     
    end
end

end