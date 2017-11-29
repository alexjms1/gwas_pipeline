function mergeResults(root, filenames, testQvals, cov, append)
% This function takes results by individual chromosomes for specified
% phenotypes (filenames) and merges them into a single tab-delimited file
% named 'C:\FastLMM\CPP_MKL\Output\' root '\Merged\' root '.' filenames{phenotype}
% '.all.fastlmm.txt'. Resultant files can be loaded into Haploview for
% visualization.  Additionally, since Q-values in individual
% chromosomes' results files are not able to control for false-discovery
% rate (FDR) over the entire genome but rather only control within a single
% chromosome (since they do not have access to the entire set of data when
% created), we compute them using the entire genome's results here.
% Summary results containing the most significant results per phenotype are
% generated ('.LowestPsQs' / '.LowestPs') for quick assessment of results. 
% 
% INPUTS:
%   root: main directory where FaST-LMM results are found, code as
%       'C:\FastLMM\CPP_MKL\Output\' root '\' root '.' filenames.
%   filenames: filenames as above, specifying phenotypes, in the form  root
%       '.'  filename[s] '.Chrom' chr '.fastlmm.txt']
%   testQvals: a variable indicating whether to test Q values across all
%       chromosomes (testQvals > 0), and if so, how extensive those results
%       should be.  testQvals = 1: Q-values, Story FDR, BH FDR.  testQvals = 2:
%       additionally reports Story FDR polynomial, Story Q polynomial.
%       testQvals = 3 additionally reports BH FDR polynomial.
%   cov: a boolean specifying whether covariate(s) were used, needed to 
%       create output files with the correct columns
%   append: a boolean specifying whether, in generating '.LowestPs' and
%       '.LowestPsQs' reports, whether results should be appended to existing
%       files or overwrite them.
if nargin ~= 5 % note to self: should have used repmat for changes to Q reports
    error('Not enough arguments passed to function. 5 are required');
end
if testQvals > 3 || testQvals < 0
    error('Q value test integer not within range of 0 - 3');
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

for phenotype=1:numFiles
    chrom=1;
    currFile = ['C:\FastLMM\CPP_MKL\Output\' root '\' root '.'  filenames{phenotype} '.Chrom' num2str(chrom) '.fastlmm.txt'];
    if(exist(currFile,'file'))
        fileID = fopen(currFile);
        if ~cov
            columnFormat = repmat('%s',1,25);
        else
            columnFormat = repmat('%s',1,27);
        end
        dataHeaders = textscan(fileID,columnFormat,1,'Delimiter','\t');
        fclose(fileID);
        dataHeaders{1}{1} = 'CHR';
        dataHeaders{2}{1} = 'SNP';
    else
        disp(['File ' root '.'  filenames{phenotype} ' for chromosome ' num2str(chrom) ' not found']);
        continue;
        %return
    end
    data = cell(20,1);
    if ~cov
        columnFormat = '%s %u8 %u8 %u %s %.10f %.10f %u16 %u16 %f %f %f %f %s %f %f %f %f %f %f %f %f %f %u64 %u64';
    else
        columnFormat = '%s %u8 %u8 %u %s %.10f %.10f %u16 %u16 %f %f %f %f %s %f %f %f %f %f %f %f %f %f %f %f %u64 %u64';
    end
    while exist(currFile,'file')
        fileID=fopen(currFile);
        data{chrom} = textscan(fileID,columnFormat,'HeaderLines',1,'Delimiter','\t');
        fclose(fileID);
        chrom=chrom+1;
        currFile = ['C:\FastLMM\CPP_MKL\Output\' root '\' root '.' filenames{phenotype} '.Chrom' num2str(chrom) '.fastlmm.txt']; 
    end
    if ~(chrom==21)
        disp(['File ' root '.'  filenames{phenotype} ' for chromosome ' num2str(chrom) ' not found']);
        continue;
    end
    pvals = data{1}{6};
    pMin=1;
    [pMinSrc, pMinSrcIndx] = min(data{1}{6}); %NEW CODE
    if pMinSrc < pMin
        pMin = pMinSrc;
        pIdx = pMinSrcIndx;
        pChr = 1;
    end
    for i=2:20
        [pMinSrc, pMinSrcIndx] = min(data{i}{6});
        if pMinSrc < pMin
            pMin = pMinSrc;
            pIdx = pMinSrcIndx;
            pChr = i;
        end
        newLen = length(data{i}{6});
        pvals(end+1:end+newLen)=data{i}{6};
    end
    lowestP_chr(phenotype) = pChr;
    lowestP_pos(phenotype) = data{pChr}{4}(pIdx);
    lowestP_rs{phenotype} = data{pChr}{1}{pIdx};
    pvals=pvals';
    lowestPs(phenotype) = min(pvals);
    switch testQvals
        case 1
            additionalCols= cell(3,1);
        case 2
            additionalCols=cell(5,1);
        case 3
            additionalCols=cell(6,1);
    end
          
    if testQvals
       additionalCols{1} = 'all_St_pFDR';
       additionalCols{2} = 'all_St_Q';
       additionalCols{3} = 'all_BH_FDR';
       [fdr, q] = mafdr(pvals);
       lowestStQs(phenotype)=min(q);
       lowestStFDRs(phenotype)=min(fdr);
       fdrBH = mafdr(pvals,'BHFDR',true);
        lowestBHFDRs(phenotype)=min(fdrBH);  
        if testQvals > 1
            additionalCols{4} = 'all_St_FDR_Poly';
            additionalCols{5} = 'all_St_Q_Poly';
            [fdrStPoly, qStPoly] = mafdr(pvals,'Method','Polynomial'); 
        end 
        if testQvals > 2 
            additionalCols{6} = 'all_BH_FDR_Poly';
            fdrBHPoly = mafdr(pvals,'BHFDR',true,'Method','Polynomial');
        end
    end
    fileID=fopen(['C:\FastLMM\CPP_MKL\Output\' root '\Merged\' root '.' filenames{phenotype} '.all.fastlmm.txt'],'w');
    numHeaders = size(dataHeaders,2)-1;
    for i=1:numHeaders
        fprintf(fileID,'%s\t',dataHeaders{i}{1});
        if i==7  && testQvals
            if testQvals == 1
                fprintf(fileID,'%s\t%s\t%s\t%s\t',additionalCols{1},additionalCols{2},additionalCols{3});
            elseif testQvals == 2
                fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t',additionalCols{1},additionalCols{2},additionalCols{3},additionalCols{4},additionalCols{5});
            else
                fprintf(fileID,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t',additionalCols{1},additionalCols{2},additionalCols{3},additionalCols{4},additionalCols{5},additionalCols{6});
            end 
       end              
    end
    fprintf(fileID,'%s\n',dataHeaders{i+1}{1});
    if ~testQvals 
      if ~cov
        columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
      else
        columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';  
      end
      dataSize = size(data,1);
        for i=1:dataSize
            perChromSize = size(data{i}{1},1);
            for j=1:perChromSize
                if ~cov
                    fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
                else
                    fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j), data{i}{26}(j), data{i}{27}(j));
                end
            end
        end
    else
          count=0;
        dataSize = size(data,1);
        if testQvals==1
            if ~cov
                columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
            else
                columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
            end
            for i=1:dataSize
                perChromSize=size(data{i}{1},1);
                for j=1:perChromSize
                        count = count+1;
                        if ~cov
                            fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
                        else
                            fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j), data{i}{26}(j), data{i}{27}(j));
                        end
                end
            end
        elseif testQvals==2
            if ~cov
                columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
            else
                columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
            end
           for i=1:dataSize
                perChromSize=size(data{i}{1},1);
                for j=1:perChromSize
                    count = count+1;
                    if ~cov
                        fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), fdrStPoly(count), qStPoly(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
                    else
                        fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), fdrStPoly(count), qStPoly(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j),  data{i}{26}(j),  data{i}{27}(j));
                    end
                end
            end

        else
            if ~cov
                columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
            else
                columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
            end
           for i=1:dataSize
                perChromSize=size(data{i}{1},1);
                for j=1:perChromSize
                    count = count+1;
                    if ~cov
                        fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), fdrStPoly(count), qStPoly(count), fdrBHPoly(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
                    else
                        fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), fdrStPoly(count), qStPoly(count), fdrBHPoly(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j), data{i}{26}(j), data{i}{27}(j));
                    end
                end
           end
        end

    end
    fclose(fileID);    
end
if append
    writePriv = 'a';
else
    writePriv = 'w';
end
if numFiles
    if testQvals
        summaryFile = fopen(['C:\FastLMM\CPP_MKL\Output\' root '\Merged\for.' root '.LowestPsQs.txt'],writePriv);
        fprintf(summaryFile,'Phen\tPval\tpFDR\tQval\tBHFDR\tbestChr\tPos\tbestRS\n');
        for cnt=1:numFiles
            fprintf(summaryFile,'%s\t%.10f\t%.15f\t%.15f\t%.15f\t%.0u\t%u\t%s\n',filenames{cnt}, lowestPs(cnt), lowestStFDRs(cnt), lowestStQs(cnt), lowestBHFDRs(cnt), lowestP_chr(cnt), lowestP_pos(cnt), lowestP_rs{cnt});
        end
    else
        summaryFile = fopen(['C:\FastLMM\CPP_MKL\Output\' root '\Merged\for.' root '.LowestPs.txt'],writePriv);
        fprintf(summaryFile,'Phen\tPval\tbestChr\tPos\tbestRS\n');
        for cnt=1:numFiles
            fprintf(summaryFile,'%s\t%.10f\t%.0u\t%u\t%s\n',filenames{cnt}, lowestPs(cnt), lowestP_chr(cnt), lowestP_pos(cnt), lowestP_rs{cnt});
        end
    end 
end
fclose(summaryFile);
format;
end