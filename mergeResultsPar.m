function mergeResultsPar(root, filenames, testQvals, append)
% A parallel processing version of mergeResults - see this function for a
% complete description.  Currently unused / depreciated (no performance
% benefits?).
if nargin ~= 4 %shoulda used repmat for changes to Q reports
    error('Not enough arguments passed to function. 4 are required');
end
if testQvals > 3
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
%   if testQvals > 1
%        lowestStQsPoly = zeros(1,numFiles);
%        lowestStFDRPoly = zeros(1,numFiles);
%    end 
%    if testQvals > 2 
%        lowestStFDRsPoly = zeros(1,numFiles);
%        lowestBHFDRsPoly = zeros(1,numFiles);
%    end




lowestP_chr = zeros(1,numFiles);
lowestP_pos = zeros(1,numFiles);
lowestP_rs = cell(1,numFiles);

parfor phenotype=1:numFiles
    chrom=1;
    currFile = ['C:\FastLMM\CPP_MKL\Output\' root '\' root '.'  filenames{phenotype} '.Chrom' num2str(chrom) '.fastlmm.txt'];
    if(exist(currFile,'file'))
        fileID = fopen(currFile);
        columnFormat = repmat('%s',1,25);
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
    columnFormat = '%s %u8 %u8 %u %s %.10f %.10f %u16 %u16 %f %f %f %f %s %f %f %f %f %f %f %f %f %f %u64 %u64';
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
        %return
    end
    pvals = data{1}{6};
    pMin=1;
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
    fileID=fopen(['C:\FastLMM\CPP_MKL\Output\' root '\Merged\parfor\' root '.' filenames{phenotype} '.all.fastlmm.txt'],'w');
    numHeaders = size(dataHeaders,2)-1;
    for i=1:numHeaders
        fprintf(fileID,'%s\t',dataHeaders{i}{1});
        if i==7  && testQvals %shoulda used repmat
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
        columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
        dataSize = size(data,1);
        for i=1:dataSize
            perChromSize = size(data{i}{1},1);
            for j=1:perChromSize
                fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
            end
        end
    else
        count=0;
        dataSize = size(data,1);
        if testQvals==1
            columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';            
            for i=1:dataSize
                perChromSize=size(data{i}{1},1);
                for j=1:perChromSize
                    count = count+1;
                    fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
                end
            end
        elseif testQvals==2
            columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
           for i=1:dataSize
                perChromSize=size(data{i}{1},1);
                for j=1:perChromSize
                    count = count+1;
                    fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), fdrStPoly(count), qStPoly(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
                end
            end

        else
            columnFormat = '%.0u\t%s\t%u\t%.0f\t%s\t%.10f\t%.10f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.15f\t%.0u\t%.0u\t%f\t%f\t%f\t%f\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%.0u\t%.0u\n';
           for i=1:dataSize
                perChromSize=size(data{i}{1},1);
                for j=1:perChromSize
                    count = count+1;
                    fprintf(fileID,columnFormat,data{i}{2}(j), data{i}{1}{j},  data{i}{3}(j), data{i}{4}(j), data{i}{5}{j}, data{i}{6}(j), data{i}{7}(j), fdr(count), q(count), fdrBH(count), fdrStPoly(count), qStPoly(count), fdrBHPoly(count), data{i}{8}(j), data{i}{9}(j), data{i}{10}(j), data{i}{11}(j), data{i}{12}(j), data{i}{13}(j), data{i}{14}{j}, data{i}{15}(j), data{i}{16}(j), data{i}{17}(j), data{i}{18}(j), data{i}{19}(j), data{i}{20}(j), data{i}{21}(j), data{i}{22}(j), data{i}{23}(j), data{i}{24}(j), data{i}{25}(j));
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
        summaryFile = fopen(['C:\FastLMM\CPP_MKL\Output\' root '\Merged\parfor.' root '.LowestPsQs.txt'],writePriv);
        fprintf(summaryFile,'Phen\tPval\tpFDR\tQval\tBHFDR\tbestChrPos\t\tbestRS\n');
        for cnt=1:numFiles
            fprintf(summaryFile,'%s\t%.10f\t%.15f\t%.15f\t%.15f\t%.0u\t%u\t%s\n',filenames{cnt}, lowestPs(cnt), lowestStFDRs(cnt), lowestStQs(cnt), lowestBHFDRs(cnt), lowestP_chr(cnt), lowestP_pos(cnt), lowestP_rs{cnt});
        end
    else
        summaryFile = fopen(['C:\FastLMM\CPP_MKL\Output\' root '\Merged\parfor.' root '.LowestPs.txt'],writePriv);
        fprintf(summaryFile,'Phen\tPval\tbestChrPos\t\tbestRS\n');
        for cnt=1:numFiles
            fprintf(summaryFile,'%s\t%.10f\t%.0u\t%u\t%s\n',filenames{cnt}, lowestPs(cnt), lowestP_chr(cnt), lowestP_pos(cnt), lowestP_rs{cnt});
        end
    end 
end
fclose(summaryFile);
format;

end