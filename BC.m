function [imported, lambdas] = BC(filename)
%Box-cox data normalization, ignoring entries that are '-9', the
%marker for missing data in FaST-LMM
%INPUTS: 
%   filename, a tab-delimited file containing data in columns 3 onward
%       (see data file example)
%OUTPUTS: 
%   imported, a box-cox transformed dataset (transformed individually
%       by column
%   lambdas, lambda values for each column used during transformation
%FILES: 
%   filename_BCtrans.txt, a tab-delimited file of box-cox transformed data
%   filename_BClambdas.txt, a tab-delimited file of box-cox lambdas
    
if nargin ~= 1
    error('Not enough arguments passed to function. 1 is required');
end
textExtension = strfind(filename,'.txt');
if textExtension 
    fileWOext = filename(1:end-4);
else
    fileWOext = filename;
    filename = [filename '.txt'];
end

if exist(filename,'file')
    imported = tdfread(filename);
    headers = fieldnames(imported);
    numCols = length(headers);
   numRows = length(imported.(headers{1}));

    %numrows = size(imported.data,1);
    %numcols = size(imported.data,2);
    
   % newData = zeros(numRows,numCols);
    lambdas = zeros(numCols,1);
 %   currData = zeros(numRows,1);   
    for x = 3:numCols
        currData = imported.(headers{x});
        currData=currData((currData~=-9));
        minVal = min(currData);
        if minVal <= 0
            currData = currData + (-1*minVal) + 1;
        end
     %   imported.(headers{x})(find(imported.(headers{x})==-9)) = NaN;
        [currData, lambdas(x)] = boxcox(currData);
        indCnt = 1;
        for y=1:numRows
            if imported.(headers{x})(y)~= -9 
                imported.(headers{x})(y)=currData(indCnt);
                indCnt = indCnt+1;
            end
        end
    end
  %  newData(find(newData==NaN)) = -9;
    tdfwrite([fileWOext '_BCtrans.txt'],imported);
    fID = fopen([fileWOext '_BClambdas.txt'],'w');
    for x=1:numCols-1
        fprintf(fID,'%s\t',headers{x});
    end
    fprintf(fID,'%s\n',headers{end});
    for x=1:numCols-1
        fprintf(fID,'%.15f\t',lambdas(x));
    end
    fprintf(fID,'%.15f\n',lambdas(end));
    fclose(fID);
else
    disp(['file ' filename ' not found.']);
end
end



