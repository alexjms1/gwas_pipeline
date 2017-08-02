function generateFullSummary(testQvals,regerenatePsQs,varargin)
% This function (re-)generates best result summary files phenotypes from
% output ...LowestPs/PsQs.txt files.  It can also regenerate the Ps/Qs from
% output files.  Summary files will be named 'Output\Summary.' mmmdd_yyyy
% '.PsQs.txt'.
%INPUTS:
%   testQvals: a boolean, specifying whether summaries should contain Q
%       values
%   regeneratePsQs: a boolean, specifying whether the Ps/Qs should be
%       regenerated when creating the summary.  This will call
%       'generateLowestPsQs' with this function's testQvals argument = 2
%   varargin: can contain a cell array of roots (directories) in which merged
%       FaST-LMM files are found, within their merged directory; if none is
%       specified, uses all directories found in '\Output'.  This will call
%       
if nargin == 2
    roots = dir('Output\');
    roots = roots([roots.isdir]);
elseif nargin == 3
    roots = struct([]);
    for i=1:length(varargin)
        roots(i).name = varargin{i};
    end
elseif nargin > 3
    error('Passed too many arguments - put root names in a single cell array');
else
    error('Not enough arguments passed; fuction requests testQval (0 or 1), regenerate Ps / Qs (0 or 1); and optionally a cell array of output roots / directories ');
end
for i=1:length(roots)
    if strcmpi(roots(i).name,'Prev')  
        roots(i)=[];
        break
    end
end
for i=1:length(roots)
    if strcmp(roots(i).name, '.')
        roots(i)=[];
        break;
    end
end
for i=1:length(roots)
    if strcmp(roots(i).name, '..')
        roots(i)=[];
        break;
    end
end

numRoots=length(roots);
data=cell(numRoots,1);
if regerenatePsQs
    for i=1:numRoots
       filename = dir(['Output\' roots(i).name '\merged\*.txt']);
       for j=1:length(filename)
           if strfind(filename(j).name,'LowestPs')
               filename(j) = [];
               j=j-1;
           end
       end
       filenamesOnly=cell(length(filename),1);
       for j=1:length(filename)
           filenamesOnly{j}=filename(j).name;
       end
       filenamesOnly=strrep(filenamesOnly,'.all.fastlmm.txt','');
       filenamesOnly=strrep(filenamesOnly,[roots(i).name '.'],'');
       generateLowestPsQs(roots(i).name,filenamesOnly,testQvals,2)  
    end
end
for i=1:numRoots
   filename = dir(['Output\' roots(i).name '\Merged\*LowestPs*.txt']);
   if length(filename)>1
       disp(['More than one summary file found for the root: ' roots(i).name '. Proceding with first: ' filename(1).name '.'])
   elseif isempty(filename)
       disp(['No summary file found for the root: ' roots(i).name '.'])
   end
   if testQvals
       columnFormat = '%s %.10f %.15f %.15f %.15f %.0u %f %s';
   else
       columnFormat = '%s %.10f %.0u %f %s';
   end
   fileID = fopen(['Output\' roots(i).name '\Merged\' filename(1).name]);
   if fileID<1
        display(['Error opening ' filename(1).name]'.');
   else
   	    data{i} = textscan(fileID,columnFormat,'Delimiter','\t');
   end
   fclose(fileID);
end
fileID=fopen(['Output\Summary.' datestr(now,'mmmdd_yyyy') '.PsQs.txt'],'w');
if testQvals
    fprintf(fileID,'Phen\tPval\tpFDR\tQval\tBHFDR\tbestChrPos\t\tbestRS\n');
    columnFormat = '%s\t%.10f\t%.15f\t.15f\t.15f\t%.0u\t%f\t%s\n';
else
    fprintf(fileID,'Phen\tPval\tbestChrPos\t\tbestRS\n');
    columnFormat = '%s\t%.10f\t%.0u\t%f\t%s\n';
end
for i=1:numRoots
    numPhenos=length(data{i}{1});
    if testQvals
        for j=1:numPhenos
            fprintf(fileID,columnFormat,[roots(i).name  ': ' data{i}{1}{j}],data{i}{2}(j),data{i}{3}(j),data{i}{4}(j),data{i}{5}(j),data{i}{6}(j), data{i}{7}(j), data{i}{8}{j});
        end
    else
        for j=1:numPhenos
            fprintf(fileID,columnFormat,[roots(i).name  ': ' data{i}{1}{j}],data{i}{2}(j),data{i}{3}(j),data{i}{4}(j),data{i}{5}(j),data{i}{6}{j});
        end
    end
end
fclose(fileID);
end

