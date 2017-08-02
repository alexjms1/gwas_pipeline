% helper function to convert xls files to text tab delimited.
function xlsToTxt(path)
format LONGG;
if path(end) == '\'
    path = path(1:end-1);
end

remapFiles = dir([path '\remapped_merged*.txt.xls']);
filenamesOnly=cell(length(remapFiles),1);
parfor j=1:length(remapFiles)
   filenamesOnly{j}=remapFiles(j).name;
end
filenamesOnly = sortn(filenamesOnly);
numRemapped = length(filenamesOnly);
raw=cell(numRemapped,1);

parfor j=1:numRemapped
    [~,~,raw{j}]=xlsread([path '\' filenamesOnly{j}]);
end
for j=1:numRemapped
    save([path '\' filenamesOnly{j}(1:end-4)],raw{j},'-ASCII','-double','-tabs');
end