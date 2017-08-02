% A quick wrapper script for redoing the calculation of lowest Ps/Qs for a
% set of analyses previously conducted.  Roots contains all the different
% phenotype analyses we want to redo lowets Ps/Qs calculations for.
roots = cell(10,1);
roots{1} = 'bctrans50';
roots{2} = 'bctrans50_no47';
roots{3} = 'bctransAll';
roots{4} = 'bctransNBL';
roots{5} = 'untrans50';
roots{6} = 'untrans50_no47';
roots{7} = 'untransAll';
roots{8} = 'untransAll_47NL';
roots{9} = 'untransNBL';
roots{10} = 'untransNBL_47NL';

for i=1:size(roots,1)
   filenames = dir(['Output\' roots{i} '\merged\*.txt']);
   for j=1:length(filenames)
       if strfind(filenames(j).name,'LowestPs')
           filenames(j) = [];
           break
       end
   end
   filenamesOnly=cell(length(filenames),1);
   for j=1:length(filenames)
       filenamesOnly{j}=filenames(j).name;
   end
   clear filenames
   generateLowestPsQs(roots{i},filenamesOnly,1)  
end