function result = duplicates(array)
%Returns an array of duplicate SNPs by chromosome from an input of a cell
%array where array{1} is chromosome number and array{2} is SNP IDs
if iscell(array) 
    numChr=length(unique(array{1}));
    result=cell(numChr,1);
    for i=1:numChr
          [b, ~, c] = unique(array{2}(array{1}==i));
          result{i}=b(accumarray(c(:),1)>1);
          clear b c
    end        
else    
    [b, ~, c] = unique(array);
    result=b(accumarray(c(:),1)>1);
end

end
