function Dmat = cell2num(C)
%
% This only works if each cell is the same length, and is a vector of
% numbers.
%

Dmat = cellfun(@(x) reshape(x,1,1,length(x)),C,'Uni',0);
Dmat = cell2mat(Dmat);
Dmat = squeeze(Dmat);


end