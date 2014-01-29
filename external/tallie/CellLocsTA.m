function [x,y] = CellLocsTA(map)
%
% Allows you to choose the locations of cells on a chosen slice map (that
% was previously saved as a .mat) and kicks out a .tif and a .mat of all
% that info, indexed by the time and date.
%

load(map); % e.g. load('ACdv.mat');
figure, imshow(ACdv)
[x,y] = ginput;
ans1 = input('Please give the cell indices as a matrix: ');
hold on, scatter(x,y,'filled')
for k = 1:length(x), text(x(k)+3,y(k)+3,num2str(ans1(k))), end
saveas(gcf,['CellLocs_' datestr(clock,'yymmdd_HHMMSS')],'tif')
s = strsplit(map,'.');
save(['CellLocs_' datestr(clock,'yymmdd_HHMMSS')],s{1},'x','y')

end

