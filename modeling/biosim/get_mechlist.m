function [list,files]=get_mechlist(DBPATH)
if nargin<1
  DBPATH = '/space/mdeh3/9/halgdev/projects/jsherfey/code/modeler/database';
end
if ~exist(DBPATH,'dir')
  DBPATH = 'C:\Users\jsherfey\Desktop\My World\Code\modelers\database';
end
d=dir(DBPATH);
list = {d(cellfun(@(x)any(regexp(x,'.txt$')),{d.name})).name};
files = cellfun(@(x)fullfile(DBPATH,x),list,'unif',0);