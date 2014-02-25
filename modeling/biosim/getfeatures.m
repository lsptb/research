function [expfeatures,simfeatures] = getfeatures(datafiles,simfiles)
% Purpose: extract features sets from a list of experimental and simulated
% data files.

expfeatures = {};
simfeatures = {};
files = cat(2,datafiles,simfiles);
for i=1:length(files)
  file = files{i};
  [fpath,fname] = fileparts(file);
  setfile = fullfile(fpath,[fname '_features.mat']);
  if exist(setfile,'file')
    load(setfile,'features');
  else
    [features,parms] = CharacterizeCells(file);
    save(setfile,'features','parms','file');
  end
  if ismember(file,datafiles), expfeatures{end+1} = features; end
  if ismember(file,simfiles), simfeatures{end+1} = features; end
end

%{
PSEUDOCODE:
% analyze single data sets
foreach file in (datafiles, simfiles)
     setfile = [file '_features'];
     if exist(setfile)
          load(setfile,'features');
     else:
          [features,parms] = characterize(file)
          save(setfile, 'features','parms','file');
     if file in datafiles: expfeatures{file} = features % one feature set
     if file in simfiles: simfeatures{file} = features % N feature sets
%}