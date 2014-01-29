function [ref_EEG_coords] = ts_extract_hpts_from_fiff(datafile,hptsfile);
%function [ref_EEG_coords] = ts_extract_hpts_from_fiff(datafile,hptsfile);
%
% uses Uutela's meg_pd fiff access toolbox hpipoints function
%
% Last Mod: 08/05/09 by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;

tmp = strfind(datafile,'.fif');
if isempty(tmp)
  fprintf('%s: error: datafile must have .fif extension\n',mfilename);
  return;
end;
stem = datafile(1:tmp(end)-1);
tmp = strfind(datafile,'/');
stem = stem(tmp(end)+1:end);
if ~exist('hptsfile','var')
  hptsfile = sprintf('%s.hpts',stem);
end;

ref_EEG_coords=[];

if ~exist(datafile,'file')
  error('file %s not found',datafile);
end;

[coords,kind,num] = hpipoints(datafile);

points = [];
for p=1:length(num)
  switch kind(p)
    case 1
      points(p).name = 'cardinal';
    case 2
      points(p).name = 'hpi';
    case 3
      points(p).name = 'eeg';
    case 4
      points(p).name = 'extra';
    otherwise
      points(p).name = [];
  end;
  points(p).num = num(p);
  points(p).x = coords(1,p);
  points(p).y = coords(2,p);
  points(p).z = coords(3,p);
end


fid=fopen(hptsfile,'wt');
npoints = length(points);
for p=1:npoints
  switch points(p).name
    case {'cardinal','extra'}
      fprintf(fid,'%s %03d %f %f %f\n',...
        points(p).name,...
        points(p).num,...
        points(p).x,...
        points(p).y,...
        points(p).z);
  end;
end;
fclose(fid);

eeg_i = find(strcmp('eeg',{points.name}));
eeg_points = points(eeg_i);
eeg0_i = find(cell2mat({eeg_points.num})==0);

if isempty(eeg0_i)
  fprintf('%s: reference EEG not found\n',mfilename);
  return;
end;

ref_EEG_coords=[eeg_points(eeg0_i).x,eeg_points(eeg0_i).y,eeg_points(eeg0_i).z];


