function ts_plot_fields(avg_data,varargin);
% ts_plot_fields - plots MEG fields on helmet surface using
%   vectorview_field_plot
%   plots multiple conditions in a ring
%
% Usage:
%   ts_plot_fields(avg_data,'key1', value1,...);
%
% equired input:
%  avg_data - average data structure (see avg_fif_data)
%
% Optional parameters:
%  vvfpfile - file name for presaved vectorview_field_plot mat file
%    {default: 'vectorview_field_plot.mat'}
%  badchanfile - name of text file containing bad channel labels
%    {default: []}
%  time0 - start time of averaged data range (msec)
%    {default: 60}
%  time1 - end time of averaged data range (msec)
%    {default: 70}
%  view_angle - 2x1 vector containing azimuth and elevation angles
%    {default: [0 20]}
%  scale_max - max value for color scale (fT ?)
%    {default: 100}
%  usemag_flag - [1|0] toggle whether to use magnetometer data
%    {default: 1}
%  conditions - vector of condition numbers (not event codes) to display
%    {default: []} (if empty, display all in a grid)
%
%  created:       05/17/06   by Don Hagler
%  last modified: 08/18/06   by Don Hagler
%

DEFAULT_VVFPFILE = 'vectorview_field_plot.mat';
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 60;
DEFAULT_TIME1 = 70;
DEFAULT_VIEW_ANGLE = [0 20];
DEFAULT_SCALE_MAX = 70;
DEFAULT_USEMAG_FLAG = 1;

if nargin < 1
  help ts_plot_fields;
  return;
end;

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), opt=struct(options{:}); 
  else opt = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

try, opt.vvfpfile;      catch, opt.vvfpfile = DEFAULT_VVFPFILE; end;
try, opt.badchanfile;   catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;         catch, opt.time0 = DEFAULT_TIME0; end;
try, opt.time1;         catch, opt.time1 = DEFAULT_TIME1; end;
try, opt.view_angle;    catch, opt.view_angle = DEFAULT_VIEW_ANGLE; end;
try, opt.scale_max;     catch, opt.scale_max = DEFAULT_SCALE_MAX; end;
try, opt.usemag_flag;   catch, opt.usemag_flag = DEFAULT_USEMAG_FLAG; end;
try, opt.conditions;    catch, opt.conditions = []; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'vvfpfile' 'badchanfile' 'time0' 'time1' 'view_angle'...
         'scale_max' 'conditions' 'usemag_flag'...
   },;
   otherwise, error(['ts_plot_fields: unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
vvfpfile = opt.vvfpfile;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
view_angle = opt.view_angle;
scale_max = opt.scale_max;
usemag_flag = opt.usemag_flag;
conditions = opt.conditions;
clear opt optfields options;

load(vvfpfile);

grad_chans = find(strncmp('grad',lower({avg_data.sensor_info.typestring}),...
             length('grad')));
mag_chans  = find(strcmp('mag',lower({avg_data.sensor_info.typestring})));
all_conds = [1:length(avg_data.averages)];
if isempty(conditions)
  conditions = all_conds;
end;
conditions = intersect(conditions,all_conds);
if isempty(conditions)
  fprintf('%s: no valid conditions to plot\n',mfilename);
  return;
end;

% set start and end samples
sfreq = avg_data.sfreq;
t_trigger = avg_data.averages(1).time(1)*1000;
t0 = round((time0 - t_trigger)*sfreq/1000);
t1 = round((time1 - t_trigger)*sfreq/1000);

color_range = [-scale_max,scale_max];

% get badchans from avg_data
bad_channels = find(cell2mat({avg_data.sensor_info.badchan})==1);
% read badchan file
labels = {avg_data.sensor_info.label}; 
if ~isempty(badchanfile)
  badchan_i = ts_read_txt_badchans(badchanfile,labels);
else
  badchan_i = [];
end;
bad_channels = unique([bad_channels,badchan_i]);
chans = union(grad_chans,mag_chans); % only use meg
badchan_i = intersect(badchan_i,chans); % exclude non-meg channels from bad chans
chans = setdiff(chans,badchan_i); % exclude bad chans

nconds = length(conditions);
nr = floor(sqrt(nconds));
nc = ceil(nconds/nr);

hold on;
for j=1:nconds
  k=conditions(j);
  subplot(nr,nc,j);
  B = mean(avg_data.averages(k).data(1:306,t0:t1),2);
  chann_index = ones(306,1);
  chann_index(badchan_i) = 0;
  if usemag_flag
    chann_index(mag_chans) = 0;
  end;
  vectorview_field_plot(B,chann_index,color_range,view_angle,vv_field_plot);
  title(sprintf('condition %d',k));
  axis off;
  colorbar off;
end;

%colormap bluehot;
%colormap(colormap_blueblackred);
colormap(colormap_cyanblueblackredyellow);

return;

