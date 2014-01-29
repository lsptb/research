function h=ts_fieldplot_avg(avg_data,varargin);
% ts_fieldplot_avg - plots MEG fields from avg_data structure on helmet surface
%   using vectorview_field_plot
%
% Usage:
%   ts_fieldplot_avg(avg_data,'key1', value1,...);
%
% Required input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  condition - condition number (not event code) used to index avg_data.averages
%    {default: 1}
%  vvfpfile - file name for presaved vectorview_field_plot mat file
%    {default: 'vectorview_field_plot.mat'}
%  badchanfile - name of text file containing bad channel labels
%    {default: []}
%  time0 - start time of averaged data range (msec)
%    {default: 50}
%  time1 - end time of averaged data range (msec)
%    {default: 60}
%  view_angle - 2x1 vector containing azimuth and elevation angles
%    {default: [0 20]}
%  scale_max - max value for color scale (fT ?)
%    {default: 100}
%  usemags - ['on'|'off'] toggle whether to use magnetometer data
%    {default: 'on'}
%
%  created:       06/05/06   by Don Hagler
%  last modified: 07/31/06   by Don Hagler
%

DEFAULT_VVFPFILE = 'vectorview_field_plot.mat';
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 60;
DEFAULT_TIME1 = 70;
DEFAULT_VIEW_ANGLE = [0 20];
DEFAULT_SCALE_MAX = 50;
DEFAULT_CONDITION = 1;
DEFAULT_USEMAGS = 'on';

GRAD1 = 1:3:306;
GRAD2 = 2:3:306;
MAG   = 3:3:306;


if nargin < 1
  help(mfilename);
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

try, opt.condition;     catch, opt.condition = DEFAULT_CONDITION; end;
try, opt.vvfpfile;      catch, opt.vvfpfile = DEFAULT_VVFPFILE; end;
try, opt.badchanfile;   catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;         catch, opt.time0 = DEFAULT_TIME0; end;
try, opt.time1;         catch, opt.time1 = DEFAULT_TIME1; end;
try, opt.view_angle;    catch, opt.view_angle = DEFAULT_VIEW_ANGLE; end;
try, opt.scale_max;     catch, opt.scale_max = DEFAULT_SCALE_MAX; end;
try, opt.usemags;       catch, opt.usemags = DEFAULT_USEMAGS; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'vvfpfile' 'badchanfile' 'time0' 'time1' 'view_angle'...
         'scale_max' 'condition' 'usemags'...
   },;
   otherwise, error([mfilename ': unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
vvfpfile = opt.vvfpfile;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
view_angle = opt.view_angle;
scale_max = opt.scale_max;
condition = opt.condition;
usemags = opt.usemags;
clear opt optfields options;

load(vvfpfile);

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
chans = [1:306]; % only use meg
badchan_i = intersect(badchan_i,chans); % exclude non-meg channels from bad chans
chans = setdiff(chans,badchan_i); % exclude bad chans

hold on;
B = mean(avg_data.averages(condition).data(1:306,t0:t1),2);
chann_index = ones(306,1);
chann_index(badchan_i) = 0;
if ~strcmp('usemags','on');
  chann_index(MAG) = 0;
end;
vectorview_field_plot(B,chann_index,color_range,view_angle,vv_field_plot);
axis off;
colorbar off;

colormap(colormap_blueblackred);

return;

