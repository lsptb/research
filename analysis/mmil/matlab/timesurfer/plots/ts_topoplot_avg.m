function ts_topoplot_avg(avg_data,varargin);
% ts_topoplot_avg - plots eeg potentials from avg_data structure on cartoon
%   head using fieldtrips's topoplotER
%
% Usage:
%   ts_topoplot_avg(avg_data,'key1', value1,...);
%
% equired input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  condition - condition number (not event code) used to index avg_data.averages
%    {default: 1}
%  chantype - channel type (e.g. 'mag','grad1','grad2','eeg','gradpow','other')
%    'gradpow' is the hypotenuse of grad1 and grad2 pairs (power)
%    {default: 'eeg'}
%  time0 - start time of averaged data range (msec)
%    {default: 50}
%  time1 - end time of averaged data range (msec)
%    {default: 60}
%  scale_max - max value for scale (uVolts, fT, or fT/cm)
%    {default: uses max/min values}
%  'badchanfile' - name of text file containing bad channel labels
%    {default: []}
%
%  created:       06/06/06   by Don Hagler
%  last modified: 07/31/06   by Don Hagler
%

DEFAULT_CHANTYPE = 'eeg';
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 50;
DEFAULT_TIME1 = 60;
DEFAULT_SCALE_MAX = 0;
DEFAULT_CONDITION = 1;

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

try, opt.chantype;      catch, opt.chantype = DEFAULT_CHANTYPE; end;
try, opt.badchanfile;   catch, opt.badchanfile = DEFAULT_BADCHANFILE; end;
try, opt.time0;         catch, opt.time0 = DEFAULT_TIME0; end;
try, opt.time1;         catch, opt.time1 = DEFAULT_TIME1; end;
try, opt.scale_max;     catch, opt.scale_max = DEFAULT_SCALE_MAX; end;
try, opt.condition;     catch, opt.condition = DEFAULT_CONDITION; end;

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'chantype' 'badchanfile' 'time0' 'time1'...
         'scale_max' 'condition'...
   },;
   otherwise, error([mfilename ': unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
chantype = opt.chantype;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
scale_max = opt.scale_max;
condition = opt.condition;
clear opt optfields options;

% convert to field trip format
switch chantype
case {'mag' 'grad1' 'grad2' 'eeg' 'gradpow' 'other'},;
otherwise
  help(mfilename)
  fprintf('%s: unsupported chantype (%s)\n',mfilename);
  return;
end;
if strcmp(chantype,'gradpow')
  FT_grad1_data = ts_avg2fieldtrip(avg_data,'condition',condition,...
    'chantype','grad1','badchanfile',badchanfile);
  FT_grad2_data = ts_avg2fieldtrip(avg_data,'condition',condition,...
    'chantype','grad2','badchanfile',badchanfile);
  FT_data = FT_grad1_data;
  FT_data.avg = sqrt(FT_grad1_data.avg.^2 + FT_grad2_data.avg.^2);
else
  FT_data = ts_avg2fieldtrip(avg_data,'condition',condition,...
    'chantype',chantype,'badchanfile',badchanfile);
end;

if isempty(FT_data)
  return;
end;

switch chantype
case {'grad1' 'grad2' 'gradpow' 'other'}
  units = 10^-13;
case 'mag'
  units = 10^-15;
case 'eeg'
  units = 10^-6;
end; 

% set configuration for topoplot
cfg = [];
cfg.xlim = [time0/1000 time1/1000];
if scale_max > 0
  cfg.zlim = [-scale_max,scale_max]*units;
end;
cfg.fontsize = 1;
cfg.interactive = 'yes';
cfg.colorbar = 'no';
%cfg.showindex   = 'yes';
%cfg.baseline = 'yes';
%cfg.baselinetype = 'relative';
cfg.comment = [];

% plot data
ts_topoplotER(cfg,FT_data);

colormap(cmap_cyanblueblackredyellow);

return;

