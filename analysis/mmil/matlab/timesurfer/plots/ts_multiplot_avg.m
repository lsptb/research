function ts_multiplot_avg(avg_data,varargin);
% ts_multiplot_avg - plots MEG/EEG sensor waveforms
%
% Usage:
%   ts_multiplot_avg(avg_data,'key1', value1,...);
%
% equired input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  'conditions' - vector of condition numbers (not event code)
%    used to index avg_data.averages
%    {default: [1]}
%  'chantype' - channel type (e.g. 'mag','grad1','grad2','eeg','grad','gradpow')
%    'grad' is both grad1 and grad2 channels
%    'gradpow' is the hypotenuse of grad1 and grad2 pairs (power)
%    {default: 'eeg'}
%  'scale_max' - max value for scale (uVolts, fT, or fT/cm)
%    {default: 0 -> will use max/min}
%  'linewidth' - trace line width
%    {default: 1.5}
%  'badchanfile' - name of text file containing bad channel labels
%    {default: []}
%
%  created:       06/06/06   by Don Hagler
%  last modified: 02/04/09   by Jason Sherfey (added catch for 'event')
%

%% todo: make labels an option
%% todo: make a way to show a subset of channels

DEFAULT_CHANTYPE = 'eeg';
DEFAULT_BADCHANFILE = [];
DEFAULT_TIME0 = 50;
DEFAULT_TIME1 = 60;
DEFAULT_SCALE_MAX = 0;
DEFAULT_CONDITIONS = [1];
DEFAULT_LINEWIDTH = 1.5;

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
try, opt.conditions;    catch, opt.conditions = DEFAULT_CONDITIONS; end;
try, opt.linewidth;     catch, opt.linewidth = DEFAULT_LINEWIDTH; end;

if isfield(opt,'event')
  opt.conditions = find(ismember([avg_data.averages.event_code],opt.event));
  opt = rmfield(opt,'event');
end

optfields = fieldnames(opt);
for index=1:length(optfields)
   switch optfields{index}
   case {'chantype' 'badchanfile' 'time0' 'time1'...
         'scale_max' 'conditions' 'linewidth'...
   },;
   otherwise, opt = rmfield(opt,optfields{index}); %error([mfilename ': unrecognized option: ''' optfields{index} '''' ]);
   end;
end;

% discard options structure;
chantype = opt.chantype;
badchanfile = opt.badchanfile;
time0 = opt.time0;
time1 = opt.time1;
scale_max = opt.scale_max;
conditions = opt.conditions;
linewidth = opt.linewidth;
clear opt optfields options;

FT_data = [];
for c=1:length(conditions)
  % convert to field trip format
  switch chantype
  case {'mag' 'grad1' 'grad2' 'eeg' 'gradpow' 'grad'},;
  otherwise
    help(mfilename);
    fprintf('\n%s: unsupported chantype (%s)\n',mfilename);
    return;
  end;
  if strcmp(chantype,'gradpow')
    FT_grad1_data = ts_avg2fieldtrip(avg_data,'condition',conditions(c),...
      'chantype','grad1','badchanfile',badchanfile);
    FT_grad2_data = ts_avg2fieldtrip(avg_data,'condition',conditions(c),...
      'chantype','grad2','badchanfile',badchanfile);
    FT_data{c} = FT_grad1_data;
    FT_data{c}.avg = sqrt(FT_grad1_data.avg.^2 + FT_grad2_data.avg.^2);
  else
    FT_data{c} = ts_avg2fieldtrip(avg_data,'condition',conditions(c),...
      'chantype',chantype,'badchanfile',badchanfile);
  end;
end;

switch chantype
case {'grad1' 'grad2' 'gradpow' 'grad'}
  units = 10^-13;
case 'mag'
  units = 10^-15;
case 'eeg'
  units = 10^-6;
end; 

% set configuration for multiplot
cfg = [];
if scale_max > 0
  cfg.ylim = [-scale_max,scale_max]*units;
end;
cfg.colorbar = 'no';
cfg.showlabels = 'yes';
cfg.linewidth = linewidth;

% plot data
ts_multiplotER(cfg,FT_data{:});

return;

