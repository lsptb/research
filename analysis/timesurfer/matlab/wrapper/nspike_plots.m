function nspike_plots(data,varargin)
% function nspike_plots(data,varargin)
%
% Purpose: 
% To generate online plots of average waveforms.
%
% Required input: 
%   epoch_data      -- may contain multiple conditions
%
% Optional input:
%   channels        -- vector of channel indices to be displayed
%   cfg             -- configuration structure with preprocessing options
%
% The configuration can contain
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see CHANNELSELECTION for details
%   cfg.trials             = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.latency            = [begin end] in seconds, or 'minperlength', 'maxperlength', 
%                            'prestim', 'poststim' (default = 'maxperlength')
%   cfg.lpfilter      = 'no' or 'yes'  lowpass filter
%   cfg.hpfilter      = 'no' or 'yes'  highpass filter
%   cfg.bpfilter      = 'no' or 'yes'  bandpass filter
%   cfg.lnfilter      = 'no' or 'yes'  line noise removal using notch filter
%   cfg.dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   cfg.medianfilter  = 'no' or 'yes'  jump preserving median filter
%   cfg.lpfreq        = lowpass  frequency in Hz
%   cfg.hpfreq        = highpass frequency in Hz
%   cfg.bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   cfg.lnfreq        = line noise frequency in Hz, default 50Hz
%   cfg.lpfiltord     = lowpass  filter order
%   cfg.hpfiltord     = highpass filter order
%   cfg.bpfiltord     = bandpass filter order
%   cfg.lnfiltord     = line noise notch filter order
%   cfg.medianfiltord = length of median filter
%   cfg.lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.lnfilttype    = digital filter type, 'but' (default) or 'fir'
%   cfg.detrend       = 'no' or 'yes'
%   cfg.blc           = 'no' or 'yes'
%   cfg.blcwindow     = [begin end] in seconds, the default is the complete trial
%   cfg.hilbert       = 'no' or 'yes'
%   cfg.rectify       = 'no' or 'yes'
%   cfg.implicitref   = 'label' or empty, add the implicit EEG reference as zeros
%   cfg.boxcar        = 'no' or number in seconds
%   cfg.derivative    = 'no' or 'yes', compute the derivative
%   cfg.absdiff       = 'no' or 'yes', compute absolute value of derivative
%
% Created by Jason Sherfey on 29-Oct-2008

if nargin==2 || nargin==3
    if isstruct(varargin{1})
        cfg         = varargin{1}; 
    elseif isnumeric(varargin{1})
        channels    = varargin{1};
    end        
end
if nargin==3
    if isstruct(varargin{2})
        cfg         = varargin{2}; 
    elseif isnumeric(varargin{2})
        channels    = varargin{2};
    end        
end

if ~exist('cfg','var'),         
    cfg=[];
    cfg.keeptrials='no';
    cfg.detrend='yes';
    cfg.blc='yes';
end

if isfield(cfg,'opt')
	opt = getfield(cfg,'opt');
	cfg = rmfield(cfg,'opt');
	if isfield(opt,'conditions'), 
		conditions=opt.conditions; 
		opt = rmfield(opt,'conditions');
	end;
	optargs = mmil_parms2args(opt);
else 
	optargs = {};
end;

if ~exist('conditions','var'),		conditions=1:length(data.epochs); end
if ~exist('channels','var'),    channels=[1:data.num_sensors];    end

for i=1:length(conditions)
    fprintf ('Preprocessing and averaging event: %d\n',conditions(i));  
	ftdat     = ts_data2fieldtrip(data,'condition',conditions(i),'channels',channels,'dimord','chan_time');
	ftavg{i}  = timelockanalysis(cfg,ftdat);
end
avg_data = ts_fieldtrip2data(ftavg,'averages',data);

%% Plotting avg_data
args = {};
badchanfile = [];
alpha = [];

args{end+1} = 'channels';
args{end+1} = {avg_data.sensor_info(channels).label};

badchans = [setdiff(1:avg_data.num_sensors,channels)];
args{end+1} = 'badchans';
args{end+1} = {avg_data.sensor_info(badchans).label};

plot_cond = 1:length(avg_data.averages);
plotevents = [avg_data.averages(plot_cond).event_code];

args{end+1} = 'conditions';
args{end+1} = plot_cond;


%% Set up titles, legends, ect.

for j=1:length(plotevents)
    condition_name{j}   = num2str(avg_data.averages(plot_cond(j)).event_code);
    condition(j)        = {num2str(avg_data.averages(plot_cond(j)).event_code)};  
end

foot_note = sprintf('Sampling Rate: %d\n',avg_data.sfreq);               

args{end+1} = 'conditionnames';
args{end+1} = condition_name;

args{end+1} = 'footnotes';
args{end+1} = foot_note;

f_plot='subject_task_';
for j=1:length(condition);
    f_plot = strcat(f_plot,'event',condition{j},'_');
end

tag = '';

args{end+1} = 'figname';
args{end+1} = strcat(f_plot,tag,'avgs.eps');

args{end+1} = 'title';
args{end+1} = sprintf('Average Data Plot');

args_list = args;
args = {avg_data args_list{:}};

ts_iEEG_Plot_Waveforms (args{:},optargs{:});






