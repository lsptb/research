function [data_out,varargout] = ts_vis_reject (data_in,varargin)
 
% Usage: [epoch_data] = ts_vis_reject (epoch_data,'option',value....
%        [epoch_data,badchannels] = ts_vis_reject (epoch_data,'option',value....
%        [epoch_data,badchannels,badtrials] = ts_vis_reject (epoch_data,'option',value....
%
% Perform rejection processing on epoch data using fieldtrip's rejectvisual
% function.  The output data will contain NaN's for any bad channels across
% all conditions.  The trials selected as bad will be removed from the
% output data set.  It would be necessary to recalculate the noise matrix
% in the calling function based on the noise time period selected by the
% user.
%
% Required input:
%
%  epoch_data - TimeSurfer epoch data structure.
%
% Optional input:
% 
%  chantype- which channels do you want to look at 
%            'mag','grad1','grad2','eeg','grad', 'all'
%            default = 'all' (all channels but 'other')
%  method  - describes how the data should be shown
%            'summary' - show a single number for each channel and trial
%            (default)
%            'channel' - show the data per channel, all trials at once
%            'trial'   - show the data per trial, all channels at once
%  metric  - metric that should be computed in summar mode for each channel
%            in each trial
%            'var'    - variance within each channel (default)
%            'min'    - minimum value in each channel
%            'max'    - maximum value each channel
%            'absmax' - maximum absolute value in each channel
%            'range'  - range from min to max in each channel
%  alim    - determines amplitude scaling for the channel and trial
%            display, if empty then the scaling is automatic (default = [])
%  latency - [begin end] in seconds or 'maxperlength' (default), 'minperlength', 
%            'prestim' ,'poststim'
%  logfile -
%  logid   -
%  
%  Option Preprocessing Options - These preprocessing steps will be applied
%  to data before viewing but will not be applied to the output data.
%
%  bandpass        - 0/1 turn bandpass filter on or off {default = 0}
%  bandpass_low_cf - the value for the high-pass filter {default = 0}
%  bandpass_high_cf- the value for the low-ass filter {default = 100}
%  detrend_events  - 0/1 to detrend the data {default = 0}
%  baseline_sub    - 0/1 to perform baseline correction {default = 0}
%  baseline_start  - start of baesline in ms {default = -80}
%  baseline_end    - end of baseline in ms {default = -5}
%  highpass        - 0/1 to turn on high pass filter {default = 0}
%  low_cf          - value of high pass filter {default = 0}
%  lowpass         - 0/1 to turn on low pass filter {default = 0}
%  high_cf         - value of low pass filter {default = 100}
%  linefilt        - 0/1 to turn on line filter {default = 0}
%  linefreq        - freq of line noise {default = 60}
%
% Required Output: 
%
%   epoch_data  - epoch_data TimeSurfer data structure
%                 trials selected as bad will be removed from data
%                 channels selected as bad will be replaced by NaN's across
%                 all conditions
%
% Optional Output:
%
%   badchannels - index of bad channels
%   badtrials   - cell array the length of the no. of conditions each with a
%                 list of indices for the trials that were removed for a given condition
%                 (indices are in ref. to original data_in)
%
% Created:       11/14/2007 Rajan Patel
% Last Modified: 05/15/2007 Rajan Patel
%
% See also: rejectvisual 


%% Verify Inputs

if (~mmil_check_nargs(nargin, 1))
    return;
end

opt = mmil_args2parms(varargin,...
                      {'chantype','all',{'all','mag','grad1','grad2','eeg','grad'},...  
                       'method','summary',{'summary','channel','trial'},...        % reject visual options
                       'metric','var',{'var','min','max','absmax','range'},...
                       'alim',[],[],...
                       'keepchannel','no',{'no','yes','nan'},...
                       'feedback','no',[],...
                       'latency','maxperlength',[],...
                       'bandpass',0,[],...          % preprocessing options
                       'bandpass_low_cf',0,[],...
                       'bandpass_high_cf',100,[],...
                       'detrend_events',0,[],...
                       'baseline_sub',0,[],...
                       'baseline_start',-80,[],...
                       'baseline_end',-5,[],...
                       'highpass',0,[],...
                       'lowpass',0,[],...
                       'low_cf',[],[],...
                       'high_cf',100,[],...
                       'linefilt',0,[],...
                       'linefreq',60,[],...
                       'channel','all',[],...
                       'logfile',      [],       [],...
                       'logfid',       1,        [] ...
                       'lpfilter','no',{'yes','no'},...
                       'hpfilter','no',{'yes','no'},...
                       'bpfilter','no',{'yes','no'},...
                       'bsfilter','no',{'yes','no'},...
                       'lnfilter','no',{'yes','no'},...
                       'dftfilter','no',{'yes','no'},...
                       'medianfilter','no',{'yes','no'},...
                       'lpfreq',[],[],...
                       'hpfreq',[],[],...
                       'bpfreq',[],[],...
                       'lnfreq',[],[],...
                       'dftfreq',[],[],...              
                       'lpfiltord',[],[],...
                       'hpfiltord',[],[],...
                       'bpfiltord',[],[],...
                       'bsfiltord',[],[],...
                       'lnfiltord',[],[],...
                       'lpfilttype','but',{'but','fir'},...
                       'hpfilttype','but',{'but','fir'},...
                       'bpfilttype','but',{'but','fir'},...
                       'bsfilttype','but',{'but','fir'},...
                       'lpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'hpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'bpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'bsfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
                       'medianfiltord',[],[],...
                       'blc','no',{'yes', 'no'},...
                       'blcwindow',[-inf 0],[],...
                       'detrend','no',{'yes','no'},...     
                       'polyremoval','no',{'yes','no'},...
                       'polyorder',2,[],...
                       'hilbert','no',{'no','abs','complex','real','imag','absreal','absimag','angle'},...
                       'rectify','no',{'yes','no'},...                       
                      },...
                      false);
   
                    
errors = ts_checkdata(data_in);   

if (length(errors)~=0)
 mmil_error(opt, 'Errors supplied data: %s.', sprintf('\t%s\n', errors{:}));
end;
                    
if ~isfield(data_in,'epochs')
 mmil_error(opt,'Input data must be epoch_data.\n');
end                    
    
error(nargoutchk(1,3,nargout,'string'));

reject = rmfield(opt,{'chantype','alim','bandpass','bandpass_low_cf','bandpass_high_cf',...
  'detrend_events','baseline_sub','baseline_start','baseline_end','highpass','lowpass',...
  'low_cf','high_cf','linefilt','linefreq'});

% reject.method      = opt.method;
% reject.metric      = opt.metric;
% reject.alim        = opt.alim;
% reject.keepchannel = 'no';
% reject.latency     = opt.latency;
% reject.channel     = opt.channel;
% reject.feedback    = 'no';

% processing info

if opt.bandpass == 1
 reject.bpfilter = 'yes';
 reject.bpfreq   = [opt.bandpass_low_cf opt.bandpass_high_cf];
end

if opt.detrend_events == 1
  reject.detrend = 'yes';
end

if opt.baseline_sub == 1
  reject.blc       = 'yes';
  reject.blcwindow = [opt.baseline_start opt.baseline_end]/1000;  
end

if opt.highpass == 1
  reject.hpfilter = 'yes';
  reject.hpfreq   = opt.low_cf;
end

if opt.lowpass == 1
  reject.lpfilter = 'yes';
  reject.lpfreq   = opt.high_cf;
end

if opt.linefilt == 1
  reject.lnfilter = 'yes';
  reject.lnfreq   = opt.linefreq;
end

%% Intiallize Data Out

data_out = data_in;
for i = 1:length(data_out.epochs)
  data_out.epochs(i).data = [];   % clear memory
end

badchannels = [];
badtrials   = [];
goodtrials  = [];

%% Visual Rejection

for j=1:length(data_in.epochs)                                                                          % one condition at a time
  mmil_logstr(opt,'Displaying condition %s.',num2str(j));
  ft_epochs = ts_epoch2fieldtrip(data_in,'condition',j,'dimord','chan_time','chantype',opt.chantype); 
  [dum,unproc_idx] = setdiff({data_in.sensor_info.label},ft_epochs.label); 
  unproc_idx       = sort(unproc_idx);                                                                  % index of channels that don't go through rejection
  data_in_badchan  = find([data_in.sensor_info.badchan] == 1);                                          % what channels are already bad
  unproc_idx       = setdiff(unproc_idx,data_in_badchan);                                               % account for already bad channels as being rejected
  reply = 'n'; 
  while ~strcmpi(reply,'Y')
   ft_rej   = rejectvisual(reject,ft_epochs);                                                           % artifact rejection
   reply = input ('Are you happy with your rejection? Y/N: ','s');
   if ~isempty(reply), reply = reply(1); end
   while ~strcmpi(reply,'Y') && ~strcmpi(reply,'N')
     reply = input ('Are you happy with your rejection? Y/N: ','s');
     if ~isempty(reply), reply = reply(1); end
   end
   close% all
  end
  mmil_logstr(opt,'Please wait - Distributing data to TimeSurfer structure.');
  clear ft_epochs
  ft_epochs = ft_rej;
  clear ft_rej
  % distribute data correctly to time surfer data out %
  data_out.epochs(j).num_rejects.manual = data_in.epochs(j).num_rejects.manual + ...
                          (data_in.epochs(j).num_trials - length(ft_epochs.trial));                     % add to manual rejected trials
  data_out.epochs(j).num_trials = length(ft_epochs.trial);                                              % new number of trials
  goodchan_idx = [];
  badchan_idx  = [];
  [dum,goodchan_idx] = intersect ({data_in.sensor_info.label}, ft_epochs.label);                        % channels processed and not bad
  [dum,badchan_idx ] = setdiff   ({data_in.sensor_info.label}, ft_epochs.label);                        % bad channels and those not processed
  goodchan_idx = sort(goodchan_idx);
  badchan_idx  = sort(badchan_idx);                                                                     
  goodtrials{j} = [];
  badtrials{j}  = [];
  for trial = 1:length(ft_epochs.trial)
    data_out.epochs(j).data(goodchan_idx,:,trial) = ft_epochs.trial{trial};
    data_out.epochs(j).data(badchan_idx,:,trial)  = nan;                                                % bad channels and those not processed become nan
    for oldtrial = 1:size(data_in.epochs(j).data,3)                                                     % find the trials that stayed
      if all(data_in.epochs(j).data(goodchan_idx,:,oldtrial) == ft_epochs.trial{trial})
        goodtrials{j} = cat(1,goodtrials{j},oldtrial);                                                  % add a good trial to set of good trials
      end
    end
  end
  badtrials{j}  = setdiff(1:size(data_in.epochs(j).data,3),goodtrials{j});                              % set badtrials
  mmil_logstr(opt,'%s bad trials: %s.',num2str(length(badtrials{j})), num2str(badtrials{j}));
  if ~isempty(unproc_idx)
    data_out.epochs(j).data(unproc_idx,:,:) = data_in.epochs(j).data(unproc_idx,:,goodtrials{j});       % add back the good trials from the channels not processed
    data_in.epochs(j).data = [];
  end
  
  for chans = 1:length(badchan_idx)                                                                         % channels not in field trip structure
    if ~(ismember(badchan_idx(chans),unproc_idx)) && (data_in.sensor_info(badchan_idx(chans)).badchan ~= 1) % if processed and not yet marked as bad
      mmil_logstr(opt,'Marking channel %s as bad.',data_in.sensor_info(badchan_idx(chans)).label);
      data_in.sensor_info(badchan_idx(chans)).badchan  = 1;                                                 % so subsequent conditions have this channel removed during rejection
      data_out.sensor_info(badchan_idx(chans)).badchan = 1;                                                 % mark it in the data out as bad
    end
  end
  data_in.epochs(j).data = []; % clear some memory
  clear goodchan_idx badchan_idx unproc_idx ft_epochs
  % end of putting data back to time surfer %
end

%% Check Bad Channels across Data

mmil_logstr(opt,'Checking that bad channels are consistent across all conditions...');

badchannels = find([data_out.sensor_info.badchan] == 1);                                               % index of bad channels

if ~isempty(badchannels)
  for i = 1:length(data_out.epochs)                                                                    % double check all the conditions and force nan for bad channels
    data_out.epochs(i).data(badchannels,:,:) = nan;
  end
end

%% Recalculate Covariance Matrix

% Can we determine the noise time period based on the data in the
% structure?

%% Assign Variable Outputs

if nargout > 1, varargout{1} = badchannels; end
if nargout > 2, varargout{2} = badtrials;   end
  