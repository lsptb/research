function [avg_data] = ts_postprocess_avg(avg_data,varargin);
% [avg_data_post] = ts_postprocess_avg(avg_data,[options]);
%
% Usage:
%  [avg_data_post] = ts_postprocess_avg(avg_data, 'key1', value1,...);
%
% Required input:
%  avg_data - averaged data structure (output from avg_fif_data.m)
%
% Optional parameters:
%  'stim_delay' - offset applied to time zero (msec)
%     (will be subtracted from time vector)
%     can be used to adjust for stimulus onset delay after trigger
%     { default: 0 }
%  'bandpass'    - [0|1] Toggle bandpass frequency filter
%     { default: 0 }
%  'bandpass_low_cf'  - low cutoff frequency (high-pass filter) (Hz)
%     { default: 0 }
%  'bandpass_low_tb'  - low cutoff transition band (Hz)
%     { default: 0 }
%  'bandpass_high_cf' - high cutoff frequency (low-pass filter) (Hz)
%     { default: 100 }
%  'bandpass_high_tb' - high cutoff transition band (Hz)
%     { default: 0 }
%  'dsfact'  - downsampling factor -- must be an integer
%     data is downsampled to lower sampling frequency
%     e.g. original sampling freq = 1000, dsfact = 4,
%         resulting sampling freq = 250
%     { default: 1 (no downsampling) }
%  'detrend_events' - [0|1] Toggle linear trend removal
%    { default: 1 }
%  'baseline_sub' - [0|1] Toggle subtract mean of baseline period
%     { default: 1 }
%  'baseline_start' - start time of baseline period (msec)
%     relative to time zero
%     { default: -80 }
%  'baseline_end'   - end time of baseline period (msec)
%     relative to time zero
%     { default: -5 }
%  'badchans'    - vector of bad channel numbers
%     { default: [] }
%  'badchanfile' - name of text file containing bad channel labels
%    {default: []}
%  'rm_badchans' - [0|1] Toggle complete removal of badchans
%     from data matrices -- otherwise just set to zero
%     { default: 0 }
%
% Output:
%   avg_data - structure containing the following fields:
%      num_sensors    (int)
%      sfreq          (double)
%      sensor_info    (struct)
%        typestring   (string)
%        label        (string)
%        loc          (4x4 matrix of doubles)
%        badchan      (int: 0 or 1)
%        type         (int)
%        kind         (int)
%        lognum       (int)
%      coor_trans     (struct)
%        device2head  (4x4 matrix of doubles) (or empty)
%        mri2head     (4x4 matrix of doubles) (or empty)
%      averages       (struct)
%        event_code   (int)
%        num_trials   (int)
%        num_rejects  (struct)
%          mag        (int)
%          grad       (int)
%          eeg        (int)
%          eog        (int)
%          manual     (int)
%        time         (vector of doubles) (in msec)
%        data         (matrix of doubles) (num_sensors x num samples)
%      noise          (struct)
%        num_trials   (int)
%        scale_fact   (double)
%        covar        (matrix of doubles) (num_sensors x num_sensors)
%
%
% Note on order of operations:
%  When specified (i.e. bandpass=1), filtering is performed first.
%     If filtering is done, detrending is done first regardless of whether
%     detrend_events=1.
%  Downsampling is then performed (if dsfact>1),
%    followed by detrending (if detrend_events=1).
%    and finally baseline correction (if baseline_sub=1),
%
%
%  created:       04/26/06   by Don Hagler
%  last modified: 03/18/07   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DEF_STIM_DELAY = 0;
DEF_BASELINE_START = -80;
DEF_BASELINE_END = -5;
DEF_BP_LOW_CF = 0;
DEF_BP_LOW_TB = 0;
DEF_BP_HIGH_CF = 100;
DEF_BP_HIGH_TB = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
   help(mfilename);
   return;
end

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

try, g.baseline_sub;    catch, g.baseline_sub   = 1;                  end;
try, g.baseline_start;  catch, g.baseline_start = DEF_BASELINE_START; end;
try, g.baseline_end;    catch, g.baseline_end   = DEF_BASELINE_END;   end;
try, g.detrend_events;  catch, g.detrend_events = 1;                  end;
try, g.stim_delay;      catch, g.stim_delay     = DEF_STIM_DELAY;     end;
try, g.badchans;        catch, g.badchans       = [];                 end;
try, g.badchanfile;     catch, g.badchanfile    = [];                 end;
try, g.rm_badchans;     catch, g.rm_badchans    = 0;                  end;
try, g.bandpass;        catch, g.bandpass       = 0;                  end;
try, g.bandpass_low_cf; catch, g.bandpass_low_cf  = DEF_BP_LOW_CF;    end;
try, g.bandpass_low_tb; catch, g.bandpass_low_tb  = DEF_BP_LOW_TB;    end;
try, g.bandpass_high_cf;catch, g.bandpass_high_cf = DEF_BP_HIGH_CF;   end;
try, g.bandpass_high_tb;catch, g.bandpass_high_tb = DEF_BP_HIGH_TB;   end;
try, g.dsfact;          catch, g.dsfact         = 1;                  end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'stim_delay' 'badchans' 'badchanfile' 'rm_badchans'...
         'detrend_events' 'baseline_start' 'baseline_end' 'baseline_sub'...
         'bandpass' 'bandpass_low_cf' 'bandpass_low_tb'...
         'bandpass_high_cf' 'bandpass_high_tb' 'dsfact'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
baseline_sub = g.baseline_sub;
baseline_start = g.baseline_start;
baseline_end = g.baseline_end;
detrend_events = g.detrend_events;
stim_delay = g.stim_delay;
badchans = g.badchans;
badchanfile = g.badchanfile;
rm_badchans = g.rm_badchans;
bandpass = g.bandpass;
bandpass_low_cf = g.bandpass_low_cf;
bandpass_low_tb = g.bandpass_low_tb;
bandpass_high_cf = g.bandpass_high_cf;
bandpass_high_tb = g.bandpass_high_tb;
dsfact = g.dsfact;

% reject bad parameters
if length(stim_delay)~=1 | ~isnumeric(stim_delay)
  fprintf('%s: stim delay must be a single number\n',mfilename);
  return;
end;
if baseline_sub
  if length(baseline_start)~=1 | ~isnumeric(baseline_start)
    fprintf('%s: baseline_start must be a single number\n',mfilename);
    return;
  end;
  if length(baseline_end)~=1 | ~isnumeric(baseline_end)
    fprintf('%s: baseline_end must be a single number\n',mfilename);
    return;
  end;
end;
if length(dsfact)~=1 | ~isnumeric(dsfact) | dsfact < 1 | ~isint(dsfact)
  fprintf('%s: dsfact must be a single integer >= 1\n',mfilename);
  return;
end;
if bandpass 
  if bandpass_low_cf > bandpass_high_cf
    fprintf('%s: low cutoff freq must be less than high cutoff freq\n',mfilename);
    return;
  end;
  if bandpass_high_tb < 0, bandpass_high_tb=0; end;
  if bandpass_low_tb < 0, bandpass_low_tb=0; end;
end;

% check bad chans are in bounds
if ~isempty(badchans)
  badchans = badchans(union(find(badchans>0),find(badchans<=avg_data.num_sensors)));
end;

% read badchan file
labels = {avg_data.sensor_info.label}; 
if ~isempty(badchanfile)
  badchan_i = ts_read_txt_badchans(badchanfile,labels);
else
  badchan_i = [];
end;
badchans = unique([badchans,badchan_i]);

% make sure there are some good channels left
goodchans = setdiff([1:avg_data.num_sensors],badchans);
if isempty(goodchans)
  fprintf('%s: no good channels... quitting\n',mfilename);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nconds = length(avg_data.averages);

% remove or zero bad chans
if ~isempty(badchans)
  if rm_badchans
    fprintf('%s: removing badchans...\n',mfilename);
    avg_data.sensor_info = avg_data.sensor_info(goodchans);
    avg_data.num_sensors = length(goodchans);
    for j=1:nconds
      avg_data.averages(j).data = avg_data.averages(j).data(goodchans,:);
    end;
    avg_data.noise.covar = avg_data.noise.covar(goodchans,:);
    avg_data.noise.covar = avg_data.noise.covar(:,goodchans);
  else
    fprintf('%s: zeroing badchans...\n',mfilename);
    for b=1:length(badchans)
      k = badchans(b);
      avg_data.sensor_info(k).badchan = 1;
    end;
    for j=1:nconds
      avg_data.averages(j).data(badchans,:) = 0;
    end;
    avg_data.noise.covar(badchans,:) = 0;
    avg_data.noise.covar(:,badchans) = 0;
  end;
end;

% filter
if bandpass
  fprintf('%s: applying bandpass filter...\n',mfilename);
  for j=1:nconds
    % first remove linear trend
    avg_data.averages(j).data = detrend(avg_data.averages(j).data')';
    if bandpass_high_tb > 0 || bandpass_low_tb > 0
      avg_data.averages(j).data = ...
        ts_filtfft_tband(avg_data.averages(j).data,avg_data.sfreq,bandpass_low_cf,...
        bandpass_high_cf,bandpass_low_tb,bandpass_high_tb);
    else
      %% todo: modify filtfft to use tbands
      %   filtfft is based on EEGLAB's eegfiltfft and is ~ 2x faster than
      %   filtfft_tband, which is based on Tao Song's fif_filter_ver2
      avg_data.averages(j).data = ...
        ts_filtfft(avg_data.averages(j).data,avg_data.sfreq,bandpass_low_cf,...
        bandpass_high_cf);
    end;
  end;
end

% adjust times with stim delay
if stim_delay
  for j=1:nconds
    avg_data.averages(j).time = ...
      avg_data.averages(j).time - stim_delay/1000;
  end;
end

% downsample
if dsfact~=1
  fprintf('%s: downsampling...\n',mfilename);
  for j=1:nconds
    avg_data.averages(j).data = resample(avg_data.averages(j).data',1,dsfact,0)';
    avg_data.averages(j).stdev = resample(avg_data.averages(j).stdev',1,dsfact,0)';
    avg_data.averages(j).time = resample(avg_data.averages(j).time,1,dsfact,0);
  end;
  avg_data.sfreq = avg_data.sfreq/dsfact;
end;

% linear trend removal
if detrend_events
  fprintf('%s: removing linear trend...\n',mfilename);
  for j=1:nconds
    avg_data.averages(j).data = detrend(avg_data.averages(j).data')';
  end;
end;

% subtract baseline
if baseline_sub
  fprintf('%s: subtracting baseline...\n',mfilename);
  for j=1:nconds
    % could assume that all conditions have the same timing,
    % but why not be flexible?
    stim_delay = -avg_data.averages(j).time(1)*1000;
    num_samples = length(avg_data.averages(1).time);
    baseline_start_samp = round((baseline_start + stim_delay) * avg_data.sfreq/1000);
    baseline_end_samp   = round((baseline_end + stim_delay) * avg_data.sfreq/1000);
    baseline_start_samp = min(max(1,baseline_start_samp),num_samples);
    baseline_end_samp = min(max(1,baseline_end_samp),num_samples);
    mean_baseline = ...
      mean(avg_data.averages(j).data(:,baseline_start_samp:baseline_end_samp),2);
    avg_data.averages(j).data = ...
      avg_data.averages(j).data - mean_baseline*ones(1,num_samples);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: finished\n',mfilename);
