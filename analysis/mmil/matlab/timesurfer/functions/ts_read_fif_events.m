function [hdr, events, trig_onset, noise_onset] = ts_read_fif_events(filename,trigchan,offset,trig_minduration)
%function [hdr, events, trig_onset] = ts_read_fif_events(filename,[trigchan],[offset],[trig_minduration])
%
% Purpose: reads all events from a Neuromag fif file
%   also reads and returns header info
%
% Usage:
%   [hdr, events, trig_onset] = ts_read_fif_events(filename,trigchan,offset)
%
% Required Parameters:
%   filename: full path of input raw fif file
%
% Optional Parameters
%   trigchan: name of trigger channel with the event codes
%     Can be cell array of trigger channels
%      -- if so, will treat as binary and convert to decimal
%     {default = 'STI101'}
%   offset: event code offset (subtracted from raw trigger values)
%     {default = minimum trigger value}
%   trig_minduration: minimum duration, in samples, of a trigger to be
%     "on".  Any signal seen to be "on" for less than this # samples is 
%     considered noise.
%     {defalt = 5 samples}
%
% Output:
%   hdr: struct containing fif file header info (e.g. sensor information)
%   events: strut array containing events
%     events.type      = string
%     events.latency   = expressed in samples, first sample of file is 1
%     events.condition = numeric event code
%     events.duration  = duration of event, expressed in samples
%   trig_onset: full length vector with the event codes at the trigger onset,
%     zeros otherwise
%   noise_onset: onset latency of trigger changes we deemed as "noise"
%
% NOTE: if "skips" are found in datafile, regular length buffers (i.e. 1 sec)
%       of zeros are inserted; this is reflected in event latencies
% 
% After reading the events structure, you can use the following tricks to
% extract information about those events in which you are interested.
%
% Determine the different events types:
%   unique({events.type})
%
% Get the index of all skip events:
%   find(strcmp('skip', {events.type}))
%
% Make a vector with all triggers:
%   [events(find(strcmp('trigger', {events.type}))).condition]
%
% See also ts_read_fif_data, ts_read_fif_header
%
% based on FieldTrip's read_fcdc_event
% Copyright (C) 2003-2005, F.C. Donders Centre
%
% created:  03/10/06 Don Hagler
% last mod: 08/05/09 Don Hagler
%

% start with an empty events structure
events = [];
trig_onset = [];

if (~mmil_check_nargs(nargin,1)) return; end;

% test whether the file exists
if ~exist(filename,'file')
  error('file %s not found',filename);
end

% set default value of trigchan
if ~exist('trigchan','var') | isempty(trigchan)
  trigchan = 'STI101';
elseif ~isstr(trigchan) && ~iscell(trigchan)
  error('trigchan must be a string or cell array of strings');
end
if ~iscell(trigchan)
  trigchan = {trigchan};
end;

% set default value of offset
if ~exist('offset','var'), offset=0; end;

% set default value of trig_minduration
if (~exist('trig_minduration','var'))
  trig_minduration = 5;
elseif ~isnumeric(trig_minduration)
  error('trig_minduration must be a number');
end;

% read header
hdr = ts_read_fif_header(filename,0); % do not get nBuffs

% check all trigchans if cell array
for i=1:length(trigchan)
  if ~ismember(trigchan{i},hdr.sensors.label)
    error('trigchan %s not found',trigchan{i});
  end;
end

% load trigger data
if length(trigchan)>1
  % loop over all trigchans
  B_all = [];
  for i=1:length(trigchan)
    % read events
    [B_trig,hdr.nBuffs,hdr.tlast,events,sfreq] = ts_read_fif_chan(filename,trigchan(i));
    if i==1
      % put skips in header
      hdr.skips = [];
      if ~isempty(events)
        for j=1:length(events)
          hdr.skips = [hdr.skips events(j).latency];
        end;
      end;
      B_all = zeros(length(B_trig),length(trigchan));
    end;
    B_all(:,i) = (B_trig>0);
  end;
  % convert from binary to decimal
  B_trig = bi2de(B_all)';
else
  % read events
  [B_trig,hdr.nBuffs,hdr.tlast,events,sfreq] = ts_read_fif_chan(filename,trigchan);
  % put skips in header
  hdr.skips = [];
  if ~isempty(events)
    for j=1:length(events)
      hdr.skips = [hdr.skips events(j).latency];
    end;
  end;
end;

% create event structure
trigger_diff = diff(B_trig);
% check for "noise"--trigger changes that
% last less than some # of samples
change_idx = find(trigger_diff~=0);
onset_idx = find(trigger_diff >0);
change_lengths = diff(change_idx);
% Noise is where we see:
% - onset of some value
% - offset of SAME VALUE less than X samples later
% Otherwise, we just let it through; too complex
% to separate noise from fast async trigger value changes
noise_candidate_idx = change_idx(find(change_lengths < trig_minduration));
noise_onset = [];
for i=1:length(noise_candidate_idx)
  trig_onset_idx = noise_candidate_idx(i);
  trig_code      = trigger_diff(trig_onset_idx);
  next_samples   = trigger_diff(trig_onset_idx+1:trig_onset_idx+trig_minduration);
  trig_offset_idx= trig_onset_idx + find(next_samples == -trig_code) - 1;
  
  % Noise found; remove it from our 'differences'
  if (~isempty(trig_offset_idx) )
    trigger_diff(trig_onset_idx) = 0;
    trigger_diff(trig_offset_idx) = 0;
    noise_onset(end+1) = trig_onset_idx;
  end;
end;
% Refresh our 'onset_idx' array now that
% we've stripped out noise
onset_idx = find(trigger_diff >0);
% NOTE: onset_idx indexes the *trigger_diff* matrix.
% That is 1 element smaller than the actual ONSET matrix,
% with that element essentially being taken off the front.
% Thus, we have to add one to get the trig_onset
% and latencies...
trig_onset = onset_idx+1;
events = struct('type', 'trigger', ...
                'duration', 1, ...
                'condition', num2cell(trigger_diff(onset_idx)), ...
                'latency',   num2cell(trig_onset) );
if ~isempty(events)
  % sort the events by latency
  [foo, indx] = sort([events.latency]);
  events = events(indx);
else
  fprintf('%s: WARNING: no events found in %s\n',mfilename,filename);
end
