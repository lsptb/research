
% configures for hospital hardware
%
%
hw_defn_nyusom;

% global settings
%
%

% default recording path
experiment.recording.recording_path = '/mnt/raid/adam';

% default sample rate for saving
experiment.recording.sample_rate = 30000;

% feature settings
%
%

% enable data glove support
experiment.data_glove.enable = 1;
% number of nspike 30khz clock ticks between dataglove samples
experiment.data_glove.ticks_between_samples = 300;

% channel configuration
%
%

% configures channels to map 1-256 -> 1-256
nview_configure_channels_256;

for i = 1:256
    experiment.channels(i).name = num2str(i);
    experiment.channels(i).lowpass_cutoff = 11000;
    experiment.channels(i).highpass_cutoff = 1;
end


% event support
%
%

experiment.event(1).name         = 'Cue Close Hand';
experiment.event(1).start_offset = 200;
experiment.event(1).end_offset   = 1500;
experiment.event(1).source       = 'dio';
experiment.event(1).source_port  = 'D';
experiment.event(1).source_code  = 6;

experiment.event(2).name         = 'Cue Open Hand';
experiment.event(2).start_offset = 200;
experiment.event(2).end_offset   = 1500;
experiment.event(2).source       = 'dio';
experiment.event(2).source_port  = 'D';
experiment.event(2).source_code  = 8;
