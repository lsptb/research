
% configures for hospital hardware
%
%
%hw_defn_nyusom;

% global settings
%
%

% default recording path
experiment.recording.recording_path_base = '/mnt/esata/nspike';

% default sample rate for saving
experiment.recording.sample_rate = 30000;

% parameters for processed iEEG data
experiment.processing.ieeg.format = 'float';
experiment.processing.ieeg.byte_size = 4;
experiment.processing.ieeg.sample_rate = 2000;
experiment.processing.ieeg.lowpass = 800;
experiment.processing.ieeg.linefilter.enable = 1;
experiment.processing.ieeg.linefilter.frequencies = [60,120,180,300];

% feature settings
%
%

% disenable data glove support
experiment.data_glove.enable = 0;
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

experiment.event(1).name         = 'Tone';
experiment.event(1).start_offset = -200;
experiment.event(1).stop_offset   = 800;
experiment.event(1).source       = 'dio';
experiment.event(1).source_port  = 'C';
experiment.event(1).source_code  = 255;

