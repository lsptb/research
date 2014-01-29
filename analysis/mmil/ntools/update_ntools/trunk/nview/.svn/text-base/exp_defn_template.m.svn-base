
% global settings
%
%

% default sample rate for saving
experiment.recording.sample_rate = 10000;
experiment.recording.comedi_num_channels_to_write = 2;
experiment.recording.nspike_num_channels_to_write_low = 256;
experiment.recording.nspike_num_channels_to_write_high = 0;
experiment.recording.nspike_high_channels_offset = 0;

%  default filtering for highpass and lowpass
for i=1:256
  experiment.channels(i).name = sprintf('iEEG%03d',i);
  experiment.channels(i).lowpass_cutoff = 4000;
  experiment.channels(i).highpass_cutoff = 1;
end

% parameters for recording comedi data
experiment.recording.comedi.sample_rate = 1e4;
experiment.recording.comedi.sample_size = 2;
experiment.recording.comedi.num_channels = 8;
experiment.recording.comedi.sample_format = 'ushort'

experiment.recording.comedi.videosync_channel_index = 1;
experiment.recording.comedi.audio_channel_index = 2;

%parameters for processing comedi data
% for downsampling the saved files
experiment.processing.comedi.sample_rate = 1e4;  
% for subsetting channels in the saved files
experiment.processing.comedi.channels = [1:8];  
% for downsampling the saved files
experiment.processing.comedi.audio.sample_rate = 1e4;  
experiment.processing.comedi.audio.sample_format = 'short';


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


% event support
%
%

experiment.event(1).name         = 'Tone';
experiment.event(1).start_offset = -200;
experiment.event(1).stop_offset   = 800;
experiment.event(1).source       = 'dio';
experiment.event(1).source_port  = 'C';
experiment.event(1).source_code  = 255;

