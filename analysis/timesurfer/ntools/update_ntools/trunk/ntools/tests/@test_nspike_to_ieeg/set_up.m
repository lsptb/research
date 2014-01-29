function self = set_up(self)
    global experiment
    load('test-sin-waves-experiment');
    
    self.NUM_CHANNELS = 256;
    self.recording_filename_root = '256-sin-waves-30kHz-int16.nspike.dat';

    self.fid_nspike = fopen(self.recording_filename_root, 'r');
    %self.fid_ieeg = fopen('256-sin-waves-30kHz-int16.ieeg.dat', 'r');
    self.orig_data = load('int16-sin-waves');
end
