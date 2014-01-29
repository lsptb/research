function self = set_up(self)
    global experiment
    load('test-sin-waves-experiment');
    load('test-sin-waves-recording');
    
    self.NUM_CHANNELS = 256;
    self.NUM_EVENTS = 2;
    self.SFREQ = 2000;
    self.ieeg_filename = '256-sin-waves-30kHz-int16.ieeg.dat';
    self.recording = recording;
    
    self.fid_ieeg = fopen(self.ieeg_filename, 'r');
    load('int16-sin-waves');
    self.orig_data = data;
end
