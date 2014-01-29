function self = test_nspike_to_ieeg(name)
    self.NUM_CHANNELS = [];
    self.recording_filename_root = [];
    self.fid_nspike = [];
    self.fid_ieeg = [];
    self.orig_data = [];
    
    tc = test_case(name);
    self = class(self, 'test_nspike_to_ieeg', tc);

end