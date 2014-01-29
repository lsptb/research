function self = test_convert(self)
    ntools_procIEEG(self.recording_filename_root);
    
    fid_ieeg = fopen('256-sin-waves-30kHz-int16.ieeg.dat', 'r');

    assert_not_equals(fid_ieeg, -1);
    
    fclose(fid_ieeg);
end