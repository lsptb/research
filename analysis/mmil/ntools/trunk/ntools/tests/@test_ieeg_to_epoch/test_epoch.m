function self = test_epoch(self)
    epoch_data = ntools_gen_epoch_data('ieeg', self.recording);

    assert_equals(self.NUM_CHANNELS, epoch_data.num_sensors, 'Should have 256 channels in data');
    assert_equals(self.SFREQ, epoch_data.sfreq, 'Should be 2000 Hz sampling rate');    
    assert_equals(self.NUM_EVENTS, length(epoch_data.epochs), 'Should have two event types');

    assert_equals(1, epoch_data.epochs(1).num_trials, 'Should have exactly 1 event');
    assert_equals(4, epoch_data.epochs(2).num_trials, 'Should have exactly 4 events');

    assert_equals(2000, size(epoch_data.epochs(1).data, 2), 'Should be entire recording, 2000 samples');
    assert_equals(400, size(epoch_data.epochs(2).data, 2), 'Should be 200 ms at 2000 Hz = 400 samples');
end