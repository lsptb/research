function self = test_data_corr_w_original(self)
    recording = self.recording;
    epoch_data = ntools_gen_epoch_data('ieeg', recording);

    upsampled_data = resample(double(epoch_data.epochs(1).data)', 15, 1)';
    for i = 1:self.NUM_CHANNELS
        assert(corr(single(self.orig_data(i,:))', upsampled_data(i,:)') > .95, ...
            sprintf('Corr for sin wav %d should be very high', i));
    end

    % set up time windows
    start_offset = recording.event(2).start_offset / 1000 * self.SFREQ;
    stop_offset = recording.event(2).stop_offset / 1000 * self.SFREQ;
    for i = 1:4
        time_win(i,:) = [(recording.event(2).timestamps(i) + start_offset):(recording.event(2).timestamps(i) + stop_offset)];
    end
    
    % should all be highly correlated
    for i = 1:self.NUM_CHANNELS
        for j = 1:4
            assert(corr(single(self.orig_data(i,time_win(j,:)))', upsampled_data(i,time_win(j,:))') > .95, ...
                sprintf('Corr for sin wav %d should be very high for epoch %d', i, j));
        end
    end

    % at least some in middle-high freqs should have weak correlations
    for i = 1:self.NUM_CHANNELS
        for j = 1:4
            corrs(i,j) = corr(single(self.orig_data(i,time_win(j,:) - 200))', upsampled_data(i,time_win(j,:))');
        end
    end
    assert(any(corrs(:) < .95), 'Some frequencies should have weak correlations when the time window is shifted');
end