pkt       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
%[cellfun(@length,{peaks.pospeak})' cellfun(@length,{peaks.negpeak})']

% skip channels if no paired peaks were detected
skipchans = {peaks(cellfun(@length,{peaks.negpeak})==0).label};
testchans = setdiff(testchans,skipchans);
[sel,jnk] = match_str({peaks.label},testchans);
peaks     = peaks(sel);
nchan     = length(testchans);

for ch    = 1:nchan
  tmp     = peaks(ch).pospeak; peaks(ch).pospeak = tmp(pkt(tmp)>=tlim(1)&pkt(tmp)<=tlim(2));
  tmp     = peaks(ch).negpeak; peaks(ch).negpeak = tmp(pkt(tmp)>=tlim(1)&pkt(tmp)<=tlim(2));
end
[sel,jnk] = match_str({data.sensor_info.label},testchans);
X         = data.epochs.data(sel,:);  % select channel data for complete time series
T         = data.epochs.time;         % complete time vector
data      = ts_data_selection(data,'chanlabel',testchans,'toilim',tlim);
