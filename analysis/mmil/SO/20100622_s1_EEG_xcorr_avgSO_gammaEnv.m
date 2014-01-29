% label = 'FCz';
% ptype = 'negpeak';
% tlim  = [2000 2500]; [500 1100] [1700 2700]
% pad   = 1000;           % [ms] epoch padding
% sofrq = [.01 4];
% gmfrq = [20 100];
% npts  = 50;             % # points for smoothing the gamma envelope

label = 'F3';
ptype = 'negpeak';
tlim  = [500 2700];
pad   = 1000;           % [ms] epoch padding
sofrq = [.01 4];
gmfrq = [30 50];
smtflag = 1;            % whether to smooth gamma envelope before epoching
npts    = 50;           % # points for smoothing the gamma envelope
blcflag = 1;            % whether to baseline correct gamma epochs before averaging
blcwin  = [];%-.4 0];
detflag = 0;
pairs   = 1;

if detflag
  parms.detectionflag     = 0;
  parms.bpfilter          = 1;
  parms.blc               = 'yes';
  parms.decimate          = 0;
  parms.smooth            = 1;
  parms.hilbertpeaks      = 1;
  parms.derivpeaks        = 1;
  parms.onlysinglepeaks   = 0;              % 1-skip cycles w/ multiple peaks; 0-select peak with largest amplitude.
  parms.zerocross         = 1;
  parms.monotonic         = 1;
  parms.gtmedian          = 1;
  parms.return_zerocross  = 0;
  parms.bpfreq            = [.01 4];%       % Hz
  parms.blcwindow         = [];             % sec, [] = entire time series; [begin end]
  parms.decimate_factor   = [];
  parms.smooth_window     = .05;            % sec, window used for moving average smoothing
  parms.zero2zero_limits  = [.25 1];        % sec, time bw zero crossings must be within these limits
  parms.mindist           = 2*min(parms.zero2zero_limits);%.2;             % sec, (0=skip)(minimum distance bw detections of same polarity)
  peaks = SO_detection(ts_data_selection(eeg,'chanlabel',label),parms);
else
  peaks = eeg_peaks;
end
if pairs
  peaks = select_peakpairs(peaks,1);
end
ref   = strmatch(label,{peaks.label});
pks   = peaks(ref).(ptype);
t     = peaks(ref).tstart:1/peaks(ref).sfreq:peaks(ref).tstop;
pks   = pks(t(pks)>=tlim(1) & t(pks)<=tlim(2));
npks  = length(pks);
sodat = ts_preproc(ts_data_selection(eeg,'channels',ref),'bpfilter','yes','bpfreq',sofrq,'bandpass_detrend_flag',0);
epoch = SO_epochs(sodat,{pks},pad);
avg   = SO_average(epoch,npks);

gamma = ts_preproc(ts_data_selection(eeg,'channels',ref),'lnfilter','yes','bpfilter','yes','bpfreq',gmfrq,'bandpass_detrend_flag',0);
gamma = calc_analytic_amplitude(gamma);
if smtflag, gamma.epochs.data = smooth(gamma.epochs.data',npts,'lowess')'; end
gamma = SO_epochs(gamma,{pks},pad);
if blcflag, gamma = ts_preproc(gamma,'blc','yes','blcwindow',blcwin); end
gamma = SO_average(gamma,npks);

if strcmp(ptype,'negpeak')
  [c,lags] = xcorr(avg.averages.data',gamma.averages.data','coeff');
else
  [c,lags] = xcorr(avg.averages.data',gamma.averages.data','coeff');
end

figure
subplot(3,1,1),plot(avg.averages.time,avg.averages.data), axis tight, title(sprintf('%s: SO-locked average (n=%g)',label,npks)); vline(0,'k'); hline(0,'k');
subplot(3,1,2),plot(gamma.averages.time,gamma.averages.data), axis tight, title(sprintf('%gpt-smoothed, blc gamma envelope (%g-%gHz)',npts,gmfrq));
set(gca,'ylim',[min(0,1.1*min(gamma.averages.data)) 1.1*max(gamma.averages.data)]); vline(0,'k'); hline(0,'k');
xlabel('latency (s)');
subplot(3,1,3),plot(lags,c); hline(0,'k'); vline(lags(round(length(lags)/2)),'k');
maxcorr = lags(abs(c)==max(abs(c))); vline(maxcorr,'k'); title(sprintf('normalized cross correlation (peak at %g)',maxcorr));
xlabel('lag (ms)')

