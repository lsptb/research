% data = grad1 & grad2

% peaks for grad1 & grad2
% load s1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat
% ORIGPeaks = peaks;
% tic
tlim  = [2000 2500];
pad   = 1000;           % [ms] epoch padding
sofrq = [.01 4];
gmfrq = [10 100];
smtflag = 1;            % whether to smooth gamma envelope before epoching
npts    = 6;           % # points for smoothing the gamma envelope
blcflag = 1;            % whether to baseline correct gamma epochs before averaging
blcwin  = [-.1 .1];
detflag = 0;
pairs   = 1;

labels  = {data.sensor_info(1:21).label}; %{'MEG 0143','MEG 0142','MEG 0213','MEG 0212'};%,'MEG 0713','MEG 0712','MEG 0723','MEG 0722'};
fignum  = 1;
for l = 1:length(labels)
  label = labels{l};
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
    peaks = SO_detection(ts_data_selection(data,'chanlabel',label),parms);
  else
    peaks = ORIGPeaks(strmatch(label,{ORIGPeaks.label}));
  end
  if pairs
    peaks = select_peakpairs(peaks,1);
  end
  ref    = strmatch(label,{peaks.label});
  ptypes = {'negpeak','pospeak'};
  for type = 1:length(ptypes)
    ptype = ptypes{type};
    pks   = peaks(ref).(ptype);
    t     = peaks(ref).tstart:1/peaks(ref).sfreq:peaks(ref).tstop;
    pks   = pks(t(pks)>=tlim(1) & t(pks)<=tlim(2));
    npks  = length(pks);
    sodat = ts_preproc(ts_data_selection(data,'chanlabel',label),'bpfilter','yes','bpfreq',sofrq,'lnfilter','yes','bandpass_detrend_flag',1);
    epoch = SO_epochs(sodat,{pks},pad);
    avg   = SO_average(epoch,npks);

    gamma = ts_preproc(ts_data_selection(data,'chanlabel',label),'bpfilter','yes','bpfreq',gmfrq,'lnfilter','yes','bandpass_detrend_flag',1);
    gamma = calc_analytic_amplitude(gamma);
    if smtflag, gamma.epochs.data = smooth(gamma.epochs.data',npts,'lowess')'; end
    gamma = SO_epochs(gamma,{pks},pad);
    if blcflag, gamma = ts_preproc(gamma,'blc','yes','blcwindow',blcwin); end
    gamma = SO_average(gamma,npks);
%     if smtflag, gamma.averages.data = smooth(gamma.averages.data',npts,'lowess')'; end
%     gamma = calc_analytic_amplitude(gamma);
    if smtflag, gamma.averages.data = smooth(gamma.averages.data',npts,'lowess')'; end        
%     if blcflag
      gamma = ts_preproc(gamma,'blc','yes','blcwindow',[]); 
      avg   = ts_preproc(avg,'blc','yes','blcwindow',[]); 
%     end
%     gamma = calc_analytic_amplitude(gamma);
%     if smtflag, gamma.averages.data = smooth(gamma.averages.data',npts,'lowess')'; end    
    if strcmp(ptype,'negpeak')
      tmpavg = avg;
      tmpavg.averages.data = -avg.averages.data;
%       if blcflag, 
        tmpavg = ts_preproc(tmpavg,'blc','yes','blcwindow',[]); 
%       end
      [c,lags] = xcorr(tmpavg.averages.data',gamma.averages.data','coeff');
      clear tmpavg
    else
      [c,lags] = xcorr(avg.averages.data',gamma.averages.data','coeff');
    end

    figure(fignum)
    subplot(3,1,1),plot(avg.averages.time,avg.averages.data), axis tight, title(sprintf('%s: SO-locked average (n=%g)',label,npks)); vline(0,'k'); hline(0,'k');
    subplot(3,1,2),plot(gamma.averages.time,gamma.averages.data), axis tight, title(sprintf('%gms-smoothed, blc gamma envelope (%g-%gHz)',round(1000*npts/data.sfreq),gmfrq));
    set(gca,'ylim',[min(0,1.1*min(gamma.averages.data)) 1.1*max(gamma.averages.data)]); vline(0,'k'); hline(0,'k');
    xlabel('latency (s)');
    subplot(3,1,3),plot(lags,c); hline(0,'k'); vline(lags(round(length(lags)/2)),'k');
    truncix = find(lags>=-pad/2 &  lags<=pad/2);
    ctrunc  = c(truncix);
    ctrunc  = smooth(ctrunc,2*npts,'lowess');
    cpeaks  = crossing(diff(ctrunc)); % find peaks
    cpeaks  = cpeaks(ctrunc(cpeaks)>0); % only keep maxima
    cpeaks  = cpeaks(ctrunc(cpeaks)>=.5*max(ctrunc(cpeaks))); % select maxima greater than 25% max
    cpeaks  = cpeaks(nearest(lags(truncix(cpeaks)),0)); % find maxima nearest to 0
    maxcorr = lags(truncix(cpeaks));%lags(c==max(c)); %maxcorr = lags(abs(c)==max(abs(c)));
    vline(maxcorr,'k'); title(sprintf('normalized cross correlation (peak at %g)',maxcorr));
    xlabel('lag (ms)'); axis tight; % set(gca,'ylim',[-.75 .75])
    fignum = fignum + 1;
  end
end
% toc
% tic
parms.toilim   = tlim;
grad_flip = calc_flip_matrix(parms,data,ORIGPeaks);
% toc % 50min

% peaks=ORIGPeaks; save('s1_SO_flip_matrix_filt0.01-4Hz_toi500-2700_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_22-Jun-2010.mat','flip','parms','peaks');


parms.refindex = 1;
tmpdat  = ts_preproc(data,'bpfilter','yes','bpfreq',[.01 4],'bandpass_detrend_flag',0);
newflip = calc_flip_matrix(parms,tmpdat,ORIGPeaks);

% [neg closer to lag=0] => [pos-DOWN] = [not neg-DOWN] = -1
negDOWN = [-1 1 -1 -1 -1 -1 -1 1 1 1 1 -1 -1 1 1 -1 1 1 1 1 1];
% pos/neg-flip wrt 1 ref (MEG 0113)
refind  = 1;
flipmat = newflip;
for k = 1:2
  flipvec = flipmat(k).matrix(refind,1:length(negDOWN));
  flipvec = negDOWN(refind)*flipvec;
  matches = bsxfun(@eq,negDOWN,flipvec);
  fprintf('%g%% agreement (%s)\n',agree(refind,k),flipmat(k).peaktype);
end

% pos/neg flip wrt all refs
refs  = 1:length(negDOWN);
agree = zeros(length(refs),2);
flipmat = grad_flip;
for refind = refs
  for k = 1:2
    flipvec = flipmat(k).matrix(refind,1:length(negDOWN));
    flipvec = negDOWN(refind)*flipvec;
    matches = bsxfun(@eq,negDOWN,flipvec);
    agree(refind,k) = 100*sum(matches)/length(flipvec);
%     fprintf('%g%% agreement (%s)\n',agree(refind,k),flipmat(k).peaktype);
  end
end
{flipmat(1).sensor_info(find(agree(:,1)>70 | agree(:,2)>70)).label}

npos = cellfun(@length,{ORIGPeaks(1:21).pospeak});
nneg = cellfun(@length,{ORIGPeaks(1:21).negpeak});
% [agree(:,1) npos' agree(:,2) nneg']






