loadflag  = 0; 
outpath   = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s1';
subject   = 's1';  
fid       = 1;

if loadflag
  findex    = 1:5;
    if strcmp(subject,'s1')
      badlabels = {'C1','Cz','C6','CPz','CP4'};
      fiffiles  = {...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_2_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_3_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_4_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_5_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_6_raw.fif'};% ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_DC_s1_7_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_DC_s1_8_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_9_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_10_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_11_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_12_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_13_raw.fif' ...
  %       '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_14_raw.fif'};
    elseif strcmp(subject,'s8')
      fiffiles  = {...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_1_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_2_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_3_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_4_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_5_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_6_nb01_060808.fif' ...
    %     '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/emptyroom_nb01_060808.fif' ...
        };  
    else
      error('Files have not been specified for %s',subject);
    end      
  % read fif files, convert data to timesurfer format, and save mat files
  matfiles = {};
  chantype = {'grad1','grad2'};
  for f = 1:length(fiffiles)
    fif = fiffiles{f};
    [fpath,fname,fext]  = fileparts(fif);
    outfile             = sprintf('%s/matfiles/%s_grad.mat',outpath,fname);
    matfiles{end+1}     = outfile;
    if exist(outfile,'file') % never overwrite (param independent)
      fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
      continue
    else
      fprintf(fid,'Reading FIF file: %s\n',fif);
    end
    data = ts_loadfif(fif,chantype,'epochs');
    fprintf(fid,'Saving MAT file: %s\n',outfile);
    save(outfile,'data');
    clear fif data
  end
  % read mat files and combine data
  fprintf(fid,'Loading MAT files:\n');
  for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
  data   = SO_combine_matfiles(matfiles(findex));
%     for f = 1:length(fiffiles)
%       fif = fiffiles{f};
%       [fpath,fname,fext]  = fileparts(fif);
%       outfile             = sprintf('%s/matfiles/%s_eeg.mat',outpath,fname);
%       matfiles{end+1}     = outfile;
%       if exist(outfile,'file') % never overwrite (param independent)
%         fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
%         continue
%       else
%         fprintf(fid,'Reading FIF file: %s\n',fif);
%       end
%       data = ts_loadfif(fif,{'eeg'},'epochs');
%       fprintf(fid,'Saving MAT file: %s\n',outfile);
%       save(outfile,'data');
%       clear fif data
%     end    
%     findex = length(fiffiles) + findex;
%     for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
%     eeg    = SO_combine_matfiles(matfiles(findex));    
  % peaks for grad1 & grad2
  load s1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat
  ORIGPeaks = peaks;    
  % data = grad1 & grad2
end

%% PSD

tic
tlim    = [500 2700];%[2000 2500];
sofrq   = [.01 4];
pxfrq   = [.1 110];       % [Hz] bandpass freqs before periodogram calc
pad     = 250;            % [ms] epoch padding
foi     = 10:2:100;          % [Hz] freqs for periodogram
window  = 'hann';         % window function for periodogram
lnfilt  = 'no';
% smtflag = 1;            % whether to smooth gamma envelope before epoching
% npts    = 6;           % # points for smoothing the gamma envelope
% blcflag = 1;            % whether to baseline correct gamma epochs before averaging
% blcwin  = [-.1 .1];
detflag = 0;
pairs   = 1;

labels  = {data.sensor_info(1:21).label}; %{'MEG 0143','MEG 0142','MEG 0213','MEG 0212'};%,'MEG 0713','MEG 0712','MEG 0723','MEG 0722'};
fignum  = 1;
clear ALLMPSD ALLPSD
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
%     sodat = ts_preproc(ts_data_selection(data,'chanlabel',label),'bpfilter','yes','bpfreq',sofrq,'lnfilter','yes','bandpass_detrend_flag',1);
%     epoch = SO_epochs(sodat,{pks},pad);
%     avg   = SO_average(epoch,npks);

    dat   = ts_preproc(ts_data_selection(data,'chanlabel',label),'bpfilter','yes','bpfreq',pxfrq,'lnfilter',lnfilt,'bandpass_detrend_flag',1);
    epo   = SO_epochs(dat,{pks},pad);
    psd   = calc_periodogram(epo,foi,window); psd.peaktype  = ptype; psd.epochs.event_code    = type;
    mpsd  = SO_average(psd,npks);             mpsd.peaktype = ptype; mpsd.averages.event_code = type;
    
    PSD(type)   = psd;
    MPSD(type)  = mpsd;
    clear psd mpsd
    
%     gamma = ts_preproc(ts_data_selection(data,'chanlabel',label),'bpfilter','yes','bpfreq',gmfrq,'lnfilter','yes','bandpass_detrend_flag',1);
%     gamma = calc_analytic_amplitude(gamma);
%     if smtflag, gamma.epochs.data = smooth(gamma.epochs.data',npts,'lowess')'; end
%     gamma = SO_epochs(gamma,{pks},pad);
%     if blcflag, gamma = ts_preproc(gamma,'blc','yes','blcwindow',blcwin); end
%     gamma = SO_average(gamma,npks);
% %     if smtflag, gamma.averages.data = smooth(gamma.averages.data',npts,'lowess')'; end
% %     gamma = calc_analytic_amplitude(gamma);
%     if smtflag, gamma.averages.data = smooth(gamma.averages.data',npts,'lowess')'; end        
%     if blcflag
%       gamma = ts_preproc(gamma,'blc','yes','blcwindow',[]); 
%       avg   = ts_preproc(avg,'blc','yes','blcwindow',[]); 
% %     end
% %     gamma = calc_analytic_amplitude(gamma);
% %     if smtflag, gamma.averages.data = smooth(gamma.averages.data',npts,'lowess')'; end    
%     if strcmp(ptype,'negpeak')
%       tmpavg = avg;
%       tmpavg.averages.data = -avg.averages.data;
% %       if blcflag, 
%         tmpavg = ts_preproc(tmpavg,'blc','yes','blcwindow',[]); 
% %       end
%       [c,lags] = xcorr(tmpavg.averages.data',gamma.averages.data','coeff');
%       clear tmpavg
%     else
%       [c,lags] = xcorr(avg.averages.data',gamma.averages.data','coeff');
%     end
% 
%     figure(fignum)
%     subplot(3,1,1),plot(avg.averages.time,avg.averages.data), axis tight, title(sprintf('%s: SO-locked average (n=%g)',label,npks)); vline(0,'k'); hline(0,'k');
%     subplot(3,1,2),plot(gamma.averages.time,gamma.averages.data), axis tight, title(sprintf('%gms-smoothed, blc gamma envelope (%g-%gHz)',round(1000*npts/data.sfreq),gmfrq));
%     set(gca,'ylim',[min(0,1.1*min(gamma.averages.data)) 1.1*max(gamma.averages.data)]); vline(0,'k'); hline(0,'k');
%     xlabel('latency (s)');
%     subplot(3,1,3),plot(lags,c); hline(0,'k'); vline(lags(round(length(lags)/2)),'k');
%     truncix = find(lags>=-pad/2 &  lags<=pad/2);
%     ctrunc  = c(truncix);
%     ctrunc  = smooth(ctrunc,2*npts,'lowess');
%     cpeaks  = crossing(diff(ctrunc)); % find peaks
%     cpeaks  = cpeaks(ctrunc(cpeaks)>0); % only keep maxima
%     cpeaks  = cpeaks(ctrunc(cpeaks)>=.5*max(ctrunc(cpeaks))); % select maxima greater than 25% max
%     cpeaks  = cpeaks(nearest(lags(truncix(cpeaks)),0)); % find maxima nearest to 0
%     maxcorr = lags(truncix(cpeaks));%lags(c==max(c)); %maxcorr = lags(abs(c)==max(abs(c)));
%     vline(maxcorr,'k'); title(sprintf('normalized cross correlation (peak at %g)',maxcorr));
%     xlabel('lag (ms)'); axis tight; % set(gca,'ylim',[-.75 .75])
%     fignum = fignum + 1;
  end
   
  warning('off','MATLAB:divideByZero');
  ALLPSD(l)   = ts_combine_data(PSD);
  ALLMPSD(l)  = ts_combine_data(MPSD);
  warning('on','MATLAB:divideByZero');
  clear PSD MPSD

  figure; 
  frq     = ALLMPSD(l).averages(1).time;
  logfrq  = 10*log(frq);
  pow1    = ALLMPSD(l).averages(1).data;
  pow2    = ALLMPSD(l).averages(2).data;
  logpow1 = 10*log10(ALLMPSD(l).averages(1).data);
  logpow2 = 10*log10(ALLMPSD(l).averages(2).data);
  subplot(1,2,1),plot(frq,logpow1,'r',frq,logpow2,'b'); xlabel('freq'); ylabel('dB'); axis tight; 
  title(sprintf('%s (n=%g, peak+/-%gms)',label,npks,pad))
  subplot(1,2,2),plot(logfrq,logpow1,'r',logfrq,logpow2,'b'); xlabel('log freq'); axis tight
  legend(ptypes{1},ptypes{2})
  
  flim = [30 50];
  fprintf('%s, %g to %gHz (neg-pos): %g dB\n',label,flim,mean(logpow1(frq>=flim(1)&frq<=flim(2)))-mean(logpow2(frq>=flim(1)&frq<=flim(2))));
  fprintf('%s, %g to %gHz (neg-pos): %g (abs)\n',label,flim,mean(pow1(frq>=flim(1)&frq<=flim(2)))-mean(pow2(frq>=flim(1)&frq<=flim(2))));
  flim = [20 60];
  fprintf('%s, %g to %gHz (neg-pos): %g dB\n',label,flim,mean(logpow1(frq>=flim(1)&frq<=flim(2)))-mean(logpow2(frq>=flim(1)&frq<=flim(2))));
  fprintf('%s, %g to %gHz (neg-pos): %g (abs)\n',label,flim,mean(pow1(frq>=flim(1)&frq<=flim(2)))-mean(pow2(frq>=flim(1)&frq<=flim(2))));
  flim = [15 80];
  fprintf('%s, %g to %gHz (neg-pos): %g dB\n',label,flim,mean(logpow1(frq>=flim(1)&frq<=flim(2)))-mean(logpow2(frq>=flim(1)&frq<=flim(2))));
  fprintf('%s, %g to %gHz (neg-pos): %g (abs)\n',label,flim,mean(pow1(frq>=flim(1)&frq<=flim(2)))-mean(pow2(frq>=flim(1)&frq<=flim(2))));
  flim = [10 100];
  fprintf('%s, %g to %gHz (neg-pos): %g dB\n',label,flim,mean(logpow1(frq>=flim(1)&frq<=flim(2)))-mean(logpow2(frq>=flim(1)&frq<=flim(2))));
  fprintf('%s, %g to %gHz (neg-pos): %g (abs)\n',label,flim,mean(pow1(frq>=flim(1)&frq<=flim(2)))-mean(pow2(frq>=flim(1)&frq<=flim(2))));
  
end
toc

visualizer(ALLPSD(end)); % view FFT power epochs

% % tic
% parms.toilim   = tlim;
% grad_flip = calc_flip_matrix(parms,data,ORIGPeaks);
% % toc % 50min
% 
% % peaks=ORIGPeaks; save('s1_SO_flip_matrix_filt0.01-4Hz_toi500-2700_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_22-Jun-2010.mat','flip','parms','peaks');
% 
% 
% parms.refindex = 1;
% tmpdat  = ts_preproc(data,'bpfilter','yes','bpfreq',[.01 4],'bandpass_detrend_flag',0);
% newflip = calc_flip_matrix(parms,tmpdat,ORIGPeaks);
% 
% % [neg closer to lag=0] => [pos-DOWN] = [not neg-DOWN] = -1
% negDOWN = [-1 1 -1 -1 -1 -1 -1 1 1 1 1 -1 -1 1 1 -1 1 1 1 1 1];
% % pos/neg-flip wrt 1 ref (MEG 0113)
% refind  = 1;
% flipmat = newflip;
% for k = 1:2
%   flipvec = flipmat(k).matrix(refind,1:length(negDOWN));
%   flipvec = negDOWN(refind)*flipvec;
%   matches = bsxfun(@eq,negDOWN,flipvec);
%   fprintf('%g%% agreement (%s)\n',agree(refind,k),flipmat(k).peaktype);
% end
% 
% % pos/neg flip wrt all refs
% refs  = 1:length(negDOWN);
% agree = zeros(length(refs),2);
% flipmat = grad_flip;
% for refind = refs
%   for k = 1:2
%     flipvec = flipmat(k).matrix(refind,1:length(negDOWN));
%     flipvec = negDOWN(refind)*flipvec;
%     matches = bsxfun(@eq,negDOWN,flipvec);
%     agree(refind,k) = 100*sum(matches)/length(flipvec);
% %     fprintf('%g%% agreement (%s)\n',agree(refind,k),flipmat(k).peaktype);
%   end
% end
% {flipmat(1).sensor_info(find(agree(:,1)>70 | agree(:,2)>70)).label}
% 
% npos = cellfun(@length,{ORIGPeaks(1:21).pospeak});
% nneg = cellfun(@length,{ORIGPeaks(1:21).negpeak});
% % [agree(:,1) npos' agree(:,2) nneg']






