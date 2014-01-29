% load /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s1/s1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat; init_peaks = peaks;
% 
peaks     = init_peaks;
peaktype  = 'pospeak';
w         = .1;           % sec, window size
s         = .01;%1/peaks(1).sfreq;          % sec, step size

Fs    = peaks(1).sfreq;
t     = peaks(1).tstart:1/Fs:peaks(1).tstop;
tshift= nearest(t,data.epochs.time(1));
for k = 1:length(peaks)
  pix = t(peaks(k).pospeak)>=data.epochs.time(1) & t(peaks(k).pospeak)<=data.epochs.time(end);
  peaks(k).pospeak = peaks(k).pospeak(pix) - tshift;
  nix = t(peaks(k).negpeak)>=data.epochs.time(1) & t(peaks(k).negpeak)<=data.epochs.time(end);
  peaks(k).negpeak = peaks(k).negpeak(nix) - tshift;
  peaks(k).tstart  = data.epochs.time(1);
  peaks(k).tstop   = data.epochs.time(end);
end
t     = peaks(1).tstart:1/Fs:peaks(1).tstop;
pks   = {peaks.(peaktype)};   
pks   = [pks{:}];
pks   = sort(pks);        % peak indices (sorted)
tk    = t(pks);           % peak times (sorted)
L     = floor(w*Fs/2);    % padding around window centers
n     = floor(s*Fs);      % step size in indices
nsmp  = length(t);        % # of time points
npks  = length(tk);       % # of peaks

c     = L+1:n:nsmp-L;     % window centers in indices
tc    = t(c);             % window centers in sec
nstep = length(c);
count = zeros(1,nstep);
for k = 1:nstep
  count(k) = length(find((pks>=c(k)-L)&(pks<=c(k)+L)));
end

npt         = 100;
optlevel    = 1;
ind         = data.num_sensors+1;
smoothcount = smooth(count,npt);
optima      = crossing(diff(smoothcount));
noptima     = length(optima);
thresh      = mean(smoothcount)+std(smoothcount); 
ignore      = smoothcount < thresh;

% figure
% t0  = 1100;
% tf  = 1200;
% ind = tc>=t0 & tc<=tf;
% subplot(3,1,1),plot(tc(ind),count(ind),'-'), axis tight; hline(thresh,'k')
% subplot(3,1,2),plot(tc(ind),smooth(count(ind),50),'-'), axis tight; hline(thresh,'k')
% subplot(3,1,3),plot(tc(ind),smooth(count(ind),100),'-'), axis tight; hline(thresh,'k')

dat   = ts_preproc(data,'bpfilter','yes','bpfreq',parms.bpfreq,'bandpass_detrend_flag',0);
dat.epochs.time = data.epochs.time(c);
dat.epochs.data = dat.epochs.data(:,c);

dat.sensor_info(ind)          = dat.sensor_info(1);
dat.sensor_info(ind).label    = 'count';
dat.epochs.data(ind,:)        = smoothcount / max(smoothcount);

dat.sensor_info(ind+1)        = dat.sensor_info(1);
dat.sensor_info(ind+1).label  = 'countopt';
dat.epochs.data(ind+1,:)      = smoothcount / max(smoothcount);
dat.epochs.data(ind+1,optima) = optlevel;

dat.sensor_info(ind+2)        = dat.sensor_info(1);
dat.sensor_info(ind+2).label  = 'countthresh';
dat.epochs.data(ind+2,:)      = smoothcount / max(smoothcount);
dat.epochs.data(ind+2,optima) = optlevel;
dat.epochs.data(ind+2,ignore) = 0;
dat.num_sensors               = length(dat.sensor_info);

dat.sensor_info(ind+3)          = dat.sensor_info(1);
dat.sensor_info(ind+3).label    = 'rawcount';
dat.epochs.data(ind+3,:)        = count;

thresh      = mean(smoothcount); 
ignore      = smoothcount < thresh;
dat.sensor_info(ind+4)        = dat.sensor_info(1);
dat.sensor_info(ind+4).label  = 'countthresh2';
dat.epochs.data(ind+4,:)      = smoothcount / max(smoothcount);
dat.epochs.data(ind+4,optima) = optlevel;
dat.epochs.data(ind+4,ignore) = 0;
dat.num_sensors               = length(dat.sensor_info);

visualizer(dat)

% save('s1_1-6_paired-peaks_SO_validation.mat','init_peaks','count');
w      = 30;
colors = 'yrkgc'; L = [.5 3 3 .5 .5];
stage2 = [1 38 60 92 99 108 158 161 164 186];
for k =1:length(stage2), vline(stage2(k)*w,'Color',colors(1),'LineWidth',L(1)); end
stage3 = [12 17 56 61 64 76 157 160 163 167];
for k =1:length(stage3), vline(stage3(k)*w,'Color',colors(2),'LineWidth',L(2)); end
stage4 = [16 18 63 75];
for k =1:length(stage4), vline(stage4(k)*w,'Color',colors(3),'LineWidth',L(3)); end
stageW = [93 182];
for k =1:length(stageW), vline(stageW(k)*w,'Color',colors(4),'LineWidth',L(4)); end
stageR = 105;
for k =1:length(stageR), vline(stageR(k)*w,'Color',colors(5),'LineWidth',L(5)); end




