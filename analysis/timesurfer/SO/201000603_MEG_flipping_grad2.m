
%% 31-May-2010. UP/DOWN state detection (Mulkovski et al, Cerebral Cortex, 2007)

data    = SO_combine_matfiles(matfiles);
toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8
data    = ts_data_selection(data,'toilim',toilim);

% labels  = {'MEG0113','MEG0143','MEG0213','MEG0222','MEG0322','MEG0343'};
labels  = {'MEG0113','MEG0122','MEG0132','MEG0143','MEG0213','MEG0222','MEG0313','MEG0322','MEG0343','MEG0413','MEG0422','MEG0432','MEG0443','MEG0723'};
% bandpass filter 20-100Hz
dat     = ts_data_selection(data,'toilim',[900 1350],'chanlabel',labels);
proc    = ts_preproc(dat,'bpfilter','yes','bpfreq',[20 100]);
dat     = ts_data_selection(dat,'toilim',[950 1300]);
proc    = ts_data_selection(proc,'toilim',[950 1300]);

% calculate standard deviation in a running 5- or 10-ms window (=RMS)
StdDevWindow = 10; % ms
SmoothWindow = 50; % ms

Fs = proc.sfreq;
dt = 1/Fs;

% convert windows from ms to indices and force odd
StdDevWindow = round(StdDevWindow/(1000*dt));
SmoothWindow = round(SmoothWindow/(1000*dt));
if ~mod(StdDevWindow,2),StdDevWindow=StdDevWindow+1; end
if ~mod(SmoothWindow,2),SmoothWindow=SmoothWindow+1; end
begntime  = length(proc.epochs.time);
endntime  = begntime - 2*floor(StdDevWindow/2);
beg_ind   = floor(StdDevWindow/2)+1;
end_ind   = begntime - floor(StdDevWindow/2);

tic
% initialize std struct and calculate std in moving window
stdproc = proc;
stdproc.epochs.time = proc.epochs.time(beg_ind:end_ind);
stdproc.epochs.data = zeros(proc.num_sensors,endntime);
for k = 1:endntime
  stdproc.epochs.data(:,k) = std(proc.epochs.data(:,k:k+StdDevWindow-1),0,2);
end
toc
% % very slow
% lo_ind = num2cell(1:end_ind);
% hi_ind = num2cell(beg_ind:begntime);
% stdproc.epochs.data(1,:) = cellfun(@(a,b)std(proc.epochs.data(1,a:b)),lo_ind,hi_ind)

% linear smoothing with a 50-ms window
METHOD = 'lowess';
SPAN   = SmoothWindow;
smoothstd = stdproc;
smoothstd.epochs.data = zeros(size(stdproc.epochs.data));
tic
for k = 1:stdproc.num_sensors
  smoothstd.epochs.data(k,:) = smooth(stdproc.epochs.data(k,:),SPAN,METHOD);
  toc
end

dat = ts_data_selection(dat,'toilim',[smoothstd.epochs.time(1) smoothstd.epochs.time(end)]);
smoothstd.epochs.event_code = 2;
dat(2) = smoothstd;
dat = ts_combine_data(dat);

tmp = ts_preproc(ts_data_selection(dat,'events',1),'bpfilter','yes','bpfreq',[.1 20]);
dat.epochs(1).data = tmp.epochs.data;

    for ch = 1:dat.num_sensors
      % normalize
      dat.epochs(2).data(ch,:) = smoothstd.epochs.data(ch,:) / max(abs(smoothstd.epochs.data(ch,:)));
      % rescale
      dat.epochs(2).data(ch,:) = dat.epochs(2).data(ch,:) * max(abs(dat.epochs(1).data(ch,:)));
    end

visualizer(dat);

%% 03-Jun-2010
tic
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

% SO detection parameters
parms     = [];
parms.fc1 = .3;   % Hz (lower cut-off freq prior to detection)
parms.fc2 = 3;    % Hz (upper cut-off freq prior to detection)
parms.monothresh = -5; % increase to be stricter; rec: stop at fc1
parms.minzero = 300; % ms (minimum distance between zero-crossings)
parms.amp_stdfact = 0; % # std > or < mean to threshold

parms.monophase_flag  = 1;
parms.surround_flag   = 1;
parms.interdet_flag   = 0;
parms.zero_flag       = 1;    
parms.zeroplus_flag   = 1; 
parms.amp_flag        = 1;    
parms.TFrej_flag      = 0;    
parms.preproc_flag    = 1;

% % LOOP OVER SUBJECTS
subjects  = {'s5','s6','s7','s8'};
subj=4;%for subj  = 4%:length(subjects)
  fprintf('Subject %g of [%s]\n',subj,num2str(1:length(subjects)));
  subject = subjects{subj};
% % subject matfiles
% if subj == 1
%   matfiles = {...

  % LOAD DATA AND FIND PEAKS
  subjects  = {'s5','s6','s7','s8'};
  subj = 4;
  fprintf('Subject %g of [%s]\n',subj,num2str(1:length(subjects)));
  subject = subjects{subj};

  grads = {'grad1','grad2'}; graphmarker = 'o+';
%     gradtype=1;%for gradtype = 1:length(grads)
%     fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
%     matfiles = {...
%     ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
%     ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
%     ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
%     };
%     grad1 = SO_combine_matfiles(matfiles);

    gradtype=2;%for gradtype = 1:length(grads)
    fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    grad2 = SO_combine_matfiles(matfiles);
    
parms     = [];
parms.fc1 = .05;   % Hz (lower cut-off freq prior to detection)
parms.fc2 = 5;    % Hz (upper cut-off freq prior to detection)
parms.monothresh = -5; % increase to be stricter; rec: stop at fc1
parms.minzero = 300; % ms (minimum distance between zero-crossings)
parms.amp_stdfact = 0; % # std > or < mean to threshold

parms.monophase_flag  = 1;
parms.surround_flag   = 1;
parms.interdet_flag   = 0;
parms.zero_flag       = 1;    
parms.zeroplus_flag   = 1; 
parms.amp_flag        = 1;    
parms.TFrej_flag      = 0;    
parms.preproc_flag    = 1;
    
toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8
  grad2    = ts_data_selection(grad2,'toilim',toilim);
  fprintf('Time = %g to %g sec\n',toilim);
  nchan       = grad2.num_sensors;
  peaks   = SO_peaks(grad2,parms);
  t       = grad2.epochs.time;
  events  = [];
  for k   = 1:length(peaks)
    events(k).label = grad2.sensor_info(k).label;
    events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
    events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
  end
  outfile = sprintf('%s/s8_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s.mat',outpath,parms.fc1,parms.fc2,toilim,grads{gradtype},'all');%[reflabels{:}]);
  save(outfile,'events','peaks');
  toc
    tic
    chans   = [1 4 5 6 10 12];
    alpha   = 1;
    bpfreq  = [5 115];    
    toi1    = [600 1350];
    toi2    = [610 1340];
    % constant spectral resolution
    foi     = 20:2:100;
    sf      = 2;
    % select data to display with TF results
    dat = ts_data_selection(grad2,'channels',chans,'toilim',toi2,'verbose',0);
    % filter data for display
    dat = ts_preproc(dat,'bpfilter','yes','bpfreq',[.1 10],'bandpass_detrend_flag',0);
    % select channels of interest
    proc    = ts_data_selection(grad2,'channels',chans,'verbose',0);
    % filter out undesired freqs before wavelet
    proc    = ts_preproc(proc,'bpfilter','yes','bpfreq',bpfreq,'bandpass_detrend_flag',0);
    % remove data padding used for filtering
    proc    = ts_data_selection(proc,'toilim',toi1,'verbose',0);
    % Morlet wavelet analysis
    tfdata  = ts_freqanalysis_fieldtrip(proc,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0);
    % remove edge artifacts produced by wavelet analysis
    tfdata  = ts_data_selection(tfdata,'toilim',toi2,'verbose',0);
    % calculate z-scores wrt entire recording
    zdata   = ts_zscore(tfdata,'baselinetype','relchange','blcwindow','all','verbose',0);
%     % correct for 1/f^2 scaling of power spectrum
%     zdata = tfdata;
%     freqs = tfdata.timefreq.frequencies;
%     fcorr = repmat(permute(freqs.^(alpha),[1 3 2]),[zdata.num_sensors length(zdata.timefreq.time) 1]);
%     zdata.timefreq.power = zdata.timefreq.power .* fcorr;
%     clear fcorr
    % average over a frequency band
    freqband= [20 100];
    tfband  = SO_freqband_average(zdata,freqband);
    % normalize band-average and scale wrt the processed time series
    % note: this is necessary for overlaying band-avg & time-domain data
    for ch = 1:tfband.num_sensors
      % de-mean
%       tfband.averages.data(ch,:) = tfband.averages.data(ch,:) - mean(tfband.averages.data(ch,:));
      % normalize
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) / max(abs(tfband.averages.data(ch,:)));
      % rescale
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) * max(abs(dat.epochs(1).data(ch,:)));
      % shift to 0
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) - min(tfband.averages.data(ch,:));
      dat.epochs(1).data(ch,:) = dat.epochs(1).data(ch,:) - min(dat.epochs(1).data(ch,:));      
    end
    % add band-average to the time-domain data structure
    dat.epochs(2)      = dat.epochs(1);
    dat.epochs(2).data = tfband.averages.data;
    % change the event code for the tfband
    dat.epochs(2).event_code = 2;
    toc
    % visualize raw time series, band-average, and TFR power
    visualizer(dat,zdata);
       
    
%% 02-Jun-2010
tic
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

% % LOOP OVER SUBJECTS
subjects  = {'s5','s6','s7','s8'};
subj      = 4;
subject   = subjects{subj};

  % LOAD DATA AND FIND PEAKS
 addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
 outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

  subjects  = {'s5','s6','s7','s8'};
  subj = 4;
  fprintf('Subject %g of [%s]\n',subj,num2str(1:length(subjects)));
  subject = subjects{subj};

  grads = {'grad1','grad2'}; graphmarker = 'o+';
    gradtype = 1;
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    grad1   = SO_combine_matfiles(matfiles); % *.epochs (actually continuous data w/ 1 trial)
    toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8  
    grad1   = ts_data_selection(grad1,'toilim',toilim);    
    peaks   = SO_peaks(grad1,parms);
    
    % filter data for subsequent cluster analysis
    data     = ts_preproc(grad1,'bpfilter','yes','bpfreq',[parms.fc1 parms.fc2],'bandpass_detrend_flag',0);
    t        = data.epochs.time;
    toc
    % Select two channels with opposite SO polarities
%     labels = {'MEG0143','MEG0213'};
    labels  = {'MEG0113','MEG0122','MEG0132','MEG0143','MEG0213','MEG0222','MEG0313','MEG0322','MEG0343','MEG0413','MEG0422','MEG0432','MEG0443','MEG0723'};
    % Initial guess about relationship bw MEG SO polarity & cortical state (based on data inspection):
    % MEG0143-positive during depolarization phase (UP state) - like surface EEG
    % MEG0213-negative during depolarization phase (UP state) - like depth EEG
    
    % Call MEG0143 the reference channel and find the clusters (n=3) to
    % which points belong in both channels at each DOWN time (assume
    % MEG0143 is DOWN when SO phase is negative)
    ncluster = 3;
    reflabel = 'MEG0143';
    refindex = strmatch(reflabel,{peaks.label},'exact');
    % get the times of negative peaks in the reference (MEG0143)
    pktypes  = {'negpeak','pospeak'};
    for typ  = 1:length(pktypes)
      type   = pktypes{typ};
      
      pkindex  = peaks(refindex).(type);
      pktime   = t(pkindex);

      % cluster the reference data
      [tmpidx,tmpctrs] = kmeans(data.epochs.data(refindex,:),ncluster);
      % sort clusters so that centers increase from 1->3
      [refctrs,ix] = sort(tmpctrs);
      refcidx = tmpidx;
      for k   = 1:ncluster,refcidx(tmpidx==k) = ix(k); end
      refcidx = refcidx(pkindex)';

      % loop over channels and record clusters at pktimes
      cluster_ndxs = zeros(length(labels),length(pkindex));
      cluster_ctrs = zeros(length(labels),ncluster);
      cluster_frac = zeros(length(labels),ncluster);
      for k  = 1:length(labels)
        this = labels{k};
        % get data for this sensor
        dat  = ts_data_selection(data,'chanlabel',this,'verbose',0);
        dat  = dat.epochs.data;
        % partition this data into ncluster=3 clusters
        [tmpidx,tmpctrs] = kmeans(dat,ncluster);
        % sort clusters so that centers increase from 1->ncluster
        [ctrs,ix] = sort(tmpctrs);
        cidx  = tmpidx;
        for j = 1:ncluster,cidx(tmpidx==j) = ix(j); end
        % store cluster indices at reference SO detection times
        cluster_cidx(k,:) = cidx(pkindex)';
        cluster_ctrs(k,:) = ctrs';
        for n = 1:ncluster, cluster_frac(k,n) = sum(cidx(pkindex)==n)/length(pkindex); end
        clear dat
      end

      figure    
      nrow  = ceil(sqrt(ncluster));
      ncol  = nrow;
      refix = strmatch(reflabel,labels,'exact');
      nchan = length(labels);
      for n = 1:ncluster
        subplot(nrow,ncol,n)
        plot(1:nchan,cluster_frac(:,n),'.',refix,cluster_frac(refix,n),'o');
        axis([0 nchan+1 0 1])
        title(sprintf('cluster %g',n))
        xlabel('channel')
        if n==1, ylabel('fraction'); end
        hline(.4,'k');
      end

      allpeaks(typ).cluster_cidx = cluster_cidx;
      allpeaks(typ).cluster_ctrs = cluster_ctrs;
      allpeaks(typ).cluster_frac = cluster_frac;
    end
  
    
    
%%
    % (Murphy et al, 03-Feb-2009. PNAS 106:5;1608-1613)
    % (Note: I've adapted their procedure to sensor space)
    % For each SO:
    % Define a 100ms window centered on sharp peak
    % Set threshold = 25% of max on this interval
    % Find maxima in other channels > threshold
    % In each channel, select maxima closest to the detection time
    % If two maxima are equally close, select the earliest
    % => each channel now has at most one maximum for this SO
    % Sort maxima (ie, channels) by the time they occurred
    % Probabilistic origin = earlest 10% of channels to have a maximum
    
    
    
    
    
    
    
    
    
    
    
    
