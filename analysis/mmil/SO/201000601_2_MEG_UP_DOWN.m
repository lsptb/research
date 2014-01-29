tic
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

% SO detection parameters
parms     = [];
parms.fc1 = .1;   % Hz (lower cut-off freq prior to detection)
parms.fc2 = 2;    % Hz (upper cut-off freq prior to detection)
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

  %% LOAD DATA AND FIND PEAKS
  grads = {'grad1','grad2'}; graphmarker = 'o+';
  gradtype=1;%for gradtype = 1:length(grads)
    fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    data = SO_combine_matfiles(matfiles);
    toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8
    data    = ts_data_selection(data,'toilim',toilim);
    fprintf('Time = %g to %g sec\n',toilim);

    proc = ts_preproc(data,'bpfilter','yes','bpfreq',[.1 8]);
    visualizer(proc)
    
    %% time-frequency analysis
    tic
    alpha   = 1;
    bpfreq  = [5 75];    
    toi1    = [1050 1175];
    toi2    = [1050 1175];
    labels  = {'MEG0113','MEG0122','MEG0132','MEG0143','MEG0213','MEG0222','MEG0313','MEG0322','MEG0343','MEG0413','MEG0422','MEG0432','MEG0443','MEG0723'};
    % constant spectral resolution
    foi     = 10:4:70;
    sf      = 2;
%     % variable spectral resolution
%     foi     = [2:2:100 105:5:300];
%     sf      = [2*ones(1,50) 5*ones(1,40)];    
    % select data to display with TF results
    dat = ts_data_selection(data,'chanlabel',labels,'toilim',toi2);
    % filter data for display
    dat = ts_preproc(dat,'bpfilter','yes','bpfreq',[.1 10]);
    % select channels of interest
    proc    = ts_data_selection(data,'chanlabel',labels);
    % filter out undesired freqs before wavelet
    proc    = ts_preproc(proc,'bpfilter','yes','bpfreq',bpfreq);
    % remove data padding used for filtering
    proc    = ts_data_selection(proc,'toilim',toi1);
    % Morlet wavelet analysis
    tfdata  = ts_freqanalysis_fieldtrip(proc,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0);
    % remove edge artifacts produced by wavelet analysis
    tfdata  = ts_data_selection(tfdata,'toilim',toi2);
    % calculate z-scores wrt entire recording
%     zdata   = ts_zscore(tfdata,'baselinetype','relchange','blcwindow','all','verbose',0);
    % correct for 1/f^2 scaling of power spectrum
    zdata = tfdata;
    freqs = tfdata.timefreq.frequencies;
    fcorr = repmat(permute(freqs.^(alpha),[1 3 2]),[zdata.num_sensors length(zdata.timefreq.time) 1]);
    zdata.timefreq.power = zdata.timefreq.power .* fcorr;
    clear fcorr
    % average over a frequency band
    freqband= [20 60];
    tfband  = SO_freqband_average(zdata,freqband);
    % normalize band-average and scale wrt the processed time series
    % note: this is necessary for overlaying band-avg & time-domain data
    for ch = 1:tfband.num_sensors
      % de-mean
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) - mean(tfband.averages.data(ch,:));
      % normalize
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) / max(abs(tfband.averages.data(ch,:)));
      % rescale
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) * max(abs(dat.epochs(1).data(ch,:)));
    end
    % add band-average to the time-domain data structure
    dat.epochs(2)      = dat.epochs(1);
    dat.epochs(2).data = tfband.averages.data;
    % change the event code for the tfband
    dat.epochs(2).event_code = 2;
    toc
    % visualize raw time series, band-average, and TFR power
    visualizer(dat,zdata);
    
    visualizer(ts_data_selection(data,'chanlabel',labels,'toilim',toi2),...
               ts_zscore(tfdata,'baselinetype','absolute','blcwindow','all','verbose',0));
    
    % multiplot power
    ts_ezplot(zdata,'zlim',[-3 3],'showlabels','yes');
    
    % overlay smoothed band-avgs for all channels
    tmp = tfband.averages.data;
    for k = 1:size(tmp,1)
      tmp(k,:) = smooth(tmp(k,:),10);
    end
    figure; plot(tfband.averages.time,tmp); axis tight

    % overlay raw data, band-avg, and filtered signal
    filtdat = dat;
    filtdat(2) = ts_data_selection(proc,'toilim',toi2);
    filtdat(2).epochs.event_code = 3;
    filtdat    = ts_combine_data(filtdat);
    visualizer(dat);
    
    % determine alpha for freq-corr if given TFR over long period
%     % note: this only works if tfdata is very long
%     ch   = 8;
%     fix  = 1:(length(freqs)-3);
%     x    = log10(freqs(fix));
%     y    = log10(squeeze(mean(tfdata.timefreq.power(ch,:,fix),2)));
%     P    = polyfit(x',y,1);
%     yfit = P(1)*x+P(2);
%     figure; plot(x,y,'-',x,yfit,'--')
%     alpha = -P(1)
    
%% UP/DOWN state detection (Mulkovski et al, Cerebral Cortex, 2007)

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

%%
tic
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

% SO detection parameters
parms     = [];
parms.fc1 = .1;   % Hz (lower cut-off freq prior to detection)
parms.fc2 = 2;    % Hz (upper cut-off freq prior to detection)
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

  %% LOAD DATA AND FIND PEAKS
  outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

  subjects  = {'s5','s6','s7','s8'};
  subj = 4;
  fprintf('Subject %g of [%s]\n',subj,num2str(1:length(subjects)));
  subject = subjects{subj};

  grads = {'grad1','grad2'}; graphmarker = 'o+';
    gradtype=1;%for gradtype = 1:length(grads)
    fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    grad1 = SO_combine_matfiles(matfiles);

    gradtype=2;%for gradtype = 1:length(grads)
    fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    grad2 = SO_combine_matfiles(matfiles);

    planar = grad1;
    planar.epochs.data = ((grad1.epochs.data.^2)+(grad2.epochs.data.^2)).^(1/2);
    
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
    
toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8
  planar    = ts_data_selection(planar,'toilim',toilim);
  fprintf('Time = %g to %g sec\n',toilim);
  nchan       = planar.num_sensors;
%   peaks   = SO_peaks(planar,parms);
%   t       = planar.epochs.time;
%   for k = 1:102
%     planar.sensor_info(k).label      = {grad1.sensor_info(k).label grad2.sensor_info(k).label};
%     planar.sensor_info(k).typestring = {grad1.sensor_info(k).typestring grad2.sensor_info(k).typestring};
%   end
%   events  = [];
%   for k   = 1:length(peaks)
%     events(k).label = planar.sensor_info(k).label;
%     events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
%     events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
%   end
%   outfile = sprintf('%s/s8_SOpeaks_filt%g-%gHz_toi%g-%g_planar_ref%s.mat',outpath,parms.fc1,parms.fc2,toilim,'all');%[reflabels{:}]);
%   save(outfile,'events','peaks');
  
    tic
    chans   = [1 4 5 6 10 12];
    alpha   = 1;
    bpfreq  = [5 115];    
    toi1    = [600 1350];
    toi2    = [610 1340];
    % constant spectral resolution
    foi     = 20:2:100;
    sf      = 2;
   for k = 1:102
      planar.sensor_info(k).label      = grad1.sensor_info(k).label;
      planar.sensor_info(k).typestring = grad1.sensor_info(k).typestring;
   end
    % select data to display with TF results
    dat = ts_data_selection(planar,'channels',chans,'toilim',toi2,'verbose',0);
    % filter data for display
    dat = ts_preproc(dat,'bpfilter','yes','bpfreq',[.1 10]);
    % select channels of interest
    proc    = ts_data_selection(planar,'channels',chans,'verbose',0);
    % filter out undesired freqs before wavelet
    proc    = ts_preproc(proc,'bpfilter','yes','bpfreq',bpfreq);
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
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) - mean(tfband.averages.data(ch,:));
      % normalize
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) / max(abs(tfband.averages.data(ch,:)));
      % rescale
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) * max(abs(dat.epochs(1).data(ch,:)));
    end
    % add band-average to the time-domain data structure
    dat.epochs(2)      = dat.epochs(1);
    dat.epochs(2).data = tfband.averages.data;
    % change the event code for the tfband
    dat.epochs(2).event_code = 2;
    toc
    % visualize raw time series, band-average, and TFR power
    visualizer(dat,zdata);
       
    
    
    

