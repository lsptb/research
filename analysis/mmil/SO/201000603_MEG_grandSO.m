
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

%% 03-Jun-2010 - grad2
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
       
%% 03-Jun-2010. Grand SO
  addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
  outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

  subjects  = {'s5','s6','s7','s8'};
  subj = 4;
  fprintf('Subject %g of [%s]\n',subj,num2str(1:length(subjects)));
  subject = subjects{subj};

  grads = {'grad1','grad2'}; graphmarker = 'o+';
    gradtype = 1;
    fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    grad1 = SO_combine_matfiles(matfiles);
    tic
    Fs      = peaks(1).sfreq;
    tstart  = peaks(1).tstart;
    tstop   = peaks(1).tstop;
    toilim  = [tstart tstop]; 
    grad1   = ts_data_selection(grad1,'toilim',toilim);
    t       = tstart:1/Fs:tstop;
    
    load /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8/s8_SOpeaks_filt0.3-3Hz_toi600-1350_grad1_refall.mat
    polarity = {'pospeak','negpeak'};
%     minsum   = 2;
    tpad     = .25; % seconds, epoch = [tk-tpad:dt:tk+tpad]
    combevnt = 1; % whether to combine grand SO edge events with indiv detections
    labels = {peaks.label};
    nchan  = length(peaks);
    ntime  = length(t);
    pad    = round(tpad*Fs);
    tbin   = zeros(nchan,ntime);
    % Loop over channels to create individual SO epochs & look across
    % channels to define the grand SO epoch
    for ch = 1:nchan
      tind  = [];
      for k = 1:length(polarity), tind = [tind [peaks(ch).(polarity{k})]]; end
      c   = num2cell(tind);
      tmp = cellfun(@(x)[x-pad:x+pad],c,'UniformOutput',false);
      tmp = [tmp{:}];
      tmp = unique(tmp);
      tbin(ch,tmp) = 1;
      clear tmp
    end
    tsosum = sum(tbin,1);
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % optional constraint on the minimum # of overlapping channels
%     tsosum(tsosum<minsum) = 0;
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%
    tsobin = tsosum > 15;
    dtsobin= diff(tsobin);
    so_beg = find(dtsobin > 0);
    so_end = find(dtsobin < 0);
    rmix   = (t(so_end)-t(so_beg))<.5;
    so_end(rmix) = [];
    so_beg(rmix) = [];
    sotime = [t(so_beg) t(so_end)];
    if combevnt
      soedge = [1*ones(1,length(so_beg)) 2*ones(1,length(so_end))] + max(events(1).type);
      for k  = 1:length(events)
        events(k).time = [events(k).time sotime];
        events(k).type = [events(k).type soedge];
      end
      save(sprintf('s8_SO_grand_combo-events_%gsecSOepochs_grad1_600-1350sec_%s.mat',tpad,[polarity{:}]),'events');
%       save(sprintf('s8_SO_grand_combo-events_%gsecSOepochs_grad1_600-1350sec_%s_sumgt%g.mat',tpad,[polarity{:}],minsum),'events');
    else
      [tmp(1:nchan).label] = deal(labels{:});
      soedge = [1*ones(length(so_beg)) 2*ones(length(so_end))];
      [tmp(1:nchan).time]  = deal(sotime);
      [tmp(1:nchan).type]  = deal(soedge);
      events = tmp;
%       save(sprintf('s8_SO_grand_events_%gsecSOepochs_grad1_600-1350sec_%s.mat',tpad,[polarity{:}]),'events');
    end

%% 03-Jun-2010. Manual flipping and grand SO.

%   % detect SOs and create events structure
%   peaks   = SO_peaks(grad1,parms);
%   t       = grad1.epochs.time;
%   events  = [];
%   for k   = 1:length(peaks)
%     events(k).label = grad1.sensor_info(k).label;
%     events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
%     events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
%   end
%   outfile = sprintf('%s/s8_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s.mat',outpath,parms.fc1,parms.fc2,toilim,grads{gradtype},'all');
%   save(outfile,'events','peaks');

  % load events structure
  load /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8/s8_SOpeaks_filt0.3-3Hz_toi600-1350_grad1_refall.mat
  
    % Create grad1 flip vector
    % note 1: lo vs hi are separated by the Sylvian fissure and are located
    %         on lower and higher parts of the dewar, respectively.
    % note 2: L vs R are separated by the interhemispheric fissure.
    L_lochannels = [1:4 55:58 63:66 73 80 82];
    R_lochannels = [51:54 81 89 95:102];
    L_hichannels = [5:25 28 59:62 67:72 74:75 78];
    R_hichannels = [31:34 36:50 76:77 83:88 90:94];
    mid_channels = [29 30 35 22 79 80];
    all_channels = sort([L_lochannels R_lochannels L_hichannels R_hichannels mid_channels]);
    posUP_chan   = sort([L_lochannels R_hichannels]);
    negUP_chan   = sort([R_lochannels L_hichannels]);
        
    tpad              = .05; % seconds, epoch = [tk-tpad:dt:tk+tpad]
    maxsum_lowerlimit = .1*nchan;
    nchan_lowerlimit  = .05*nchan;
    window            = .1; % seconds, smoothing window
    
    Fs     = peaks(1).sfreq;
    tstart = peaks(1).tstart;
    tstop  = peaks(1).tstop;
    toilim = [tstart tstop]; 
    t      = tstart:1/Fs:tstop;    
    labels = {peaks.label};
    nchan  = length(peaks);
    ntime  = length(t);
    pad    = round(tpad*Fs);
    tbin   = zeros(nchan,ntime); % binary time vector (1=withinSO)
    
    % Initialize corrected-wave structure
    flipped = grad1;
    
    % UP-state locked analysis. Force everything to be neg-UP.
    % loop over channels
    for n = 1:nchan
      y = grad1.epochs.data(n,:);
      % if pos-UP
      if ismember(n,posUP_chan)
        y                        = -y; % flip
        flipped.epochs.data(n,:) =  y;
%         events(n).time           = [events(n).time(events(n).type==2) events(n).time(events(n).type==1)];        
        ix1 = events(n).type==1;
        ix2 = events(n).type==2;
        events(n).type(ix1) = 2;
        events(n).type(ix2) = 1;
        peaktype = 'pospeak';
      else % if neg-UP or midline channel
        peaktype = 'negpeak';
      end
      % create binary time vector based on the selected peaktype
      c   = num2cell(peaks(n).(peaktype));
      tmp = cellfun(@(x)[x-pad:x+pad],c,'UniformOutput',false);
      tmp = [tmp{:}];
      tmp = unique(tmp);
      tbin(n,tmp) = 1;
      clear tmp y
      % find start of each epoch
      chanevnt(n).label = peaks(n).label;      
      begtimes = t(diff(tbin(n,:)) > 0);
      endtimes = t(diff(tbin(n,:)) < 0);
      if isempty(begtimes) || isempty(endtimes), continue; end;
      cntr     = (endtimes-begtimes)/2;
      chanevnt(n).times = [begtimes' endtimes' cntr' t(peaks(n).(peaktype))']; % (# detections) x 4      
      chanevnt(n).index = peaks(n).(peaktype);
    end
    % Combine SO epochs across channels to create the grand SO
    SO_sum = sum(tbin,1);
      % at any time, this is the # of channels that have SO epochs containing the point
      % note: SO_sum goes up with padding
    % Initial Smooth SO_sum (100ms => 100-point)
    npoint    = round(window*Fs);
    SO_sum    = smooth(SO_sum,npoint,'lowess');
    SO_possum = SO_sum > 0;
    SO_diffsum = diff(SO_possum);
    SO_beg     = find(SO_diffsum > 0);
    SO_end     = find(SO_diffsum < 0);
    SOepochs   = cellfun(@(x,y)[x,y],num2cell(SO_beg),num2cell(SO_end),'UniformOutput',false);

    % reprocess channels
    for n = 1:nchan
      % Apply constraints on sums
      % (1) repetition constraint (at most one time per polarity per chan per SO epoch)
      tmpindex            = cellfun(@(x)find(chanevnt(n).times(:,4)>t(x(1)) & chanevnt(n).times(:,4)<t(x(2))),SOepochs,'UniformOutput',false);
      tmpindex            = tmpindex(~cellfun(@isempty,tmpindex));
      if isempty(tmpindex)
        chanevnt(n).index = []; 
      else
        tmpindex          = cellfun(@(x)x(1),tmpindex);
        chanevnt(n).index = chanevnt(n).index(tmpindex);    
        chanevnt(n).times = chanevnt(n).times(tmpindex,:); 
      end
    end
    
    % (2) maxsum threshold: max(SO_sum) > threshold
    % note: bigger threshold => less diffuse detection times across chans
    for n = 1:length(SOepochs)
      idx = SOepochs{n}(1):SOepochs{n}(2);
      t1  = t(SOepochs{n}(1));
      t2  = t(SOepochs{n}(2));
      tmp = arrayfun(@(x)find((x.times(:,4)>=t1)&(x.times(:,4)<=t2)),chanevnt,'UniformOutput',false);
      tmp = cellfun(@length,tmp);
      if sum(tmp) < nchan_lowerlimit
        tbin(:,idx) = 0;
      elseif max(SO_sum(idx)) < maxsum_lowerlimit
        tbin(:,idx) = 0;
      end
    end
      
    % (3) distribution constraint
    % split bi- or multi-modal SO_sum distributions at local minima

    % Find start and stop of each grand SO
    SO_sum     = sum(tbin,1);
    SO_sum     = smooth(SO_sum,npoint,'lowess');
    SO_possum  = SO_sum > 0;
    SO_diffsum = diff(SO_possum);
    SO_beg     = find(SO_diffsum > 0);
    SO_end     = find(SO_diffsum < 0);
    
    % Create events structure for marking start & stop of grand SO
    SO_times   = [t(SO_beg) t(SO_end)];
    SO_types   = [1*ones(1,length(SO_beg)) 2*ones(1,length(SO_end))]+2;

%     [tmp(1:nchan).label] = deal(labels{:});
%     [tmp(1:nchan).time]  = deal(SO_times);
%     [tmp(1:nchan).type]  = deal(SO_types);
%     events = tmp;

    for k = 1:length(events)
      events(k).time = [events(k).time SO_times];
      events(k).type = [events(k).type SO_types+2];
    end
    save(sprintf('s8_SO_grand_events_%gsecSOepochs_maxsumlo%g_grad1_600-1350sec_%s.mat',tpad,ceil(maxsum_lowerlimit),'flipped-negUP'),'events');
    
    % Visualize the corrected waves with grand SO events
    flipped.epochs.data(end,:)   = [SO_diffsum' 0];
    flipped.epochs.data(end-1,:) = SO_possum;
    flipped.epochs.data(end-2,:) = SO_sum;
    visualizer(flipped);
    
    %% ORIGIN analysis
    
    % loop over SO epochs
    pkeep  = 10; % percent of channels to keep as probable origins for each SO
    nchan  = length(events);
    nkeep  = ceil(pkeep*nchan/100);
    target = 2; % 1-pospeak, 2-negpeak
    times  = arrayfun(@(x)(x.time(x.type==target)),events,'UniformOutput',false);
    SObeg  = t(SO_beg);
    SOend  = t(SO_end);
    nepoch = length(SObeg);
    tmp    = zeros(1,nkeep);
    clear origin
    for s = 1:nepoch
      lo = SObeg(s);
      hi = SOend(s);
      ix = cellfun(@(x)(find(x>=lo&x<=hi)),times,'UniformOutput',false);
      chanix  = find(~cellfun(@isempty,ix));
      tmptime = times(chanix);% detection times for channels in this epoch
      ix = ix(chanix);        % events for channels in this epoch
      % TMPTIME HAS TOO MANY POINTS IN SOME CASES
      % NEED TO FIGURE OUT WHAT IS GOING WRONG
      % TEMP FIX: set x(y) to x(y(1))
      % NEED TO UNDO IN THE FUTURE
      tmptime = cellfun(@(x,y)(x(y(1))),tmptime,ix); 

      labels  = {events(chanix).label};
      % NOTE: channel labels{ch} has peak type=target during epoch at tmptime{ch}
      % sort detection times
      [sorted_times, sorted_idx] = sort(tmptime);
      % rearrange labels to match sorted_times
      labels = labels(sorted_idx);
      chanix = chanix(sorted_idx);
      % store the earliest pkeep%
      origin(s).channels(1:min(length(labels),nkeep)) = chanix(1:nkeep);
      origin(s).labels(1:min(length(labels),nkeep)) = labels(1:nkeep);
      origin(s).times(1:min(length(labels),nkeep)) = sorted_times(1:nkeep);
    end
    clear tmp

    % 3D plot probable origins (earliest 10%)
    figure
    sens  = flipped.sensor_info;
    x = arrayfun(@(x)(x.loc(1,4)),sens,'UniformOutput',true)';
    y = arrayfun(@(x)(x.loc(2,4)),sens,'UniformOutput',true)';
    z = arrayfun(@(x)(x.loc(3,4)),sens,'UniformOutput',true)';
    for s = 1:length(origin)
      chanix  = origin(s).channels;
      chanix  = chanix(chanix~=0);
      [s1 s2] = match_str({sens.label},{events(chanix).label});

      % 3D plot dots at sensor locations
      clf
      subplot(2,2,1),plot3(x,y,z,'.')
      title(['t = ' num2str(min(origin(s).times))]);
      for a = 1:length(chanix)
        % convert chan index in events to chan index in grad1 struct
        b = s1(a);
        % add channel numbers of probable location to 3D plot
        text(x(b), y(b), z(b),num2str(b))
      end; grid on; vline(0,'g'); hline(0,'r'); view(3)
      subplot(2,2,2),plot3(x,y,z,'.')
      for a = 1:length(chanix)
        b = s1(a); text(x(b), y(b), z(b),num2str(b))
      end; grid on; vline(0,'g'); hline(0,'r'); view(2)
      subplot(2,2,3),plot3(x,y,z,'.')
      for a = 1:length(chanix)
        b = s1(a); text(x(b), y(b), z(b),num2str(b))
      end; grid on; vline(0,'g'); hline(0,'r'); view(0,0)
      subplot(2,2,4),plot3(x,y,z,'.')
      for a = 1:length(chanix)
        b = s1(a); text(x(b), y(b), z(b),num2str(b))
      end; grid on; vline(0,'g'); hline(0,'r'); view(180,0)      
      pause
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
    
    
    
    
    
    
    
    
    
    
    
    
