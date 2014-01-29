%% 08-Jun-2010
% waveforms 
  addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
  addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
  addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/scripts
  Sleep_SO_MEG_ParamFile
  if parms.clearall,clear all; Sleep_SO_MEG_ParamFile; end
  SensorEventFile         = parms.SensorEventFile;
  FlipFile                = parms.FlipFile;
  CorrSensorEventFile     = parms.CorrSensorEventFile;
  ClusterSensorEventFile  = parms.ClusterSensorEventFile;
  datestring = date;
  tstart     = tic;
  tstartorig = tstart;
  toilim     = parms.toilim; %[600 1350];
  writelog   = parms.writelog;
  layout     = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
  outpath    = [parms.rootoutdir '/' parms.SubjectID];
  FigPath    = sprintf('%s/images',outpath);
  logfile    = sprintf('%s/%s.log',outpath,date);
   
  subjects   = {'s5','s6','s7','s8'};
  subj       = strmatch(parms.SubjectID,subjects); % subj = 4
  subject    = subjects{subj};
  if strcmp(subject,'s8')
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
  
  if writelog
    fid = fopen(logfile,'a');
    fprintf(fid,'---------------------------------\n%s\nSlow oscillation analysis\n---------------------------------\n',datestr(now));
    if parms.clearall
      fprintf(fid,'Memory cleared.\n');
    end
  else
    fid = 1;
  end
  
  %% Load data (grad1 & grad2)
    loadflag = parms.loadflag; %   loadflag = 0;
    toiflag  = parms.toiflag;  %   toiflag  = 0;
    RSSflag  = parms.RSSflag;
  
  if loadflag
    % read fif files, convert data to timesurfer format, and save mat files
    matfiles = {};
    chantype = {'grad1','grad2'};
    findex   = parms.matfile_index;%1:3;
    fprintf(fid,'Processing %s data\n',[chantype{:}]);
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
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % read mat files and combine data
    fprintf(fid,'Loading MAT files:\n');
    for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
    data   = SO_combine_matfiles(matfiles(findex));
  end
  if toiflag
    if ischar(toilim) && strcmpi(toilim,'all')
      toilim = [data.epochs.time(1) data.epochs.time(end)];
    end
    fprintf(fid,'Selecting times %g-%g sec\n',toilim);
    data   = ts_data_selection(data,'toilim',toilim);
    fprintf(fid,'done.\n');
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  if RSSflag
    fprintf(fid,'Calculating RSS from grad1 & grad2\n');
    sens      = data.sensor_info;
    nchan     = length(sens);
    grad1_ind = strmatch('grad1',{data.sensor_info.typestring});
    grad2_ind = strmatch('grad2',{data.sensor_info.typestring});
    nloc      = length(grad1_ind);
    for k = 1:nloc
      Bx        = data.epochs.data(grad1_ind(k),:);
      By        = data.epochs.data(grad2_ind(k),:);
      data.epochs.data(nchan+k,:) = (Bx.^2 + By.^2).^(1/2);
      sens(nchan+k)               = sens(grad1_ind(k));
      sens(nchan+k).label         = [data.sensor_info(grad1_ind(k)).label data.sensor_info(grad2_ind(k)).label];
      sens(nchan+k).typestring    = 'RSS';
      sens(nchan+k).lognum        = min([data.sensor_info(grad1_ind(k)).lognum data.sensor_info(grad2_ind(k)).lognum]) - 1;
    end
    data.sensor_info = sens;
    data.num_sensors = length(sens);
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  %% Detection
  detection = parms.detectionflag; %   detection = 1;
  tmptype=unique({data.sensor_info.typestring}); 
  if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
  if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
  tmpstr = [tmpstr '_' datestring];  
  detectionstring = sprintf('filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s',parms.bpfreq,parms.toilim,[tmptype{:}],'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);
  % peak detection (grad1 & grad2)
  if detection
    if isempty(SensorEventFile) || ~exist(SensorEventFile,'file')
      SensorEventFile = sprintf('%s/%s_SO_init_peaks_%s.mat',outpath,subjects{subj},detectionstring);
    end
    clear tmpstr tmptype
    if exist(SensorEventFile,'file') && ~parms.overwrite
      fprintf(fid,'Loading sensor event file: %s\n',SensorEventFile);
      load(SensorEventFile); % events, peaks
      init_events = events; init_peaks = peaks; clear events peaks
    else
      fprintf(fid,'Resetting stopwatch timer for SO detection\n'); tstart = tic;
      args        = mmil_parms2args(parms);
      init_peaks  = SO_detection(data,args{:});
      t           = data.epochs.time;
      init_events      = [];
      for k   = 1:length(init_peaks)
        init_events(k).label = data.sensor_info(k).label;
        init_events(k).time  = [t([init_peaks(k).pospeak init_peaks(k).negpeak])];
        init_events(k).type  = [1*ones(1,length(init_peaks(k).pospeak)) 2*ones(1,length(init_peaks(k).negpeak))];
      end
      fprintf(fid,'saving %s\n',SensorEventFile);
      peaks  = init_peaks;
      events = init_events;
      save(SensorEventFile,'events','peaks');
      clear peaks events
    end
    telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  end
  
  bpfreq      = parms.bpfreq; if parms.bpfilter, bpfilter = 'yes'; else bpfilter = 'no'; end
  refavgflag  = parms.refavgflag;           %   refavgflag     = 1;
  plotflag    = parms.plotflag;             %   plotflag       = 1;
%   flipflag    = parms.flipflag;             %   flipflag       = 1;
%   noflipflag  = parms.noflipflag;           %   noflipflag     = 1;
  clusterflag = parms.clusterflag;
  flippad     = parms.FlipWindow;
  clusterpad  = parms.ClusterWindow;    
  epochpad    = parms.EpochPad*1000;        % s=>ms
  corrdetflag = parms.corrdetectionflag;    % whether to re-run detection after flipping
  
  if parms.preprocflag
    if ~RSSflag
      init_dat  = ts_data_selection(data,'chantype','grad');
    else
      init_dat  = data;
    end  

    fprintf(fid,'Preprocessing raw data\n');
    init_dat    = ts_preproc(init_dat,'bpfilter',bpfilter,'bpfreq',bpfreq,'blc',parms.blc,'bandpass_detrend_flag',parms.detrend);
    Fs          = init_dat.sfreq;                  % Hz
    t           = init_dat.epochs.time;            % sec
    [sel1 sel2] = match_str({init_peaks.label},{init_dat.sensor_info.label});
    init_peaks  = init_peaks(sel1);
    nchan       = init_dat.num_sensors;
    ntime       = length(init_dat.epochs.time);
    refpeaktype = parms.peaktype;               % ref peak type to use      
    if iscell(refpeaktype), refpeaktype = refpeaktype{1}; end
    npeak       = arrayfun(@(x)length(x.(refpeaktype)),init_peaks);
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  % Get flip matrix
  if parms.calc_flip_matrix
    if isempty(FlipFile) || ~exist(FlipFile,'file')
      FlipFile = sprintf('%s/%s_SO_flip_matrix_%s.mat',outpath,subjects{subj},detectionstring);
    end
    if exist(FlipFile,'file') && ~parms.overwrite
      fprintf(fid,'Loading flip matrix: %s\n',FlipFile);
      load(FlipFile); % flip
    else
      flip = calc_flip_matrix(parms,init_dat,init_peaks);
      fprintf(fid,'Saving flip file: %s\n',FlipFile);
      save(FlipFile,'flip');
      telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    end
    % inspect flip matrix
    type = 1; tmpbin=bsxfun(@eq,flip(type).matrix,flip(type).matrix');
    fprintf(fid,'Percent symmetric flip matrix (%s): %g\n',flip(type).peaktype,sum(tmpbin(:))/numel(tmpbin));
    type = 2; tmpbin=bsxfun(@eq,flip(type).matrix,flip(type).matrix');
    fprintf(fid,'Percent symmetric flip matrix (%s): %g\n',flip(type).peaktype,sum(tmpbin(:))/numel(tmpbin));
  end
  
  % Select a reference sensor
  REF   = parms.Ref;
  if ischar(REF)
    if strcmp(REF,'max')                  % channel w/ the max # of detections
      REF = find(max(npeak)==npeak);
    elseif strcmp(REF,'all')              % all channels
      REF = 1:init_dat.num_sensors;
    end
  elseif iscell(REF)                      % cell array of channel labels
    [REF,jnk] = match_str({init_peaks.label},REF);
  elseif isnumeric(REF)
    % already indices
  else
    fprintf(fid,'reference not understood\n');
    return
  end

  for refcnt = 1:length(REF)
    ref = REF(refcnt);
    reflabel    = init_peaks(ref).label;        % channel label of reference
    refpeaks    = init_peaks(ref).(refpeaktype);% ref peaks to use for flipping
    clear flip_dat flip_peaks flip_events cluster_peaks cluster_events clusters
    fprintf(fid,'Using %ss in reference channel: %s (%g of %g)\n',refpeaktype,reflabel,refcnt,length(REF));
    fprintf('Using %ss in reference channel: %s (%g of %g)\n',refpeaktype,reflabel,refcnt,length(REF));
    % Select a reference sensor (ref) and its pos or neg peaks (refpeaktype)
    if refavgflag
  %     if length(REF) > 1, error('Only one reference is supported right now; swap REF & peaktype for-loops to fix'); end    
  %     ref         = REF;                          % index of ref in init_peaks
      fprintf(fid,'Resetting stopwatch timer for flipping based on ref-averaged slow wave (%s-%s)\n',reflabel,refpeaktype);
      tstart  = tic;
      % Flip init_dat wrt polarity of ref-averaged SO waves over refpeaktype (all grads) => flip_dat
      if exist('flip','var')
        fprintf(fid,'Preparing flip vector (%s-%s) using flip matrix\n',reflabel,refpeaktype);
        if strcmp(refpeaktype,'pospeak'), flip_type = 1; else flip_type = 2; end
        flip_ind = strmatch(reflabel,{flip(flip_type).sensor_info.label});
        flip_vec = flip(flip_type).matrix(flip_ind,:);
        flip_dat = init_dat;
        flip_dat.epochs.data = init_dat.epochs.data .* repmat(flip_vec',1,length(t));
        flip_vec(flip_vec==1)  = 0;
        flip_vec(flip_vec==-1) = 1;
      else    
        % find simultaneous detections
        % - convert reference peaks to a cell array
        cellrefpeaks = num2cell(refpeaks);  % one element per refpeak
        cellepochs   = {};                  % will have one element per channel
        % - for each chan, find the refpeaks with which it is invovled
        for k = 1:nchan
          % combine pos & neg peaks in this channel & compare all to refpeaks
          tk = []; tk = init_peaks(k).pospeak; tk = sort([tk init_peaks(k).negpeak]);
          % find detections tk in k that occur near detections in cellrefpeaks
          involved      = find(cellfun(@(y)(any((tk>=y-flippad*Fs)&(tk<=y+flippad*Fs))),cellrefpeaks));
          cellepochs{k} = refpeaks(involved); % = indices to refpeaks with which chan k is involved
        end    
        % generate averages based on simultaneous detections    
        npeak           = cellfun(@length,cellepochs);
        epochdata       = SO_epochs(init_dat,cellepochs,epochpad);
        if strcmp(parms.blc,'yes'), epochdata = ts_preproc(epochdata,'blc',parms.blc); end
        averagedata     = SO_average(epochdata,npeak);
        % create flip vector & flip_dat based on average over epochs from simultaneous detections
        flip_dat  = init_dat;
        flip_vec  = zeros(1,nchan);
        tt        = averagedata.averages.time;
        tix       = nearest(tt,-flippad):nearest(tt,flippad); % flip window
        tmpdat    = averagedata.averages.data(:,tix);
        ref0      = averagedata.averages.data(ref,nearest(tt,0)); % reference value            
        refpol    = ref0 > 0; % 1 if pos, 0 if neg, reference polarity
        for k  = 1:nchan
          % select window tk+/-fpad for this channel around {tk}ref
          tmp  = tmpdat(k,:);
          % find peaks in this window in this channel
          ix   = crossing(diff(tmp));
          % find the max peak
          amp  = max(abs(tmp(ix)));
          % threshold at 25% of the max
          thsh = .25*amp;
          % get indices to peaks above the threshold
          ix   = ix(abs(tmp(ix))>thsh);
          % find the remaining peak closest to t=0
          ix   = ix(nearest(tt(tix(ix)),0));
          x0   = averagedata.averages.data(k,tix(ix));
          xpol = x0 > 0;
          if xpol ~= refpol
            % opposite polarity => flip the channel
            flip_dat.epochs.data(k,:) = -init_dat.epochs.data(k,:);
            flip_vec(k)               = 1;
          end
          clear tmp
        end
        clear tmpdat
        tmpflip = num2cell(flip_vec);
        [flip_dat.sensor_info(1:nchan).flip] = deal(tmpflip{:});
        clear tmpflip
      end
      telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    end

    %% Get grad1 & grad2 peaks (pos & neg) from flip_dat => flip_peaks
    if corrdetflag
%       if isempty(CorrSensorEventFile) || ~exist(CorrSensorEventFile,'file')
        CorrSensorEventFile = sprintf('%s/%s_SO_flip_peaks_Ref%s-%s_%s.mat',outpath,subjects{subj},reflabel,refpeaktype,detectionstring);
%       end
      if exist('flip','var')
        fprintf(fid,'Using flip vector to prepare flip-corrected peaks & events structures\n');
        % ues flip_vec to modify peaks
        flip_peaks = init_peaks;
        [flip_peaks(1:nchan).pospeak]   = deal([]);
        [flip_peaks(1:nchan).negpeak]   = deal([]);
        [flip_peaks(flip_vec==0).pospeak] = init_peaks(flip_vec==0).pospeak; % no flip  (pos=>pos)
        [flip_peaks(flip_vec==0).negpeak] = init_peaks(flip_vec==0).negpeak; % no flip  (neg=>neg)
        [flip_peaks(flip_vec==1).pospeak] = init_peaks(flip_vec==1).negpeak; % flip     (neg=>pos)
        [flip_peaks(flip_vec==1).negpeak] = init_peaks(flip_vec==1).pospeak; % flip     (pos=>neg)
        % create flip_events from flip_peaks
        flip_events  = [];
        for k   = 1:length(flip_peaks)
          flip_events(k).label = flip_dat.sensor_info(k).label;
          flip_events(k).time  = [t([flip_peaks(k).pospeak flip_peaks(k).negpeak])];
          flip_events(k).type  = [1*ones(1,length(flip_peaks(k).pospeak)) 2*ones(1,length(flip_peaks(k).negpeak))];
        end      
        peaks = flip_peaks; events = flip_events;
        fprintf(fid,'Saving flip event file: %s\n',CorrSensorEventFile);
        save(CorrSensorEventFile,'events','peaks','flip_vec','flip');
        clear peaks events
        telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);      
      elseif exist(CorrSensorEventFile,'file') && ~parms.overwrite
        fprintf(fid,'Loading flip event file: %s\n',CorrSensorEventFile);
        load(CorrSensorEventFile); % events, peaks, flip_vec
        flip_events = events; flip_peaks = peaks; clear events peaks
      else
        fprintf(fid,'Re-running detection on flipped data\n');
        args        = mmil_parms2args(parms);
        flip_peaks  = SO_detection(flip_dat,args{:}); % grad1/grad2, pos/neg
        flip_events  = [];
        for k   = 1:length(flip_peaks)
          flip_events(k).label = flip_dat.sensor_info(k).label;
          flip_events(k).time  = [t([flip_peaks(k).pospeak flip_peaks(k).negpeak])];
          flip_events(k).type  = [1*ones(1,length(flip_peaks(k).pospeak)) 2*ones(1,length(flip_peaks(k).negpeak))];
        end
        peaks  = flip_peaks;
        events = flip_events;
        fprintf(fid,'Saving flip event file: %s\n',CorrSensorEventFile);
        save(CorrSensorEventFile,'events','peaks','flip_vec');
        clear peaks events
        telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
      end
    end

    %% Cluster detections

    % Method 1 - arbitrary reference, fixed window size

    % Method 2 - all reference channels, fixed window size
    % - create clusters for each reference channel (see Method 1)
    % - select all cluster epochs from all reference channels
    % - combine epochs if any detections are within 20/40ms of each other
    % - if a channel has multiple detections in the same epoch: select the
    %   largest amplitude peak
    % - select all cluster epochs containing >20% chans
    % result: each epoch has at most one instance of a channel
    if ~isempty(parms.Coh_NFFTfactor),NFFT=round(parms.Coh_NFFTfactor*Fs); else NFFT=[]; end
    if parms.clusterflag
      warning off
      fprintf(fid,'Resetting stopwatch timer for clustering detections based on reference: %s (%s-flipped)\n',reflabel,refpeaktype); tstart = tic;
      angles        = ts_BetweenSensorAngles(flip_dat.sensor_info);
      epochpad      = ceil(epochpad*Fs/1000); % convert ms to indices
      clusterpad    = ceil(parms.ClusterWindow*Fs);  
      if isempty(ClusterSensorEventFile)% || ~exist(ClusterSensorEventFile,'file')
        ClusterSensorEventFile = sprintf('%s/%s_SO_cluster_flip_peaks_Ref%s-%s_%s.mat',outpath,subjects{subj},reflabel,refpeaktype,detectionstring);
      end
      if exist(ClusterSensorEventFile,'file') && ~parms.overwrite
        fprintf(fid,'Loading cluster event file: %s\n',ClusterSensorEventFile);
        load(ClusterSensorEventFile); % events, peaks, flip_vec, clusters
        cluster_events  = events; 
        cluster_peaks   = peaks;
      else
        % initialize cluster_peaks
        cluster_peaks = [];
        tmplabels     = {flip_dat.sensor_info.label};
        [cluster_peaks(1:nchan).label] = deal(tmplabels{:});
        clear tmplabels clusters
        [cluster_peaks(1:nchan).sfreq]  = deal(flip_dat.sfreq);
        [cluster_peaks(1:nchan).tstart] = deal(flip_dat.epochs.time(1));
        [cluster_peaks(1:nchan).tstop]  = deal(flip_dat.epochs.time(end));  
        % Loop over peaktypes (pos & neg)
        peaktypes  = {'negpeak','pospeak'};
        cnt        = 0;
        for c = 1:length(peaktypes)
          peaktype = peaktypes{c};
          fprintf(fid,'Clustering %ss across %g channels based on flipping wrt averages of %s detections in the reference channel (%s)\n',peaktype,nchan,refpeaktype,reflabel);
          % Get refpeaks for this peaktype
          refpeaks = flip_peaks(ref).(peaktype);
          allpeaks = {flip_peaks.(peaktype)};
          % Cluster {flip_peaks.(peaktype)} using ref windows centered on refpeaks
          [tmpepochs{1:nchan}]   = deal([]); % this will hold epochs for each chan
          [cluster_peaks(1:nchan).(peaktype)] = deal([]); % this will hold peaks for each chan
          if c == 1
            clusters(c)             = rmfield(flip_dat,'epochs');
            clusters(c).epochs      = [];
            clusters(c).delaydist_corrcoef   = [];  % one R per epoch per peaktype
            clusters(c).delaydist_corrcoef_p = [];  % one p per epoch per peaktype
            clusters(c).CohPSD_corrcoef      = [];  % one R per peaktype
            clusters(c).CohPSD_corrcoef_p    = [];  % one p per peaktype
            clusters(c).Coh_mean    = [];            % one mean Coh per epoch
            clusters(c).Coh_std     = [];            % one std Coh per epoch
            clusters(c).Coh_foilim  = parms.Coh_foilim;
            clusters(c).PSD_mean    = [];
            clusters(c).PSD_std     = [];
            clusters(c).PSD_foilim  = parms.PSD_foilim;
            clusters(c).N           = [];            % one N per epoch (# involved chans)
            clusters(c).reflabel    = reflabel;
            clusters(c).peaktype    = peaktype;
            clusters(c).time        = t;
            clusters(c).tstart      = flip_dat.epochs.time(1);
            clusters(c).tstop       = flip_dat.epochs.time(end);      
            clusters(c).matfiles    = data.epochs.matfiles;
            clusters(c).duration    = data.epochs.duration;
          else
            clusters(c)           = clusters(1);
            clusters(c).peaktype  = peaktype;
            clusters(c).epochs    = [];
            clusters(c).delaydist_corrcoef   = zeros(1,length(refpeaks));
            clusters(c).delaydist_corrcoef_p = zeros(1,length(refpeaks));
            clusters(c).Coh_mean  = zeros(1,length(refpeaks));
            clusters(c).Coh_std   = zeros(1,length(refpeaks));
            clusters(c).PSD_mean  = zeros(1,length(refpeaks));
            clusters(c).PSD_std   = zeros(1,length(refpeaks));            
            clusters(c).N         = zeros(1,length(refpeaks));
          end
          % loop over refpeaks
          for pk    = 1:length(refpeaks)
            clear tmpCoh tmpPSD
            [tmpCoh{1:nchan}]      = deal([]);
            [tmpPSD{1:nchan}]      = deal([]);            
            tmpCoh  = [];
            tmpPSD  = [];
            refpeak = refpeaks(pk);       % index into time vector for this refpeak
            cnt     = cnt + 1;            % cumulative detection count (=event_code)
            clusters(c).epochs(pk).event_code   = cnt;
            clusters(c).epochs(pk).refpeak      = refpeak;
            clusters(c).epochs(pk).channels     = [];
            clusters(c).epochs(pk).labels       = {};
            clusters(c).epochs(pk).peaks        = [];
            clusters(c).epochs(pk).times        = [];
            % find all detections in each refwindow refpeak+/-clusterpad
            detected = cellfun(@(x)(x(x>=refpeak-clusterpad & x<=refpeak+clusterpad)),allpeaks,'UniformOutput',false);
            epochtix = refpeak + [-epochpad:epochpad];
            cohtix   = refpeak + [-clusterpad:clusterpad];
            epochix  = epochtix>=1 & epochtix<=length(t);
            epochtix = epochtix(epochix);
            refepoch = flip_dat.epochs.data(ref,epochtix);            
            for ch   = 1:nchan
              if isempty(detected{ch}), continue; end
              % record the channels detected in this refwindow
              clusters(c).epochs(pk).channels   = [clusters(c).epochs(pk).channels ch];
              clusters(c).epochs(pk).labels     = {clusters(c).epochs(pk).labels{:} flip_dat.sensor_info(ch).label};
              % if a chan has multiple detections in the refwindow, keep the max based on flip_dat
              maxpeak = max(flip_dat.epochs.data(ch,detected{ch}))==flip_dat.epochs.data(ch,detected{ch});
              clusters(c).epochs(pk).peaks      = [clusters(c).epochs(pk).peaks detected{ch}(maxpeak)];
                % note each channel in each ref window now has at most 1 peak       
              tmpepochs{ch}(:,end+1)        = zeros(length([-epochpad:epochpad]),1);
              tmpepochs{ch}(epochix,end)    = flip_dat.epochs.data(ch,epochtix)';%cat(2,tmpepochs{ch},flip_dat.epochs.data(ch,epochtix)');
              cluster_peaks(ch).(peaktype)  = [cluster_peaks(ch).(peaktype) detected{ch}(maxpeak)];
              [Cxy,freq]  = mscohere(data.epochs.data(ref,cohtix)',data.epochs.data(ch,cohtix)',[],parms.Coh_Noverlap,NFFT,Fs);
              freq        = freq(~isnan(Cxy));
              Cxy         = Cxy(~isnan(Cxy));
              fix         = freq >= parms.Coh_foilim(1) & freq <= parms.Coh_foilim(2);
              if isfinite(Cxy(1))
                tmpCoh = [tmpCoh Cxy(fix)']; 
                if clusters(c).Coh_foilim(1) > min(freq(fix)), clusters(c).Coh_foilim(1) = min(freq(fix)); end
                if clusters(c).Coh_foilim(2) > max(freq(fix)), clusters(c).Coh_foilim(2) = max(freq(fix)); end                
              end
              clear Cxy freq fix       
            end % end loop over channels
            % sort each cluster by delay wrt the earliest detection in that refwindow
            % (channels, labels, peaks, times, delays, distances)
            [clusters(c).epochs(pk).peaks II]   = sort(clusters(c).epochs(pk).peaks);
            clusters(c).epochs(pk).times        = t(clusters(c).epochs(pk).peaks);
            clusters(c).epochs(pk).labels       = clusters(c).epochs(pk).labels(II);
            clusters(c).epochs(pk).channels     = clusters(c).epochs(pk).channels(II);
            % calc delay wrt earliest detection in refwindow
            T = clusters(c).epochs(pk).times - clusters(c).epochs(pk).times(1);
            clusters(c).epochs(pk).delays       = T;
            % calc distances bw each channel & chan w/ earliest detection in refwindow
            D = angles(clusters(c).epochs(pk).channels(1),clusters(c).epochs(pk).channels);
            D(isnan(D)) = 0;
            clusters(c).epochs(pk).distances    = D;
            % calc correlation coefficient R(delay,distance) => one value per refpeak
            [R,P] = corrcoef([D',T']);
            if size(R,2) > 1, R = R(1,2); P = P(1,2); end
            clusters(c).N(pk)     = length(D);
            if ~isnan(R)
              clusters(c).delaydist_corrcoef(pk)   = R; 
              clusters(c).delaydist_corrcoef_p(pk) = P; 
            end
            clear R P
            % calc PSD for this epoch
            [Pxx,freq]  = periodogram(refepoch',[],[],Fs);
            fix         = freq >= parms.PSD_foilim(1) & freq <= parms.PSD_foilim(2);
            clusters(c).epochs(pk).PSD      = Pxx(fix);
            clusters(c).epochs(pk).PSDfreq  = freq(fix);
            % calc mean PSD & coherence for this epoch within some foilim
            clusters(c).PSD_mean(pk)        = mean(Pxx(fix));
            clusters(c).PSD_std(pk)         = std(Pxx(fix));
            clusters(c).Coh_mean(pk)        = mean(tmpCoh);
            clusters(c).Coh_std(pk)         = std(tmpCoh);
            clear tmpCoh Pxx freq fix
          end % end loop over peaks of this peaktype
          % calc correlation coefficient R(coh,psd) over all epochs
          [R,P] = corrcoef([clusters(c).Coh_mean',clusters(c).PSD_mean']);
          if size(R,2) > 1, R = R(1,2); P = P(1,2); end
          if ~isnan(R)
            clusters(c).CohPSD_corrcoef   = R; 
            clusters(c).CohPSD_corrcoef_p = P; 
          end
          clear R P
          % epochdata & avgdata for this peaktype
          epochtime             = [-epochpad:epochpad]/Fs;
          clusters(c).epochtime = epochtime;
          clusters(c).epochdata = tmpepochs;
          ntmp                  = length(clusters(c).epochdata);
          clear tmpepochs
          clusters(c).avgdata   = zeros(ntmp,length(epochtime));
%           cohtix                = nearest(clusters(1).epochtime,0) + (-clusterpad:clusterpad);
          for ch = 1:ntmp
            if ~isempty(clusters(c).epochdata{ch})
              clusters(c).avgdata(ch,:) = mean(clusters(c).epochdata{ch},2); % ch: [time x peaks] => [time]
              % calculate coherence for individual channels across all epochs
              % - concatenate all epochs for one channel
              % - set NFFT to < length of an epoch
              % - set noverlap to 0 so that each epoch is analyzed separately
%               tmp = clusters(c).epochdata{ch}(cohtix,:); tmp = tmp(:);
%               mscohere(flip_dat.epochs.data(ref,cohtix)',flip_dat.epochs.data(ch,cohtix)',[],[],[],Fs);
%               clusters(
            end
          end
  %         telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
          fprintf(fid,'Found %g %s clusters.\n',length(clusters(c).epochs),clusters(c).peaktype);
          fprintf(fid,'Average %s cluster size = %g of %g channels.\n',peaktype,round(mean(clusters(c).N)),clusters(c).num_sensors);
          fprintf(fid,'Each channel is involved in %g %s clusters on average.\n',round(mean(cellfun(@length,{cluster_peaks.(peaktype)}))),peaktype);
          fprintf(fid,'Mean coherence (%g-%gHz) between %s & all channels across all %s epochs: %g\n',parms.Coh_foilim,reflabel,peaktype,mean(clusters(c).Coh_mean));
          fprintf(fid,'Mean PSD (%g-%gHz) in %s across all %s epochs: %g\n',parms.PSD_foilim,reflabel,peaktype,mean(clusters(c).PSD_mean));
          fprintf(fid,'[[ Correlation coefficient between epoch mean coherence and PSD: %g (p=%g) ]]\n',clusters(c).CohPSD_corrcoef,clusters(c).CohPSD_corrcoef_p);
          if parms.writelog,fprintf('Correlation coefficient between epoch mean coherence and PSD: %g (p=%g)\n',clusters(c).CohPSD_corrcoef,clusters(c).CohPSD_corrcoef_p);end
          if parms.plotflag
            tmp=flip_dat;tmp.averages=tmp.epochs;tmp=rmfield(tmp,'epochs');tmp.averages.data=clusters(c).avgdata;tmp.averages.time=clusters(c).epochtime;
            ts_ezplot(tmp,'layout',layout,'title',sprintf('Cluster averages: ref=%s (%s)',reflabel,peaktype),'showlabels','yes');
            clear tmp
          end
          clusters(c).parms = parms;
        end % end loop over peaktypes

        % create cluster_events
        cluster_events  = [];
        for k   = 1:length(cluster_peaks)
          cluster_events(k).label = cluster_peaks(k).label;
          cluster_events(k).time  = [t([cluster_peaks(k).pospeak cluster_peaks(k).negpeak])];
          cluster_events(k).type  = [1*ones(1,length(cluster_peaks(k).pospeak)) 2*ones(1,length(cluster_peaks(k).negpeak))];
        end      
        % NOTE: we now have one peaks struct for epoching all grads based on clusters defined wrt pos & neg refpeaks
        %	save cluster_peaks, cluster_events, clusters, flip_vec
        if parms.saveflag
          peaks  = cluster_peaks;
          events = cluster_events;          
          fprintf(fid,'Saving cluster event file: %s\n',ClusterSensorEventFile);
          save(ClusterSensorEventFile,'events','peaks','flip_vec','clusters');
          clear peaks events tmptype tmpstr
        end
      end
      % ------------------------------------------------------------
      % Plot some results
      C   = 55;           % head circumference, [cm]
      % note: need brain circumference! not head.
      CF  = (C/100)/360;  % conversion factor, (deg/s) x CF = (m/s)
%       CF  = CF/.3;        % sulcus correction (70% cortical surface in sulci)
%       fprintf(fid,'Assume spherical head with circumference = %gcm with 1/.3 sulcus correction.\n',C);    
      fprintf(fid,'Assume spherical head with circumference = %gcm.\n',C);
      Rthresh = [.4 .6];   
      for rr = 1:length(Rthresh)
        Rth = Rthresh(rr);            % R threshold
        Nth = round(.05*clusters(1).num_sensors); % N threhsold
        for c = 1:length(clusters)
    %       c  = 1;     % pktype
          R  = clusters(c).delaydist_corrcoef;
          N  = clusters(c).N;
          pk = find(abs(R)>Rth & N>Nth);    
          rmu=[]; dmu=[]; dmax=[]; vv=[];
      %     hfig = figure;
          for k = 1:length(pk)
            r  = clusters(c).epochs(pk(k)).distances; rmu(k)=mean(r);
            d  = clusters(c).epochs(pk(k)).delays;    dmu(k)=mean(d); dmax(k)=max(d);
            p  = polyfit(r,d,1); m = p(1); % [s/deg]
            v  = (1/m)*CF; vv(k) = v; % [m/s]
      %       figure(hfig); plot(r,d,'.'); lsline; hline(0,'k');
      %       title(sprintf('%s: %g of %g (R=%g, |v|=%g)',clusters(c).peaktype,k,length(pk),R(pk(k)),v));
      %       axis([0 130 0 1]); xlabel('angular distance'); ylabel('delay (s)');
      %       pause
          end
          fprintf(fid,'For [%s] %s clusters (n=%g of %g) with R(delay,distance) > %g (p<%g) & nchan > %g:\n',reflabel,clusters(c).peaktype,length(pk),length(clusters(c).epochs),Rth,max(clusters(c).delaydist_corrcoef_p(pk)),Nth);
          fprintf(fid,'Mean avg distance (deg): %g\n',mean(rmu));
          fprintf(fid,'Mean avg delay (s): %g\n',mean(dmu));
          fprintf(fid,'Mean max delay (s): %g\n',mean(dmax));
          fprintf(fid,'Mean propagation rate (m/s): %g\n',mean(vv));
        end
      end
  %     figure; plot(R,'.'); hline(0,'k'); hline(Rth,'r'); axis tight
  %     avgdat = flip_dat;
  %     for c  = 1:2
  %       avgdat.epochs(c)      = flip_dat.epochs;
  %       avgdat.epochs(c).data = clusters(c).avgdata;
  %       avgdat.epochs(c).time = clusters(c).epochtime;
  %       avgdat.epochs(c).event_code = c;
  %     end
  %     evcode = [2]; visualizer(ts_data_selection(avgdat,'events',evcode));
      %   arrayfun(@(x)(sum(x.type(x.type==1))),cluster_events)
      %   arrayfun(@(x)(sum(x.type(x.type==2))),cluster_events)

      % ------------------------------------------------------------
      
      % Store clusters if clusters from all refs will be combined later
      %  ...
      
      telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
      warning on
    end

  % Epoch waves
  % 	epoch flip_dat based on pos & neg cluster_peaks & then average given channel-dependent num_trials 
  % 		=> flip_[pos,neg]cluster_[epochs,averages]_waves (SlowWaves)

  % TF analysis & epoching TF zscores & TF waves

  % Plots: multiplot, topoplot {waves, zscores}

  % Event-related synchrony
    ClusterSensorEventFile = '';
    clusters = rmfield(clusters,'epochdata');
    SOreference(refcnt).reflabel    = reflabel;
    SOreference(refcnt).refpeaktype = refpeaktype;
    SOreference(refcnt).clusters    = clusters;
    SOreference(refcnt).events      = cluster_events;
    SOreference(refcnt).peaks       = cluster_peaks;
    SOreference(refcnt).flip_vec    = flip_vec;
    
    telapsed = toc(tstartorig); 
    fprintf(fid,'Time elapsed since beginning of run: %g min\n',telapsed/60);
    if parms.writelog,fprintf('Time elapsed since beginning of run: %g min\n',telapsed/60); end
    
  end % end loop over references
  
%   AllClusterSensorEventFile = sprintf('%s/%s_SO_cluster_results_Ref-all-%s_%s.mat',outpath,subjects{subj},refpeaktype,detectionstring);
%   if refcnt > 1 && (~exist(AllClusterSensorEventFile,'file') || parms.overwrite)
%     fprintf(fid,'Saving multi-reference (n=%g) cluster data file: %s\n',refcnt,AllClusterSensorEventFile);
%     save(AllClusterSensorEventFile,'SOreference');
%   end

  % This is where clusters from all reference channels can be combined
  % ...
  
  
  
  

% % 	% TF analysis
% % 	wavelet analysis (flip_dat) => flip_TFR
% % 	calc zscore 								=> flip_TFz
% % 	epoch flip_TFz based on pos & neg cluster_peaks & then average given channel-dependent num_trials 
% % 		=> flip_[pos,neg]cluster_[epochs,averages]_TFzscore (TFzscores)
% % 		note: maybe do TF rejection before averaging
% % 	average over freqband 
% % 		=> flip_[pos,neg]cluster_[epochs,averages]_TFzscore_averaged_X-YHz (TFwaves)
% % 	
% % 	multiplot average SlowWaves, TFwaves, TFzscores
% % 	topoplot series (dt=.1) average SlowWaves, TFwaves == TFzscores (freqband)
% % 	
% % 	% calculate event-related synchrony measures for flip_TFR epoched with cluster_peaks given channel-dependent num_trials
% % 		note: maybe do TF rejection before calculating synchrony	
% % 		- need new function
% % 		- syncplot & syncview for TF measures
% % 		- maybe mscohere (see August's code)
% % 		- maybe look at measures for individual epochs
% % 	
% % 	% average event-related synchrony measure over freqband
% % 		=> one wave for each channel pair and peaktype
% % 	how to view?
% % 	
% % 	% calculate single-trial synchrony measures... (not implemented yet)
% % 	
% % end reference selection
% % 	
% % 	
% % - for cluster averages with flipping based on simultaneous detections: if num_trials for a given channel is < threshold, call it a badchan and do not display in multiplot
% % 	
% % 
% % % epoch types: indiv chan peaks, refchan peaks, cluster peaks
% % % flip wrt all refpeaks or simultaneous detections



%             % Cluster-based epoching & averaging
%               clear pkcell
%               [sel1 sel2] = match_str({peaks.label},{dat.sensor_info.label});
%               pks         = peaks(sel1);
%               ref         = strmatch(reflabel,{pks.label});
%               fprintf(fid,'Epoching & averaging clusters after flipping %s wrt reference %s (%s), Method 1\n',gradtype,pks(ref).label,peaktype);
%               pksix              = pks(ref).(peaktype);
%               [pkcell{1:nchan}]  = deal(pksix);
%               npeak              = cellfun(@length,pkcell);
%               epoched5           = SO_epochs(dat,pkcell,epad);
%               if strcmp(parms.blc,'yes'), epoched5 = ts_preproc(epoched5,'blc',parms.blc); end
%               averaged5          = SO_average(epoched5,npeak);
%               pkcell5            = pkcell;
% 
%               % Method 2: average channel k wrt detection in the reference that occur within
%               % some period of a detection in k
%               det = num2cell(pksix);
%               clear pkcell
%               for k = 1:nchan
%                 % find detections in k that occur near detections in ref
%                 tk = []; 
%                 if isfield(pks,'pospeak'), tk = [tk pks(k).pospeak]; end
%                 if isfield(pks,'negpeak'), tk = [tk pks(k).negpeak]; end
%                 tk = sort(tk);
%                 ix = find(cellfun(@(y)(any((tk>=y-tau*Fs)&(tk<=y+tau*Fs))),det));
%                 pkcell{k} = pksix(ix);
%               end
% 
%               fprintf(fid,'Epoching & averaging clusters after flipping %s wrt reference %s (%s), Method 1\n',gradtype,pks(ref).label,peaktype);
%               npeak             = cellfun(@length,pkcell);
%               epoched6           = SO_epochs(dat,pkcell,epad);
%               if strcmp(parms.blc,'yes'), epoched6 = ts_preproc(epoched6,'blc',parms.blc); end
%               averaged6          = SO_average(epoched6,npeak);
%               pkcell6            = pkcell;
%               if saveflag
%                 % save results for this gradtype, peaktype, & reference
%                 fprintf(fid,'Saving cluster average file: %s\n',ClusterRefAverageFile);
%                 save(ClusterRefAverageFile,'averaged5','averaged6','pkcell5','pkcell6','description','flip');
%                 telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
%               end            
% 
%             if plotflag
%               if ~isempty(averaged5)
%                 prefix      = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged5);
%                 titlestring = sprintf('Clusters: flipped, wrt all ref detections (REF = %s)',reflabel);
%                 ts_ezplot(averaged5,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
%               end
%               if ~isempty(averaged6)
%                 prefix        = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged6);
%                 titlestring = sprintf('Clusters: flipped, wrt simultaneous ref detections (REF = %s)',reflabel);
%                 ts_ezplot(averaged6,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
%               end
%             end
%           end % end if clusterflag
%         end % end loop over references
%       end % end loop over peaktypes
%     end % end loop over gradtypes
%   end % end if refavgflag

  % TF analysis
  if parms.timefreqflag
    gammaband = parms.fband;%[10 100];
    fprintf(fid,'TF analysis\n');
    fprintf(fid,'Resetting stopwatch timer\n'); tstart = tic;
    tftoilim  = parms.tftoilim;%[1075 1225];
    foi       = parms.foi;%[2:4:54 66:4:100];
    sf        = parms.sf;%6;
    dat       = ts_data_selection(data,'toilim',tftoilim,'chantype',gradtype);
    for k = 1:nchan
      if flip(k)
        dat.epochs.data(k,:) = -dat.epochs.data(k,:);
      end
    end
%     clear data
    dat       = ts_preproc(dat,'blc','yes','dsfact',parms.dsfact,'bpfilter','yes','bpfreq',[min(foi)/2 max(foi)*2.5]);
    dat       = ts_data_selection(dat,'toilim',[tftoilim(1)+5 tftoilim(end)-5]);
    telapsed  = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
    fprintf(fid,'Beginning wavelet analysis\n');
    tfdata    = ts_freqanalysis_fieldtrip(dat,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0);
    tfzscore  = ts_zscore(tfdata,'baselinetype','zscore','blcwindow',[-inf inf],'skipzero',1);
    % visualizer(dat,tfzscore);
    
    tfzgamma  = SO_freqband_average(tfzscore,gammaband);
    tfzgamma.averages(2) = tfzgamma.averages;
    tfzgamma.averages(2).event_code = 2;
    tmpdat               = ts_preproc(dat,'bpfilter','yes','bpfreq',parms.bpfreq);
    tfzgamma.averages(2).data       = tmpdat.epochs.data; clear tmpdat
    tfzgamma.epochs                 = tfzgamma.averages;
    tfzgamma                        = rmfield(tfzgamma,'averages');
    % rescale gamma for overlay
    for k = 1:nchan
      tfzgamma.epochs(1).data(k,:) = tfzgamma.epochs(1).data(k,:) / max(abs(tfzgamma.epochs(1).data(k,:))); % normalize by max
      tfzgamma.epochs(1).data(k,:) = tfzgamma.epochs(1).data(k,:) * max(tfzgamma.epochs(2).data(k,:)); % rescale by waves
    end
    telapsed  = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
    visualizer(tfzgamma);

    % epoching
    load(CorrSensorEventFile); % events, peaks
    [sel1 sel2] = match_str({peaks.label},{tfdata.sensor_info.label});
    pks         = peaks(sel1);
    
    % Indep Epoch TF
    tfzepochs = SO_timefreq_epochs(tfzscore,{pks.(peaktype)},epad);
    npeaks    = cellfun(@length,{pks.(peaktype)});
    tfzavgdat = SO_timefreq_average(tfzepochs,npeaks);   
    ts_ezplot(tfzavgdat,'zlim',[-5 5],'foilim',gammaband,'showlabels','yes','layout',layout,'save',parms.saveflag,'close',parms.closeflag,'cond_labels',condlabels(pktype),'autoscale',0);

    % Ref Epoch TF
    pksix       = pks(ref).(peaktype);
    [pkcell{1:nchan}] = deal(pksix);
    tfzrefepochs = SO_timefreq_epochs(tfzscore,pkcell,epad);
    npeaks       = cellfun(@length,pkcell);
    tfzrefavgdat = SO_timefreq_average(tfzrefepochs,npeaks);
    ts_ezplot(tfzrefavgdat,'zlim',[-2 2],'foilim',gammaband,'showlabels','yes','layout',layout,'save',parms.saveflag,'close',parms.closeflag,'cond_labels',condlabels(pktype));
    
    % Cluster epoch TF
    
    
    
    % TF epoch rejection
    
    % Event-related synchrony
    % phase-locking
    % coherence
    
    % => syncplot
    % => syncview
    
    
%     % Indep Epoch TF
%     tfepochs = SO_timefreq_epochs(tfdata,{pks.(peaktype)},epad);
%     npeaks   = cellfun(@length,{pks.(peaktype)});
%     tfavgdat = SO_timefreq_average(tfepochs,npeaks);
%     % Z-score wrt pre-0 epoch interval (blue means lead; red means lag wrt ref)
%     tfz      = ts_zscore(tfepochs,'baselinetype','relchange','blcwindow',[-.8 -.4],'skipzero',1);
%     tfzavg   = SO_timefreq_average(tfz,npeaks);
%     % Indep-Average TF multiplot
%     ts_ezplot(tfavgdat,'showlabels','yes','layout',layout,'save',0,'close',0,'cond_labels',condlabels(pktype),'autoscale',1);
%     ts_ezplot(tfzavg,  'zlim',[-1.5 1.5],'showlabels','yes','layout',layout,'save',0,'close',0,'cond_labels',condlabels(pktype));
    % Ref Epoch TF
%     pksix       = pks(ref).(peaktype);
%     [pkcell{1:nchan}] = deal(pksix);
%     tfrefepochs = SO_timefreq_epochs(tfdata,pkcell,epad);
%     npeaks      = cellfun(@length,pkcell);
%     tfrefavgdat = SO_timefreq_average(tfrefepochs,npeaks);
%     % Z-score
%     tfrefz      = ts_zscore(tfrefepochs,'baselinetype','zscore','blcwindow',[-.8 -.4],'skipzero',1);
%     tfrefzavg   = SO_timefreq_average(tfrefz,npeaks);    
    % Ref-Average TF multiplot
%     ts_ezplot(tfrefavgdat,'showlabels','yes','layout',layout,'save',0,'close',0,'cond_labels',condlabels(pktype));
%     ts_ezplot(tfrefzavg,  'showlabels','yes','layout',layout,'save',0,'close',0,'cond_labels',condlabels(pktype));
    % Freq-band average
    % Phase coherence
    telapsed  = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);

    figure
    nsmooth     = 150;
    tmp         = SO_freqband_average(tfzscore,gammaband);
    tmp.epochs  = tmp.averages;
    tmp         = rmfield(tmp,'averages');
    ix = 1:length(x);
    t  = dat.epochs.time(ix);
    x  = sum(dat.epochs.data(:,ix),1);
    y  = sum(tmp.epochs.data(:,ix),1);
    X = smooth(x,nsmooth);
    Y = smooth(max(x).*y./max(y),nsmooth)';
    subplot(4,1,1),plot(t,x,'b'); axis tight; hline(0,'k'); title(['signal (sum over flipped ' num2str(parms.bpfreq(1)) '-' num2str(parms.bpfreq(2)) 'Hz filtered channels)'])
    subplot(4,1,2),plot(t,y','r'); axis tight; hline(0,'k'); title([num2str(gammaband(1)) '-' num2str(gammaband(2)) 'Hz gamma power (sum over flipped channels)'])
    subplot(4,1,3),plot(t,X,'b',t,max(X).*Y./max(Y),'r'); axis tight; hline(0,'k');
    legend('signal','gamma power'); title(['overlay (smoothed, ' num2str(nsmooth) '-point window)'])
    subplot(4,1,4),plot(xcorr(X,max(Y).*Y./max(Y)));
    title('cross-correlation of smoothed signals')

  end
  
  telapsed = toc(tstartorig); fprintf(fid,'Time elapsed since beginning of run: %g min\n',telapsed/60);


%% 17-Jun-2010   
%   Calculate coherence between a reference and all other channels at each reference 
%   peak over the interval t=[tk tk+tau] where tk = time of the k-th peak and tau = 400ms.  
%   Average coherence over the band [10-50]+[70-100] Hz and then over channels.  
%   Perform the calculation separately for positive and negative peaks.  
%   Overall plots of positive and negative peak-locked mean coherence over peaks 
%   1:min(# pos peaks, # neg peaks).  Is there a consistent difference between the two waveforms?  
%   Test with reference MEG0143.










  