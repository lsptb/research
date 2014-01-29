%% 08-Jun-2010
% waveforms

  clear all
  tstart    = tic;
  addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
  addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
  layout    = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
  
  subjects  = {'s5','s6','s7','s8'};
  subj      = 4;
  subject   = subjects{subj};
  outpath   = ['/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/' subject];
  logfile   = sprintf('%s/%s.log',outpath,date);
  writelog  = 1;
  detection = 1;
  fiffiles  = {...
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_1_nb01_060808.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_2_nb01_060808.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_3_nb01_060808.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_4_nb01_060808.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_5_nb01_060808.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_6_nb01_060808.fif' ...
%     '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/emptyroom_nb01_060808.fif' ...
    };  
  
  if writelog
    fid = fopen(logfile,'a');
    fprintf(fid,'---------------------------------\n%s\nSlow oscillation analysis\n---------------------------------\n',date);
  else
    fid = 1;
  end
  % read fif files, convert data to timesurfer format, and save mat files
  matfiles = {};
  chantype = {'grad1','grad2'};
  findex   = 1:3;
  fprintf(fid,'Processing %s data\n',[chantype{:}]);
  for f = 1:length(fiffiles)
    fif = fiffiles{f};
    [fpath,fname,fext]  = fileparts(fif);
    outfile             = sprintf('%s/matfiles/%s_grad.mat',outpath,fname);
    matfiles{end+1}     = outfile;
    if exist(outfile,'file')
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
  findex = 1:3; % which files  
  for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
  data   = SO_combine_matfiles(matfiles(findex));
%   toilim = [600 1350];
  toilim = [600 1350]; %[data.epochs.time(1) data.epochs.time(end)];  
  fprintf(fid,'Selecting times %g-%g sec\n',toilim);
  data   = ts_data_selection(data,'toilim',toilim);
  fprintf(fid,'done.\n');
  telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);

  % peak detection
  if detection
    fprintf(fid,'Resetting stopwatch timer\n'); tstart = tic;    
    parms = [];
    parms.bpfilter          = 1;
    parms.decimate          = 0;
    parms.smooth            = 0;
    parms.hilbertpeaks      = 1;
    parms.derivpeaks        = 0;
    parms.zerocross         = 1;
    parms.monotonic         = 1;
    parms.gtmedian          = 1;
    parms.return_zerocross  = 0;

    parms.bpfreq            = [.2 5];   % Hz
    parms.decimate_factor   = 4;
    parms.smooth_window     = .05;      % sec
    parms.zero2zero_limits  = [.25 1];  % sec
    parms.toilim            = toilim;
    parms.mindist           = .2;       % sec
    
    args  = mmil_parms2args(parms);
    peaks = SO_detection(data,args{:});

    t       = data.epochs.time;
    events  = [];
    for k   = 1:length(peaks)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
      events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
    end
    SensorEventFile = sprintf('%s/%s_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec_%s.mat',outpath,subjects{subj},parms.bpfreq,parms.toilim,'grad','all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,date);
    if exist(SensorEventFile,'file')
      fprintf(fid,'not overwriting %s\n',SensorEventFile);
      return
    else
      fprintf(fid,'saving %s\n',SensorEventFile);
      save(SensorEventFile,'events','peaks');
    end
    telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  end

%   load new_SO_detection_algorithm_s8_600-1350sec_grad1_bp.3-5Hz_smooth.05sec_zero2zero.25-1sec_20100608.mat
  bpfreq = parms.bpfreq;
  
  %% Create flip vector based on the SO-locked averages based on the channel with the most detections
  fprintf(fid,'Resetting stopwatch timer\n');
  tstart  = tic;

  gradtype = 'grad1';
  peaktype = 'negpeak';
  dat      = ts_data_selection(data,'chantype',gradtype);
  dat      =  ts_preproc(dat,'bpfilter','yes','bpfreq',bpfreq,'blc','no','bandpass_detrend_flag',0);
  % Find the channel with the max # of detections => Reference channel
  [sel1 sel2] = match_str({peaks.label},{dat.sensor_info.label});
  pks         = peaks(sel1);
  
  nchan = dat.num_sensors;
  ntime = length(dat.epochs.time);
  npeak = arrayfun(@(x)length(x.(peaktype)),pks);
  ref   = 4;%find(max(npeak)==npeak);
  Fs    = dat.sfreq;
  t     = dat.epochs.time;
  pad   = 1000; % ms
  
  % Epoch all channel wrt detection in the reference channel
  % Method 1: average channel k wrt all detections in the reference
  fprintf(fid,'Epoching & averaging %s wrt reference %s (%s)\n',gradtype,pks(ref).label,peaktype);
  clear pkcell
  pksix             = pks(ref).(peaktype);
  [pkcell{1:nchan}] = deal(pksix);
  npeak             = cellfun(@length,pkcell);
  epoched           = SO_epochs(dat,pkcell,pad);
  averaged          = SO_average(epoched,npeak);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  
  ts_ezplot(averaged,'showlabels','yes','layout',layout,'title','Before flipping, wrt all ref detections');
  
  % Method 2: average channel k wrt detection in the reference that occur within
  % some period of a detection in k
  tau = .2; % chan k must have detection within tj+/-tau for the j-th ref detection
  det = num2cell(pksix); % reference detections
  clear pkcell
  for k = 1:nchan
    % find detections in k that occur near detections in ref
    % combine pos & neg peak detections in indiv chans before comparison
    tk = []; tk = pks(k).pospeak; tk = sort([tk pks(k).negpeak]);
    % find ref detections for which this chan has a detection within tau 
    ix = find(cellfun(@(y)(any((tk>=y-tau*Fs)&(tk<=y+tau*Fs))),det));
    % only keep the ref detections that are shared
    pkcell{k} = pksix(ix);
  end
  
  fprintf(fid,'Epoching & averaging %s wrt reference %s (%s), Method 2:\n',gradtype,pks(ref).label,peaktype);
  npeak             = cellfun(@length,pkcell);
  epoched2           = SO_epochs(dat,pkcell,pad);
  averaged2          = SO_average(epoched2,npeak);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  
  ts_ezplot(averaged2,'showlabels','yes','layout',layout,'title','Before flipping, wrt simultaneous ref detections');
  
  % Create flip vector
  
  % SO-locked averages of SO pos and neg-locked grad1 and grad2 waveforms
  % => four multiplots over both grads for one reference channel from one
  % chantype
  avgdat = averaged;
  t0ix   = nearest(epoched.epochs.time,0);
  r0     = avgdat.averages.data(ref,t0ix); % reference value
  rpol   = r0 > 0; % 1 if pos, 0 if neg, reference polarity
  for k  = 1:nchan
    % if polarity is opposite of ref, flip the channel
    x0   = avgdat.averages.data(k,t0ix);
    xpol = x0 > 0;
    if xpol ~= rpol
      % opposite polarity => flip the channel
      dat.epochs.data(k,:) = -dat.epochs.data(k,:);
    end
  end
  
  % redo epoching and averaging
  clear pkcell
 fprintf(fid,'Re-Epoching & averaging after flipping %s wrt reference %s (%s)\n',gradtype,pks(ref).label,peaktype);
  pksix             = pks(ref).(peaktype);
  [pkcell{1:nchan}] = deal(pksix);
  npeak             = cellfun(@length,pkcell);
  epoched3           = SO_epochs(dat,pkcell,pad);
  averaged3          = SO_average(epoched3,npeak);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  
  ts_ezplot(averaged3,'showlabels','yes','layout',layout,'title','Flipped, wrt all ref detections');
  
  % Method 2: average channel k wrt detection in the reference that occur within
  % some period of a detection in k
  det = num2cell(pksix);
  clear pkcell
  for k = 1:nchan
    % find detections in k that occur near detections in ref
    tk = []; tk = pks(k).pospeak; tk = sort([tk pks(k).negpeak]);
    ix = find(cellfun(@(y)(any((tk>=y-tau*Fs)&(tk<=y+tau*Fs))),det));
    pkcell{k} = pksix(ix);
  end
  
  fprintf(fid,'Epoching & averaging %s wrt reference %s (%s), Method 2:\n',gradtype,pks(ref).label,peaktype);
  npeak             = cellfun(@length,pkcell);
  epoched4           = SO_epochs(dat,pkcell,pad);
  averaged4          = SO_average(epoched4,npeak);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  
  ts_ezplot(averaged4,'showlabels','yes','layout',layout,'title','Flipped, wrt simultaneous ref detections');
  
try save('s8_600-1350sec_grad1_neg_ref4.mat','averaged','averaged2','averaged3','averaged4'); end
% TFR


