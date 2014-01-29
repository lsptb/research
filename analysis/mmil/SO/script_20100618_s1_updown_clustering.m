%% 18-Jun-2010
% waveforms 
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/scripts
Sleep_SO_MEG_ParamFile_s1
if parms.clearall,clear all; Sleep_SO_MEG_ParamFile_s1; end
SensorEventFile         = parms.SensorEventFile;
FlipFile                = parms.FlipFile;
CorrSensorEventFile     = parms.CorrSensorEventFile;
ClusterSensorEventFile  = parms.ClusterSensorEventFile;
datestring = date;
tstart     = tic;
tstartorig = tstart;
toilim     = parms.toilim;
writelog   = parms.writelog;
layout     = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
outpath    = [parms.rootoutdir '/' parms.SubjectID];
FigPath    = sprintf('%s/images',outpath);

subjects   = {'s1','s2','s3','s4','s5','s6','s7','s8'};
subj       = strmatch(parms.SubjectID,subjects);
subject    = subjects{subj};
logfile    = sprintf('%s/%s_%s.log',outpath,date,subject);
if strcmp(subject,'s1')
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
  loadflag = parms.loadflag;
  toiflag  = parms.toiflag;
  RSSflag  = parms.RSSflag;

if loadflag
  % read fif files, convert data to timesurfer format, and save mat files
  matfiles = {};
  chantype = {'grad1','grad2'};
  findex   = parms.matfile_index;
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
    parms.toilim = toilim;
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
detection = parms.detectionflag;
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
    init_events = [];
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
  if parms.peakpairs_flag
    tmptype=unique({data.sensor_info.typestring}); 
    if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
    if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
    if parms.peakpairs_flag,  tmpstr = [tmpstr '_pairedpeaks'];       end  
    tmpstr = [tmpstr '_' datestring];  
    detectionstring = sprintf('filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s',parms.bpfreq,parms.toilim,[tmptype{:}],'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);
    SensorEventFile = sprintf('%s/%s_SO_init_peaks_%s.mat',outpath,subjects{subj},detectionstring);
    if exist(SensorEventFile,'file') && ~parms.overwrite
      fprintf(fid,'Loading paired peak sensor event file: %s\n',SensorEventFile);
      load(SensorEventFile); % events, peaks
      init_events = events; init_peaks = peaks; clear events peaks
    else
      init_peaks  = select_peakpairs(init_peaks,parms.peakpairs_tau);
      t           = data.epochs.time;
      init_events = [];
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
  end      
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
end
t = init_peaks(1).tstart:1/init_peaks(1).sfreq:init_peaks(1).tstop;

% select toi after SO detection
if ~isempty(parms.post_toilim) && toiflag
  data = ts_data_selection(data,'toilim',parms.post_toilim);
  parms.post_toilim = [data.epochs.time(1) data.epochs.time(end)];
  parms.toilim      = parms.post_toilim;
  toilim            = parms.toilim;
  for k = 1:length(init_peaks)
    init_peaks(k).pospeak = init_peaks(k).pospeak(t(init_peaks(k).pospeak)>=parms.toilim(1) & t(init_peaks(k).pospeak)<=parms.toilim(2));
    init_peaks(k).negpeak = init_peaks(k).negpeak(t(init_peaks(k).negpeak)>=parms.toilim(1) & t(init_peaks(k).negpeak)<=parms.toilim(2));
    init_peaks(k).tstart  = parms.toilim(1);
    init_peaks(k).tstop   = parms.toilim(2);
  end
  init_events = [];
  for k   = 1:length(init_peaks)
    init_events(k).label = data.sensor_info(k).label;
    init_events(k).time  = [t([init_peaks(k).pospeak init_peaks(k).negpeak])];
    init_events(k).type  = [1*ones(1,length(init_peaks(k).pospeak)) 2*ones(1,length(init_peaks(k).negpeak))];
  end
  t = data.epochs.time;
  tmptype=unique({data.sensor_info.typestring}); 
  if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
  if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
  if parms.peakpairs_flag,  tmpstr = [tmpstr '_pairedpeaks'];       end  
  tmpstr = [tmpstr '_' datestring];  
  detectionstring = sprintf('filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s',parms.bpfreq,parms.toilim,[tmptype{:}],'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);    
end

Fs       = data.sfreq;
epochpad = ceil(parms.EpochPad*Fs);

SensorEventPowerFile = sprintf('%s/%s_SO_init_peaks_with_power_%s.mat',outpath,subjects{subj},detectionstring);
if exist(SensorEventPowerFile,'file')
    fprintf(fid,'Loading sensor event file: %s\n',SensorEventPowerFile);
    load(SensorEventPowerFile); % events, peaks
    init_events = events; init_peaks = peaks; clear events peaks
else
  fprintf(fid,'Resetting stopwatch timer for epoch power calculation\n'); tstart = tic;
  % calculate power
  PSD_offset  = round(parms.PSD_offset*Fs);
  PSD_pad     = round(parms.PSD_pad*Fs);
  PSD_window  = PSD_offset + [-PSD_pad:PSD_pad];
  PSD_foi     = parms.PSD_foi;
      % evaluate periodogram over [refpeak + PSD_window]
  peaks = init_peaks;
  for ref = 1:length(peaks)
    fprintf('Calculating power for reference %g of %g\n',ref,length(peaks));
    reflabel  = peaks(ref).label;
    refindex  = strmatch(reflabel,{data.sensor_info.label});
    pos       = peaks(ref).pospeak;
    neg       = peaks(ref).negpeak;
    pks       = [pos neg];
    npos      = length(pos);
    nneg      = length(neg);
    npks      = length(pks);
    evcodes   = [1*ones(1,npos) 2*ones(1,nneg)];
    [peaks(ref).PSD_pospeak{1:npos}] = deal([]);
    [peaks(ref).PSD_negpeak{1:nneg}] = deal([]);
    peaks(ref).PSD_freq   = PSD_foi;
    peaks(ref).PSD_offset = PSD_offset;
    peaks(ref).PSD_pad    = PSD_pad;
    % pospeaks
    peaktype  = 'pospeak'; 
    PSDfield  = 'PSD_pospeak'; 
    npks      = npos;
    for k = 1:npks
      refpeak = peaks(ref).(peaktype)(k);
      psdtix  = refpeak + PSD_window;
      if psdtix(end) <= length(t)
        [Pxx,f] = periodogram(double(data.epochs.data(refindex,psdtix)),[],PSD_foi,Fs);
        peaks(ref).(PSDfield){k} = Pxx;
      end
    end
    % negpeaks
    peaktype  = 'negpeak'; 
    PSDfield  = 'PSD_negpeak'; 
    npks      = nneg;
    for k = 1:npks
      refpeak = peaks(ref).(peaktype)(k);
      psdtix  = refpeak + PSD_window;
      if psdtix(end) <= length(t)
        [Pxx,f] = periodogram(double(data.epochs.data(refindex,psdtix)),[],PSD_foi,Fs);
        peaks(ref).(PSDfield){k} = Pxx;
      end
    end
  end
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
end

% compare power
refs    = 1:3;
foilim  = [30 50];
figure
cnt     = 1;
for ref = refs
  subplot(length(refs),1,cnt)
  fix= peaks(ref).PSD_freq>=foilim(1) & peaks(ref).PSD_freq<=foilim(2);
  ix = ~cellfun(@isempty,peaks(ref).PSD_pospeak);
  y1 = cellfun(@(x)(mean(x(fix))),peaks(ref).PSD_pospeak(ix));
  L1 = length(y1);  
  ix = ~cellfun(@isempty,peaks(ref).PSD_negpeak);
  y2 = cellfun(@(x)(mean(x(fix))),peaks(ref).PSD_negpeak(ix));
  L2 = length(y2);
  plot(1:L1,y1,'b-',1:L2,y2,'r-'); axis tight
  if ref==1, legend('pos','neg'); end
  cnt = cnt + 1;
end


refs    = 1;
foilim  = [10 100];
figure
for ref = 1:refs
  set(gcf,'Name',peaks(ref).label);
  npks  = min(100,length(peaks(ref).pospeak));
  cnt   = 1;
  for k = 1:npks
    subplot(ceil(sqrt(npks)),ceil(sqrt(npks)),cnt)
    fix= peaks(ref).PSD_freq>=foilim(1) & peaks(ref).PSD_freq<=foilim(2);
    f  = peaks(ref).PSD_freq(fix);
    y1 = peaks(ref).PSD_pospeak{k}(fix);
    y2 = peaks(ref).PSD_negpeak{k}(fix);
    plot(f,y1,'b-',f,y2,'r-'); axis tight; vline(30,'k'); vline(50,'k');
    if k==1, legend('pos','neg'); end
    title(sprintf('%g/%g',t(peaks(ref).pospeak(k)),t(peaks(ref).negpeak(k))));
    axis off
    cnt = cnt + 1;
  end
  pause
end

  
peaks = peaks_psd_gtpeak;
peaks = peaks_psd_centered;












  