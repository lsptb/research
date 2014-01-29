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
  badlabels  = {};
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
    for f = 1:length(fiffiles)
      fif = fiffiles{f};
      [fpath,fname,fext]  = fileparts(fif);
      outfile             = sprintf('%s/matfiles/%s_eeg.mat',outpath,fname);
      matfiles{end+1}     = outfile;
      if exist(outfile,'file') % never overwrite (param independent)
        fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
        continue
      else
        fprintf(fid,'Reading FIF file: %s\n',fif);
      end
      data = ts_loadfif(fif,{'eeg'},'epochs');
      fprintf(fid,'Saving MAT file: %s\n',outfile);
      save(outfile,'data');
      clear fif data
    end    
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % read mat files and combine data
    fprintf(fid,'Loading MAT files:\n');
    for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
    data   = SO_combine_matfiles(matfiles(findex));
    findex = length(fiffiles) + findex;
    for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
    eeg    = SO_combine_matfiles(matfiles(findex));    
  end
  if toiflag
    if ischar(toilim) && strcmpi(toilim,'all')
      toilim = [data.epochs.time(1) data.epochs.time(end)];
      parms.toilim = toilim;
    end
    fprintf(fid,'Selecting times %g-%g sec\n',toilim);
    data   = ts_data_selection(data,'toilim',toilim);
    eeg    = ts_data_selection(eeg,'toilim',toilim);
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
  
  noneeg = [61:64]; [eeg.sensor_info(noneeg).badchan] = deal(1);
  eeg = ts_data_selection(eeg); % 70 channel EEG cap
  
  
  %% Laplacian
  % EEG Coordinates
  
  % Standard EEG labels and coords
  eegcap; % eeglabels, xx, yy, ii
  [eeg.sensor_info(1:end).label] = deal(eeglabels{:});
  
  % Remove bad channels
  [s1 s2] = match_str({eeg.sensor_info.label},badlabels);
  if ~isempty(s1),eeg = ts_data_selection(eeg,'badchans',s1); end

%   % Using location & transformation matrices from fif files
%   nchan = eeg.num_sensors;
%   T     = eeg.coor_trans.device2head;
%   pos   = zeros(nchan,3); % (x,y,z) for each channel
%   for k = 1:nchan
%     loc         = eeg.sensor_info(k).loc;
%     loc         = T*loc;
%     pos(k,1:3)  = loc(1:3,4);
%   end
%   method = 'stereographic'; % gnomic, stereographic, ortographic, inverse, polar
%   prj    = elproj(pos, method); % * [0 1; -1 0];
%             % ELPROJ makes a azimuthal projection of a 3D electrode cloud
%             %  on a plane tangent to the sphere fitted through the electrodes
%             %  the projection is along the z-axis
%   X = prj(:,1);   % x-coordinates
%   Y = prj(:,2);   % y-coordinates
%   % Scale the data to a circle with x-axis and y-axis: 0.05 to 0.95
%   offset = .01/.9; % -.5 => [-.45 .45]
%   x = 0.9*((X-min(X))/(max(X)-min(X))+offset); % y
%   y = 0.9*((Y-min(Y))/(max(Y)-min(Y))+offset); % x
%  
% figure('Name',['Subject loc * device2head (' method ' projection)']);
% for k = 1:length(X)
%   subplot('position',[x(k) y(k) .04 .04]); axis off
%   title(sprintf('%s (%g)',eeglabels{k},ii(k)));
% end

% Create Laplacian nearest-neighbor matrix
EEG_neighbor_matrix = false(eeg.num_sensors,eeg.num_sensors);
thresh = 25;
theta  = ts_BetweenSensorAngles(eeg.sensor_info,0);
for k  = 1:size(theta,1)
  ch   = find(theta(k,:)<thresh);
  [tmp,I] = sort(theta(k,ch));
  ch   = ch(I);
  if length(ch) > 4
    ch = ch(1:4);
  end
  EEG_neighbor_matrix(k,ch) = true;
%   fprintf('%g: %s\n',k,num2str(ch));
end
figure; imagesc(EEG_neighbor_matrix);

% Calculate Laplacian matrix
nsamp = length(eeg.epochs.time);
lap   = zeros(nchan,nsamp);
for k = 1:nchan
  neighb = eeg.epochs.data(EEG_neighbor_matrix(k,:),:);
  tmp    = sum(neighb,1) / size(neighb,1);
  lap(k,:) = eeg.epochs.data(k,:) - tmp;
  clear tmp neighb
end

% Update eeg structure with Laplacian data
eeg.epochs.data = lap; clear lap
eeg = ts_data_selection(eeg,'badchans',find(sum(EEG_neighbor_matrix,2)<4));

% bck=data; data=eeg; save('matfiles/sleep_s1_1-6_raw_eeg_Laplacian.mat','data'); data=bck; clear bck

% SO detection
args      = mmil_parms2args(parms);
eeg_peaks = SO_detection(eeg,args{:});

init_peaks  = eeg_peaks;
t           = eeg.epochs.time;
init_events = [];
for k   = 1:length(init_peaks)
  init_events(k).label = eeg.sensor_info(k).label;
  init_events(k).time  = [t([init_peaks(k).pospeak init_peaks(k).negpeak])];
  init_events(k).type  = [1*ones(1,length(init_peaks(k).pospeak)) 2*ones(1,length(init_peaks(k).negpeak))];
end
peaks=init_peaks; events=init_events;
save('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s1/s1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_eeg-Laplacian_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat','events','peaks');
init_peaks  = select_peakpairs(init_peaks,parms.peakpairs_tau);
init_events = [];
for k   = 1:length(init_peaks)
  init_events(k).label = eeg.sensor_info(k).label;
  init_events(k).time  = [t([init_peaks(k).pospeak init_peaks(k).negpeak])];
  init_events(k).type  = [1*ones(1,length(init_peaks(k).pospeak)) 2*ones(1,length(init_peaks(k).negpeak))];
end
peaks=init_peaks; events=init_events;
save('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s1/s1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_eeg-Laplacian_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_pairedpeaks_18-Jun-2010.mat','events','peaks');
clear peaks events

  %% Detection
  detection = parms.detectionflag;
  tmptype=unique({data.sensor_info.typestring}); 
  if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
  if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
  if parms.peakpairs_flag,  tmpstr = [tmpstr '_pairedpeaks'];       end  
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
      if parms.peakpairs_flag
        init_peaks  = select_peakpairs(init_peaks,parms.peakpairs_tau);
      end      
      t             = data.epochs.time;
      init_events   = [];
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
  
  
