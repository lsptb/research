%% 07-Jun-2010
  clear all
  tstart    = tic;
  addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
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
    '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/emptyroom_nb01_060808.fif' ...
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
    fprintf(fid,'Detecting positive and negative peaks of the slow oscillation.\n');
    fprintf(fid,'SO detection parameters:\n');
    parms                 = [];
    parms.fc1             = .3;  fprintf(fid,'fc1         = %g Hz\n',parms.fc1);   % Hz (lower cut-off freq prior to detection)
    parms.fc2             = 3;   fprintf(fid,'fc2         = %g Hz\n',parms.fc2);      % Hz (upper cut-off freq prior to detection)
    parms.monothresh      = .3;  fprintf(fid,'monothresh  = %g\n',parms.monothresh);      % increase to be stricter; rec: stop at fc1
    parms.minzero         = 250; fprintf(fid,'minzero     = %g ms\n',parms.minzero);    % ms (minimum distance between zero-crossings)
    parms.amp_stdfact     = 0;   fprintf(fid,'amp_stdfact = %g\n',parms.amp_stdfact);      % # std > or < mean to threshold

    parms.monophase_flag  = 0;   fprintf(fid,'monophase_flag  = %g\n',parms.monophase_flag);
    parms.surround_flag   = 1;   fprintf(fid,'surround_flag   = %g\n',parms.surround_flag);
    parms.interdet_flag   = 0;   fprintf(fid,'interdet_flag   = %g\n',parms.interdet_flag);
    parms.zero_flag       = 1;   fprintf(fid,'zero_flag       = %g\n',parms.zero_flag);
    parms.zeroplus_flag   = 1;   fprintf(fid,'zeroplus_flag   = %g\n',parms.zeroplus_flag);
    parms.amp_flag        = 1;   fprintf(fid,'amp_flag        = %g\n',parms.amp_flag);
    parms.TFrej_flag      = 0;   fprintf(fid,'TFrej_flag      = %g\n',parms.TFrej_flag);    
    parms.preproc_flag    = 1;   fprintf(fid,'preproc_flag    = %g\n',parms.preproc_flag);

    % detect SOs and create events structure
    fprintf(fid,'Resetting stopwatch timer\n');
    tstart  = tic;
    peaks   = SO_peaks(data,parms);
    t       = data.epochs.time;
    events  = [];
    for k   = 1:length(peaks)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
      events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
    end
    SensorEventFile = sprintf('%s/%s_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s_ampstdfact%g_minzero%gms_monothresh%gHz_%s.mat',outpath,subjects{subj},parms.fc1,parms.fc2,toilim,'grad','all',parms.amp_stdfact,parms.minzero,parms.monothresh,date);
    if exist(SensorEventFile,'file')
      fprintf(fid,'not overwriting %s\n',SensorEventFile);
      return
    else
      fprintf(fid,'saving %s\n',SensorEventFile);
      save(SensorEventFile,'events','peaks');
    end
    telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  end
  
  %% Create flip vector based on the SO-locked averages based on the channel with the most detections
  fprintf(fid,'Resetting stopwatch timer\n');
  tstart  = tic;

  gradtype = 'grad1';
  peaktype = 'negpeak';
  dat      = ts_data_selection(data,'chantype',gradtype);
  dat      =  ts_preproc(dat,'bpfilter','yes','bpfreq',[parms.fc1 parms.fc2],'blc','no','bandpass_detrend_flag',0);
  % Find the channel with the max # of detections => Reference channel
  [sel1 sel2] = match_str({peaks.label},{dat.sensor_info.label});
  pks         = peaks(sel1);
  
  nchan = dat.num_sensors;
  ntime = length(dat.epochs.time);
  npeak = arrayfun(@(x)length(x.(peaktype)),pks);
  ref   = find(max(npeak)==npeak);
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
  
  ts_ezplot(averaged,'showlabels','yes','layout',layout);
  
  % Method 2: average channel k wrt detection in the reference that occur within
  % some period of a detection in k
  tau = .3; % chan k must have detection within tj+/-tau for the j-th ref detection
  det = num2cell(pksix);
  clear pkcell
  for k = 1:nchan
    % find detections in k that occur near detections in ref
    tk = pks(k).(peaktype);
    ix = find(cellfun(@(y)(any((tk>=y-tau*Fs)&(tk<=y+tau*Fs))),det));
    pkcell{k} = pksix(ix);
  end
  
  fprintf(fid,'Epoching & averaging %s wrt reference %s (%s), Method 2:\n',gradtype,pks(ref).label,peaktype);
  npeak             = cellfun(@length,pkcell);
  epoched2           = SO_epochs(dat,pkcell,pad);
  averaged2          = SO_average(epoched2,npeak);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  
  ts_ezplot(averaged2,'showlabels','yes','layout',layout);
  
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
  
  ts_ezplot(averaged3,'showlabels','yes','layout',layout);
  
  % Method 2: average channel k wrt detection in the reference that occur within
  % some period of a detection in k
  det = num2cell(pksix);
  clear pkcell
  for k = 1:nchan
    % find detections in k that occur near detections in ref
    tk = pks(k).(peaktype);
    ix = find(cellfun(@(y)(any((tk>=y-tau*Fs)&(tk<=y+tau*Fs))),det));
    pkcell{k} = pksix(ix);
  end
  
  fprintf(fid,'Epoching & averaging %s wrt reference %s (%s), Method 2:\n',gradtype,pks(ref).label,peaktype);
  npeak             = cellfun(@length,pkcell);
  epoched4           = SO_epochs(dat,pkcell,pad);
  averaged4          = SO_average(epoched4,npeak);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  
  ts_ezplot(averaged4,'showlabels','yes','layout',layout);
  
  % TF analysis
  foi = [2:4:56 64:4:100];
  sf  = 6;
  fprintf(fid,'Wavelet analysis\n');
  tfdata   = ts_freqanalysis_fieldtrip(dat,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0);
  telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  % Epoch TF data for one reference channel
  
  % SO-locked averages of SO pos and neg-locked grad1 and grad2 TF power
  % => four multiplots over both grads for one reference channel from one
  % chantype
  
  % grand epoching to create SO epochs
  
  % Re-epoch wave and TF data using grand epochs
  
  % Next: 
  %   - delay vs distance analysis (R)
  %   - origin analysis
  %   - propagation analysis (streamlines)
  % save results
  
  % low-freq spectral analysis (FFT, TF, & TF time-average)
  % cluster analysis for origin determination
  % monte carlo stats for pos/neg differences
  % event-related and single-trial analysis of SO & gamma phase-synchrony
  % RSS planar analysis
  % source space analysis
  
  %%
    % Create grad1 flip vector
    % note 1: lo vs hi are separated by the Sylvian fissure and are located
    %         on lower and higher parts of the dewar, respectively.
    % note 2: L vs R are separated by the interhemispheric fissure.
    % 
      % IMPORTANT!!
      % THIS WILL NEED TO BE REPLACED BY A BLOCK THAT DETERMINES WHETHER OR
      % NOT TO FLIP BASED ON THE REF-BASED SO-LOCKED AVERAGES
      
    L_lochannels = [1:4 55:58 63:66 73 80 82];
    R_lochannels = [51:54 81 89 95:102];
    L_hichannels = [5:25 28 59:62 67:72 74:75 78];
    R_hichannels = [31:34 36:50 76:77 83:88 90:94];
    mid_channels = [29 30 35 22 79 80];
    all_channels = sort([L_lochannels R_lochannels L_hichannels R_hichannels mid_channels]);
    posUP_chan   = sort([L_lochannels R_hichannels]);
    negUP_chan   = sort([R_lochannels L_hichannels]);
    nchan        = length(peaks);
    
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
    ntime  = length(t);
    pad    = round(tpad*Fs);
    tbin   = zeros(nchan,ntime); % binary time vector (1=withinSO)
    grad1  = ts_data_selection(grad1,'toilim',toilim);
    
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
    GroupEventFile = sprintf('%s/%s_SO_grand_events_%gsecSOepochs_maxsumlo%g_%s_%g-%gsec%s.mat',outpath,subjects{subj},tpad,ceil(maxsum_lowerlimit),grads{gradtype},toilim,suffix);
    if exist(GroupEventFile,'file')
      fprintf('not overwriting %s\n',GroupEventFile);
      return
    else
      fprintf('saving %s\n',GroupEventFile);
      save(GroupEventFile,'events');
    end
    if 0
      % Visualize the corrected waves with grand SO events
      flipped.epochs.data(end,:)   = [SO_diffsum' 0];
      flipped.epochs.data(end-1,:) = SO_possum;
      flipped.epochs.data(end-2,:) = SO_sum;
      visualizer(flipped);
    end
   
    %% 04-Jun-2010. DELAYvsDISTANCE & ORIGIN analysis
    % Distances
    % Euclidean distance
    cfg.layout =[];
    cfg.layout = layout;
    cfg.layout = prepare_layout(cfg); % create 2d layout from 3d positions
    distance2D = dist(cfg.layout.pos(1:102,1:2)');
                  % distance(i,j) = sqrt( (x(i)-x(j)).^2 + (y(i) - y(j)).^2 );
    % Angular distance
    zshift     = .08; % shifts the Z center to almost the middle of the helmet
    theta3D    = ts_BetweenSensorAngles(flipped.sensor_info,zshift);
                  % theta(i,j) = real(acos(dot(Pi,Pj)/(norm(Pi)*norm(Pj)))*180/pi)
                  %              where Pk = (x,y,z)-coord of k-th sensor
    theta3D(isnan(theta3D)) = 0;
    
    % Sort channels by detection time for each SO epoch
    nchan  = length(events);
    target = 2; % Which peak to use for sorting. 1-pospeak, 2-negpeak.
    times  = arrayfun(@(x)(x.time(x.type==target)),events,'UniformOutput',false);
    SObeg  = t(SO_beg);
    SOend  = t(SO_end);
    nepoch = length(SObeg);
    clear slowwaves
    slowwaves.sensor_info = flipped.sensor_info;
    slowwaves.theta3D     = theta3D;
    slowwaves.distance2D  = distance2D;
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
      
      % store sorted channel labels, indices, detection times, & distances
      slowwaves.waves(s).channels   = chanix;
      slowwaves.waves(s).labels     = labels;
      slowwaves.waves(s).times      = sorted_times;
      slowwaves.waves(s).theta3D    = slowwaves.theta3D(chanix(1),chanix);
      slowwaves.waves(s).distance2D = slowwaves.distance2D(chanix(1),chanix);
      
      % calculate correlation coefficients & best fit line for delay vs dist
      dt    = slowwaves.waves(s).times - min(slowwaves.waves(s).times); % delay
      dx2d  = slowwaves.waves(s).distance2D;                            % 2D distance
      dx3d  = slowwaves.waves(s).theta3D;                               % 3D distance
      [c,r] = polyfit(dx2d,dt,1);                           % best fit line, c=[m b]=[slope,intersect]
      [R,p] = corrcoef(dx2d,dt);
      R     = R(1,2);
      slowwaves.waves(s).DelayDist2D_R_m_b = [R c(1) c(2)]; % corrcoef,slope,intersect
      [c,r] = polyfit(dx3d,dt,1);                           % best fit line, c=[m b]=[slope,intersect]
      [R,p] = corrcoef(dx3d,dt);
      R     = R(1,2);
      slowwaves.waves(s).DelayDist3D_R_m_b = [R c(1) c(2)]; % corrcoef,slope,intersect
    end
    
    if 0
      % 3D plot probable origins (earliest 10%)
      figure
      pkeep  = 10; % percent of channels to keep as probable origins for each SO
      nkeep  = ceil(pkeep*nchan/100);
      sens   = flipped.sensor_info;
      x = arrayfun(@(x)(x.loc(1,4)),sens,'UniformOutput',true)';
      y = arrayfun(@(x)(x.loc(2,4)),sens,'UniformOutput',true)';
      z = arrayfun(@(x)(x.loc(3,4)),sens,'UniformOutput',true)';
      for s = 1:length(slowwaves.waves)
        chanix  = slowwaves.waves(s).channels(1:nkeep);
        chanix  = chanix(chanix~=0);
        [s1 s2] = match_str({sens.label},{events(chanix).label});

        % 3D plot dots at sensor locations
        clf
        subplot(2,2,1),plot3(x,y,z,'.')
        title(['t = ' num2str(min(slowwaves.waves(s).times(1:nkeep)))]);
        for a = 1:length(chanix)
          % convert chan index in events to chan index in grad1 struct
          b = s1(a);
          % add channel numbers of probable location to 3D plot
          text(x(b), y(b), z(b),num2str(b))
        end; grid on; vline(0,'g'); hline(0,'r'); view(3)
        subplot(2,2,2),plot3(x,y,z,'.'),title('top')
        for a = 1:length(chanix)
          b = s1(a); text(x(b), y(b), z(b),num2str(b))
        end; grid on; vline(0,'g'); hline(0,'r'); view(2)
        subplot(2,2,3),plot3(x,y,z,'.'),title('back')
        for a = 1:length(chanix)
          b = s1(a); text(x(b), y(b), z(b),num2str(b))
        end; grid on; vline(0,'g'); hline(0,'r'); view(0,0)
        subplot(2,2,4),plot3(x,y,z,'.'),title('front')
        for a = 1:length(chanix)
          b = s1(a); text(x(b), y(b), z(b),num2str(b))
        end; grid on; vline(0,'g'); hline(0,'r'); view(180,0)      
        pause
      end
    end
    
    % Examine correlation coefficients
    R2d = arrayfun(@(x)(x.DelayDist2D_R_m_b(1)),slowwaves.waves);
    R3d = arrayfun(@(x)(x.DelayDist3D_R_m_b(1)),slowwaves.waves);
    if 0
      figure
      subplot(1,2,1),plot(R2d,'.'),axis tight,set(gca,'ylim',[-1 1]); title('R(delay,distance)');xlabel('epoch')
      subplot(1,2,2),plot(R3d,'.'),axis tight,set(gca,'ylim',[-1 1]); title('R(delay,angle)')
    end
    
    AcceptanceLevel = .5;
    GoodWaves2D     = find(R2d > AcceptanceLevel);
    GoodWaves3D     = find(R3d > AcceptanceLevel);
    GoodWaves       = intersect(GoodWaves2D,GoodWaves3D);
    GoodWaveTimes   = arrayfun(@(x)(x.times(1)),slowwaves.waves(GoodWaves));
        
    if 0
      SelWaves = [131:201]; % GoodWaves
      % sensor locations in 3D
      x = arrayfun(@(x)(x.loc(1,4)),slowwaves.sensor_info,'UniformOutput',true)';
      y = arrayfun(@(x)(x.loc(2,4)),slowwaves.sensor_info,'UniformOutput',true)';
      z = arrayfun(@(x)(x.loc(3,4)),slowwaves.sensor_info,'UniformOutput',true)';
      figure
      npp   = 5;                                 % number of waves per page
      ngood = length(SelWaves);
      npp   = min(npp,ngood);
      cnt   = 0;
      for s = 1:ngood
        cnt = cnt + 1;
        if cnt==1, clf; end
        ind = SelWaves(s);
        ref = slowwaves.waves(ind).channels(1);
        dt  = slowwaves.waves(ind).times - min(slowwaves.waves(ind).times);
        dx2d= slowwaves.waves(ind).distance2D;
        dx3d= slowwaves.waves(ind).theta3D;
        subplot(3,npp,cnt),      plot(dx2d,dt,'b.'), lsline; axis([0 100 -1 1]);
        title(sprintf('%gsec',slowwaves.waves(ind).times(1)));
        text(60,-.8,sprintf('R=%g',slowwaves.waves(ind).DelayDist2D_R_m_b(1)));
        if mod(s,npp)==1, ylabel('delay (sec)'); xlabel('2D distance'); end
        subplot(3,npp,npp+cnt),  plot(dx3d,dt,'r.'), lsline; axis([0 180 -1 1]);
        text(100,-.8,sprintf('R=%g',slowwaves.waves(ind).DelayDist3D_R_m_b(1)));
        if mod(s,npp)==1, ylabel('delay (sec)'); xlabel('3D angle');    end
        subplot(3,npp,2*npp+cnt),plot3(x,y,z,'k.'); hold on
        plot3(x(ref),y(ref),z(ref),'b*','MarkerSize',10)
        for k = 1:length(dt),b=slowwaves.waves(ind).channels(k); text(x(b), y(b), z(b),num2str(b)); end
        title(sprintf('%g = %s',ref,slowwaves.waves(ind).labels{1}))
        grid on; vline(0,'g'); hline(0,'r'); view(3); axis tight
        if mod(s,npp)==0, pause, cnt = 0; end
      end
    end
    
    SlowWaveFile = sprintf('%s/%s_SO_Waves_%s_%g-%gsec%s.mat',outpath,subjects{subj},grads{gradtype},toilim,suffix);
    slowwaves.parms.SOdetection = parms;
    slowwaves.parms.Grouping.Flipping         = 'manual';
    slowwaves.parms.Grouping.TargetPeak       = 'UPstate';
    slowwaves.parms.Grouping.Pad              = tpad;
    slowwaves.parms.Grouping.SmoothFrame      = window;
    slowwaves.parms.CorrCoef.zshift           = zshift;
    slowwaves.parms.CorrCoef.AcceptanceLevel  = AcceptanceLevel;
    slowwaves.parms.matfiles= matfiles;
    slowwaves.parms.toilim  = toilim;
    slowwaves.parms.suffix  = suffix;
    slowwaves.parms.outpath = outpath;
    slowwaves.results.SensorEventFile = SensorEventFile;
    slowwaves.results.GroupEventFile  = GroupEventFile;
    slowwaves.results.SlowWaveFile    = SlowWaveFile;
    slowwaves.results.NumWaves      = length(slowwaves.waves);
    slowwaves.results.WaveOrigins   = arrayfun(@(x)(x.labels(1)),slowwaves.waves);
    slowwaves.results.WaveTimes     = arrayfun(@(x)(x.times(1)),slowwaves.waves);
    slowwaves.results.NumGoodWaves  = length(GoodWaves);
    slowwaves.results.GoodWaveIndex = GoodWaves;
    slowwaves.filename              = SlowWaveFile;
    slowwaves.timestamp             = datestr(now);
    slowwaves.SubjID                = subjects{subj};
    
    if exist(SlowWaveFile,'file')
      fprintf('not overwriting %s\n',SlowWaveFile);
      return
    else
      fprintf('saving %s\n',SlowWaveFile);
      save(SlowWaveFile,'slowwaves','events','peaks');
    end
    toc
    
    
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
    
    
    
    
    
    
    
    
    
    
    
    
