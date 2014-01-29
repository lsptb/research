%% 08-Jun-2010
% waveforms 
  addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
  addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
  addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/scripts
  Sleep_SO_MEG_ParamFile
  if parms.clearall
    clear all; Sleep_SO_MEG_ParamFile; end
  tstart    = tic;
  toilim    = parms.toilim; %[600 1350];
  writelog  = parms.writelog;
  layout    = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
  outpath   = [parms.rootoutdir '/' parms.SubjectID];
  FigPath   = sprintf('%s/images',outpath);
  logfile   = sprintf('%s/%s.log',outpath,date);
  
  subjects  = {'s5','s6','s7','s8'};
  subj      = strmatch(parms.SubjectID,subjects); % subj = 4
  subject   = subjects{subj};
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
  
  %% Detection
  detection = parms.detectionflag; %   detection = 1;
  
  % peak detection (grad1 & grad2)
  if detection
    
    if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
    if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
    tmpstr = [tmpstr '_' date];
    SensorEventFile = sprintf('%s/%s_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s.mat',outpath,subjects{subj},parms.bpfreq,parms.toilim,'grad','all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);
    if exist(SensorEventFile,'file') && ~parms.overwrite
      fprintf(fid,'Loading sensor event file: %s\n',SensorEventFile);
      load(SensorEventFile); % events, peaks
    else
      fprintf(fid,'Resetting stopwatch timer\n'); tstart = tic;
      args    = mmil_parms2args(parms);
      peaks   = SO_detection(data,args{:});
      t       = data.epochs.time;
      events  = [];
      for k   = 1:length(peaks)
        events(k).label = data.sensor_info(k).label;
        events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
        events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
      end
      fprintf(fid,'saving %s\n',SensorEventFile);
      save(SensorEventFile,'events','peaks');
    end
    telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  end

  bpfreq = parms.bpfreq;
  if parms.bpfilter, bpfilter = 'yes'; else bpfilter = 'no'; end
  
  %% Create flip vector based on the SO-locked averages based on the channel with the most detections
  refavgflag  = parms.refavgflag; %   refavgflag     = 1;
  plotflag    = parms.plotflag;   %   plotflag       = 1;
  flipflag    = parms.flipflag;   %   flipflag       = 1;
  noflipflag  = parms.noflipflag; %   noflipflag     = 1;
  clusterflag = parms.clusterflag;
  if refavgflag
    fprintf(fid,'Resetting stopwatch timer\n');
    tstart  = tic;
    if ~iscell(parms.chantype)
      parms.chantype = {parms.chantype};
    end
    if ~iscell(parms.peaktype)
      parms.peaktype = {parms.peaktype};
    end
    condlabels = parms.peaktype;
    peaks_orig = peaks;
    % LOOP over gradtypes
    for gtype  = 1:length(parms.chantype)
      gradtype = parms.chantype{gtype}; %gradtype = 'grad1';
      % LOOP over peaktypes
      for pktype = 1:length(parms.peaktype)
        peaks    = peaks_orig;
        peaktype = parms.peaktype{pktype}; %peaktype = 'negpeak';
        dat      = ts_data_selection(data,'chantype',gradtype);
        dat      =  ts_preproc(dat,'bpfilter',bpfilter,'bpfreq',bpfreq,'blc',parms.blc,'bandpass_detrend_flag',parms.detrend);
        % Find the channel with the max # of detections => Reference channel
        [sel1 sel2] = match_str({peaks.label},{dat.sensor_info.label});
        pks         = peaks(sel1);
        nchan       = dat.num_sensors;
        ntime       = length(dat.epochs.time);
        npeak       = arrayfun(@(x)length(x.(peaktype)),pks);
        % Define references
        REF   = parms.Ref;
        if ischar(REF)
          if strcmp(REF,'max') % channel w/ the max # of detections
            REF = find(max(npeak)==npeak);
          elseif strcmp(REF,'all') % all channels
            REF = 1:dat.num_sensors;
          end
        elseif iscell(REF) % cell array of channel labels
          [REF,jnk] = match_str({dat.sensor_info.label},REF);
        elseif isnumeric(REF)
          % already indices
        else
          fprintf('reference not understood\n');
          return
        end
        if length(REF) > 1, error('Only one reference is supported right now; swap REF & peaktype for-loops to fix'); end
        epad  = parms.EpochPad*1000;  % s=>ms
        Fs    = dat.sfreq;       % Hz
        t     = dat.epochs.time; % sec
        tau   = parms.CoRefAvgWindow; % sec
        fpad  = parms.FlipWindow;
        cpad  = parms.ClusterWindow;
        % LOOP over reference channels
        dat_orig = dat;
        pks_orig = pks;
        for r = 1:length(REF)
          dat = dat_orig;
          pks = pks_orig;
          ref = REF(r);
          reflabel = pks(ref).label;
          RefAverageFile = sprintf('%s/s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s.mat',outpath,bpfreq,reflabel,peaktype,toilim,gradtype);
          if exist(RefAverageFile,'file') && ~parms.overwrite
            fprintf(fid,'Loading %s, %s, Ref-%s (%g-%gsec): %s\n',gradtype,peaktype,reflabel,toilim,RefAverageFile);
            load(RefAverageFile);
            saveflag = 0;
          else
            fprintf(fid,'Processing %s, %s, Ref-%s (%g-%gsec)\n',gradtype,peaktype,reflabel,toilim);
            averaged = []; averaged2 = []; averaged3 = []; averaged4 = [];          
            pkcell   = []; pkcell2   = []; pkcell3   = []; pkcell4   = [];  
            flip     = zeros(1,nchan);
            saveflag = 1;
          end
          description.averaged  = 'noflip_allrefpeaks';
          description.averaged2 = 'noflip_involvedrefpeaks';
          description.averaged3 = 'flipped_allrefpeaks';
          description.averaged4 = 'flipped_involvedrefpeaks';      
          description.averaged5 = 'flipped_clustered_allrefpeaks';
          description.averaged6 = 'flipped_clustered_involvedrefpeaks';
          pksix    = pks(ref).(peaktype);
          
          % Without Flipping
          if noflipflag && (isempty(averaged) || isempty(averaged2))
            saveflag = 1;
            % Epoch all channel wrt detection in the reference channel
            % Method 1: average channel k wrt all detections in the reference
            fprintf(fid,'Epoching & averaging %s wrt reference %s (%s), Method 1\n',gradtype,reflabel,peaktype);
            clear pkcell
            [pkcell{1:nchan}] = deal(pksix);
            npeak             = cellfun(@length,pkcell);
            epoched           = SO_epochs(dat,pkcell,epad);
            if strcmp(parms.blc,'yes'), epoched = ts_preproc(epoched,'blc',parms.blc); end
            averaged          = SO_average(epoched,npeak);
            pkcell1           = pkcell;

            % Method 2: average channel k wrt detection in the reference that occur within
            % some period of a detection in k
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

            fprintf(fid,'Epoching & averaging %s wrt reference %s (%s), Method 2\n',gradtype,pks(ref).label,peaktype);
            npeak             = cellfun(@length,pkcell);
            epoched2          = SO_epochs(dat,pkcell,epad);
            if strcmp(parms.blc,'yes'), epoched2 = ts_preproc(epoched2,'blc',parms.blc); end
            averaged2         = SO_average(epoched2,npeak);
            pkcell2           = pkcell;
          end % end noflipflag
          
          % Flipping 
          
          if flipflag && (isempty(averaged3) || isempty(averaged4))
            saveflag = 1;
            % Flip waves with opposite polarity wrt the reference
            if isempty(averaged)
              [pkcell{1:nchan}] = deal(pksix);
              npeak             = cellfun(@length,pkcell);
              epoched           = SO_epochs(dat,pkcell,epad);
              if strcmp(parms.blc,'yes'), epoched = ts_preproc(epoched,'blc',parms.blc); end
              avgdat            = SO_average(epoched,npeak);
            else
              avgdat = averaged;
            end
            tix    = nearest(avgdat.averages.time,-fpad):nearest(avgdat.averages.time,fpad); % +/-fpad window
            tmpdat = avgdat.averages.data(:,tix);
            t0ix   = nearest(avgdat.averages.time,0);
            r0     = avgdat.averages.data(ref,t0ix); % reference value            
            rpol   = r0 > 0; % 1 if pos, 0 if neg, reference polarity
            for k  = 1:nchan
              % select window tk+/-fpad for this channel around {tk}ref
              tmp  = tmpdat(k,:);
              ix   = crossing(diff(tmp));
              amp  = max(abs(tmp(ix)));
              thsh = .25*amp;
              ix   = ix(abs(tmp(ix))>thsh);
              ix   = ix(nearest(avgdat.averages.time(tix(ix)),0));
              x0   = avgdat.averages.data(k,tix(ix));
              xpol = x0 > 0;
              if xpol ~= rpol
                % opposite polarity => flip the channel
                dat.epochs.data(k,:) = -dat.epochs.data(k,:);
                flip(k)              = 1;
              end
            end

            % Redo epoching and averaging
            clear pkcell
            fprintf(fid,'Re-Epoching & averaging after flipping %s wrt reference %s (%s), Method 1\n',gradtype,pks(ref).label,peaktype);
            pksix             = pks(ref).(peaktype);
            [pkcell{1:nchan}] = deal(pksix);
            npeak             = cellfun(@length,pkcell);
            epoched3           = SO_epochs(dat,pkcell,epad);
            if strcmp(parms.blc,'yes'), epoched3 = ts_preproc(epoched3,'blc',parms.blc); end
            averaged3          = SO_average(epoched3,npeak);
            pkcell3            = pkcell;

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

            fprintf(fid,'Re-Epoching & averaging after flipping %s wrt reference %s (%s), Method 2:\n',gradtype,pks(ref).label,peaktype);
            npeak             = cellfun(@length,pkcell);
            epoched4           = SO_epochs(dat,pkcell,epad);
            if strcmp(parms.blc,'yes'), epoched4 = ts_preproc(epoched4,'blc',parms.blc); end
            averaged4          = SO_average(epoched4,npeak);
            pkcell4            = pkcell;
          end % end flipflag
          
          if saveflag
            % save results for this gradtype, peaktype, & reference
            fprintf(fid,'Saving reference-based SO-locked average file: %s\n',RefAverageFile);
            save(RefAverageFile,'averaged','averaged2','averaged3','averaged4','pkcell','pkcell2','pkcell3','pkcell4','description','flip');
            telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
          end
          if plotflag
            if ~isempty(averaged)
              prefix      = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,reflabel,peaktype(1:3),toilim,gradtype,description.averaged);
              titlestring = sprintf('Before flipping, wrt all ref detections (REF = %s)',reflabel);
              ts_ezplot(averaged,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
            end
            if ~isempty(averaged2)
              prefix      = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged2);
              titlestring = sprintf('Before flipping, wrt simultaneous ref detections (REF = %s)',reflabel);
              ts_ezplot(averaged2,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
            end
            if ~isempty(averaged3)
              prefix      = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged3);
              titlestring = sprintf('Flipped, wrt all ref detections (REF = %s)',reflabel);
              ts_ezplot(averaged3,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
            end
            if ~isempty(averaged4)
              prefix        = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged4);
              titlestring = sprintf('Flipped, wrt simultaneous ref detections (REF = %s)',reflabel);
              ts_ezplot(averaged4,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
            end
          end
          clear epoched epoched2 epoched3 epoched4 averaged averaged2 averaged3 averaged4 npeak pksix avgdat
          
          % process (corrected) dat
          
          %% Re-Detection
          corrdetectionflag = parms.corrdetectionflag; %   detection = 1;

          % peak detection (grad1 & grad2)
          if corrdetectionflag && flipflag
            CorrSensorEventFile = sprintf('%s/%s_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec_Ref-%s-%s-flipped_%s.mat',outpath,subjects{subj},parms.bpfreq,parms.toilim,gradtype,'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,reflabel,peaktype,date);
            if exist(CorrSensorEventFile,'file') && ~parms.overwrite
              fprintf(fid,'Loading (flipped) sensor event file: %s\n',CorrSensorEventFile);
              load(CorrSensorEventFile); % events, peaks
            else
              fprintf(fid,'Resetting stopwatch timer\n'); tstart = tic;
              args    = mmil_parms2args(parms);
              peaks   = SO_detection(dat,args{:});
              t       = data.epochs.time;
              events  = [];
              for k   = 1:length(peaks)
                events(k).label = dat.sensor_info(k).label;
                events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
                events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
              end
              fprintf(fid,'saving %s\n',CorrSensorEventFile);
              save(CorrSensorEventFile,'events','peaks','flip');
            end
            telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
          end
          
          % Cluster detections (best to do after flipping)
          if clusterflag
            ClusterFile = sprintf('%s/s8_cluster_%s_%s_Ref-%s_test1.mat',outpath,gradtype,[parms.peaktype{:}],reflabel);
            if ~exist(ClusterFile,'file') || parms.overwrite
              % define tk+/-ClusterWindow (cpad) around each ref peak
              [sel1 sel2] = match_str({peaks.label},{dat.sensor_info.label});
              pks         = peaks(sel1);           
              ref         = strmatch(reflabel,{pks.label});
              cpad        = ceil(cpad*dat.sfreq);
              pksix       = pks(ref).(peaktype);
              allpeaks    = {pks.(peaktype)};
              % calc distances
              theta       = ts_BetweenSensorAngles(dat.sensor_info);
              % cluster all detection within window (require same polarity across all channels)
              if pktype == 1
                clear SOclusters
                SOclusters.sensor_info = dat.sensor_info;
                SOclusters.sfreq       = dat.sfreq;
                SOclusters.tstart      = dat.epochs.time(1);
                SOclusters.tstop       = dat.epochs.time(end);
                SOclusters.coor_trans  = dat.coor_trans;
              end
              SOclusters.peaktype      = peaktype;
              for pk = 1:length(pksix)
                this     = pksix(pk);
                SOclusters(pktype).epochs(pk).event_code  = pktype;
                SOclusters(pktype).epochs(pk).reflabel    = reflabel;
                SOclusters(pktype).epochs(pk).refpeaktype = peaktype;
                SOclusters(pktype).epochs(pk).reftime     = dat.epochs.time(this);
                SOclusters(pktype).epochs(pk).channels    = [];
                SOclusters(pktype).epochs(pk).labels      = {};
                SOclusters(pktype).epochs(pk).indices     = [];
                SOclusters(pktype).epochs(pk).times       = [];
                allchans = cellfun(@(x)(x(find(x>=this-cpad & x<=this+cpad))),allpeaks,'UniformOutput',false);
                for chan = 1:nchan
                  if isempty(allchans{chan}), continue; end
                  SOclusters(pktype).epochs(pk).channels  = [SOclusters(pktype).epochs(pk).channels chan];
                  SOclusters(pktype).epochs(pk).labels    = {SOclusters(pktype).epochs(pk).labels{:} dat.sensor_info(chan).label};
                  % if a window has multiple peaks, only keep the largest amplitude
                  maxind = find(max(dat.epochs.data(chan,allchans{chan}))==dat.epochs.data(chan,allchans{chan}));
                  SOclusters(pktype).epochs(pk).indices   = [SOclusters(pktype).epochs(pk).indices allchans{chan}(maxind)];
                  % note: if flipflag=1, then polarity matches ref peak
                  % note each channel in each ref window now has at most 1 peak                
                end
                [SOclusters(pktype).epochs(pk).indices I] = sort(SOclusters(pktype).epochs(pk).indices);
                SOclusters(pktype).epochs(pk).labels      = SOclusters(pktype).epochs(pk).labels(I);
                SOclusters(pktype).epochs(pk).channels    = SOclusters(pktype).epochs(pk).channels(I);
                SOclusters(pktype).epochs(pk).times       = t(SOclusters(pktype).epochs(pk).indices);
                % calc delays wrt earliest detection in each window
                T = SOclusters(pktype).epochs(pk).times - SOclusters(pktype).epochs(pk).times(1);                
                SOclusters(pktype).epochs(pk).delay       = T;
                % calc distances bw each channel & channel w/ the earliest detection in the ref window
                D = theta(SOclusters(pktype).epochs(pk).channels(1),SOclusters(pktype).epochs(pk).channels);
                D(isnan(D)) = 0;
                SOclusters(pktype).epochs(pk).distance    = D;                
                % calc correlation coefficient analysis b/w distance & latency
                R = corrcoef(D,T);
                if size(R,2) > 1, R = R(1,2); end  
                SOclusters(pktype).epochs(pk).corrcoef    = R;
              end
              SOclusters(pktype).Ndetections = arrayfun(@(x)length(x.channels),SOclusters(pktype).epochs);
              SOclusters(pktype).corrcoef    = [SOclusters(pktype).epochs.corrcoef];
%               figure; Nthresh = 10; plot([SOclusters.epochs([SOclusters.epochs.corrcoefN]>Nthresh).corrcoefR],'.')
              % create peaks & events for clustered detections
              if pktype == 1, peaks = []; events = []; end
              for chan = 1:nchan
                lab                     = dat.sensor_info(chan).label;
                ind = arrayfun(@(x)x.indices(strmatch(lab,x.labels)),SOclusters(pktype).epochs,'UniformOutput',false);
                I   = ~cellfun(@isempty,ind);
                if pktype == 1
                  peaks(chan).label     = lab;
                  peaks(chan).sfreq     = dat.sfreq;
                  peaks(chan).tstart    = t(1);
                  peaks(chan).tstop     = t(end);
                end
                peaks(chan).(peaktype)  = [ind{I}];
              end            
              if pktype == length(parms.peaktype)
                events  = [];
                for k   = 1:length(peaks)
                  events(k).label = peaks(k).label;
                  events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
                  events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
                end                  
                fprintf(fid,'Saving cluster event file: %s\n',ClusterFile);
                save(ClusterFile,'events','peaks','SOclusters');
              end
            else
              fprintf(fid,'Loading cluster event file: %s\n',ClusterFile);
              load(ClusterFile); % events, peaks, SOclusters
            end
            
            ClusterRefAverageFile = sprintf('%s/s8_filt%g-%gHz_Ref-%s_%s_Clustered-SO-locked_average_%g-%gsec_%s.mat',outpath,bpfreq,reflabel,peaktype,toilim,gradtype);            
            if exist(ClusterRefAverageFile,'file')
              fprintf(fid,'Loading cluster average file: %s\n',ClusterRefAverageFile);
              load(ClusterRefAverageFile); % averaged5, averaged6, description, flip
            else
              % Cluster-based epoching & averaging
              clear pkcell
              [sel1 sel2] = match_str({peaks.label},{dat.sensor_info.label});
              pks         = peaks(sel1);
              ref         = strmatch(reflabel,{pks.label});
              fprintf(fid,'Epoching & averaging clusters after flipping %s wrt reference %s (%s), Method 1\n',gradtype,pks(ref).label,peaktype);
              pksix              = pks(ref).(peaktype);
              [pkcell{1:nchan}]  = deal(pksix);
              npeak              = cellfun(@length,pkcell);
              epoched5           = SO_epochs(dat,pkcell,epad);
              if strcmp(parms.blc,'yes'), epoched5 = ts_preproc(epoched5,'blc',parms.blc); end
              averaged5          = SO_average(epoched5,npeak);
              pkcell5            = pkcell;

              % Method 2: average channel k wrt detection in the reference that occur within
              % some period of a detection in k
              det = num2cell(pksix);
              clear pkcell
              for k = 1:nchan
                % find detections in k that occur near detections in ref
                tk = []; 
                if isfield(pks,'pospeak'), tk = [tk pks(k).pospeak]; end
                if isfield(pks,'negpeak'), tk = [tk pks(k).negpeak]; end
                tk = sort(tk);
                ix = find(cellfun(@(y)(any((tk>=y-tau*Fs)&(tk<=y+tau*Fs))),det));
                pkcell{k} = pksix(ix);
              end

              fprintf(fid,'Epoching & averaging clusters after flipping %s wrt reference %s (%s), Method 1\n',gradtype,pks(ref).label,peaktype);
              npeak             = cellfun(@length,pkcell);
              epoched6           = SO_epochs(dat,pkcell,epad);
              if strcmp(parms.blc,'yes'), epoched6 = ts_preproc(epoched6,'blc',parms.blc); end
              averaged6          = SO_average(epoched6,npeak);
              pkcell6            = pkcell;
              if saveflag
                % save results for this gradtype, peaktype, & reference
                fprintf(fid,'Saving cluster average file: %s\n',ClusterRefAverageFile);
                save(ClusterRefAverageFile,'averaged5','averaged6','pkcell5','pkcell6','description','flip');
                telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
              end            
            end
            if plotflag
              if ~isempty(averaged5)
                prefix      = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged5);
                titlestring = sprintf('Clusters: flipped, wrt all ref detections (REF = %s)',reflabel);
                ts_ezplot(averaged5,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
              end
              if ~isempty(averaged6)
                prefix        = sprintf('s8_filt%g-%gHz_Ref-%s_%s_SO-locked_average_%g-%gsec_%s_%s.mat',bpfreq,pks(ref).label,peaktype(1:3),toilim,gradtype,description.averaged6);
                titlestring = sprintf('Clusters: flipped, wrt simultaneous ref detections (REF = %s)',reflabel);
                ts_ezplot(averaged6,'showlabels','yes','layout',layout,'title',titlestring,'save',parms.saveflag,'close',parms.closeflag,'logfid',fid,'outpath',FigPath,'prefix',prefix,'cond_labels',condlabels(pktype),'autoscale',parms.autoscale,'axes',parms.allaxes);
              end
            end
          end % end if clusterflag
        end % end loop over references
      end % end loop over peaktypes
    end % end loop over gradtypes
  end % end if refavgflag

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
  
  

  
  
  