% save results here:
% /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg (s1,s2,...)

subjnum                 = 2;
prestim                 = 1; % sec
poststim                = 1; % sec
params                  = SO_params(subjnum);
params.baseline_flag    = 1;
params.bpfilter         = 'yes';
params.bpfreq           = [.2 200];
params.bandpass_low_tb  = .2;
params.bandpass_high_tb = 5;
params.dsfact           = 3;
params.matfile_index    = 1:8;

% load and preprocess data
params.chantype = 'grad'; args = mmil_parms2args(params);
grad = ts_process_data(params.datafiles(params.matfile_index),args{:});

params.chantype = 'eeg'; args = mmil_parms2args(params);
eeg = ts_process_data(params.datafiles(params.matfile_index),args{:});

% detect pos/neg peaks in all slow oscillations
params.bpfreq = [.01 4]; args = mmil_parms2args(params);
megdet = SO_detection(grad,args{:});
eegdet = SO_detection(eeg,args{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDED 21-Mar-2011:
% save MEG & EEG detections for viewing with visualizer
t       = grad.epochs.time;
events  = [];
for k   = 1:length(megdet)
  events(k).label = grad.sensor_info(k).label;
  events(k).time  = [t([megdet(k).pospeak megdet(k).negpeak])];
  events(k).type  = [1*ones(1,length(megdet(k).pospeak)) 2*ones(1,length(megdet(k).negpeak))];
end
save('SO_detections_grad.mat','events','megdet');
events  = [];
for k   = 1:length(eegdet)
  events(k).label = eeg.sensor_info(k).label;
  events(k).time  = [t([eegdet(k).pospeak eegdet(k).negpeak])];
  events(k).type  = [1*ones(1,length(eegdet(k).pospeak)) 2*ones(1,length(eegdet(k).negpeak))];
end
save('SO_detections_eeg.mat','events','eegdet');
% visualizer(grad); % load SO_detections_grad
% visualizer(eeg);  % load SO_detections_eeg
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
nmeg      = grad.num_sensors;
neeg      = eeg.num_sensors;
nchan     = nmeg + neeg;
badchans  = [];
peaktypes = {'negpeak','pospeak'};
for c = 1:2
  peaktype = peaktypes{c};

  params.method             = 'histogram';
  params.threshold          = 3;
  params.StepSize           = .01;
  params.IntegrationWindow  = .05;
  params.MinSeparation      = .05;
  params.ClusterWindow      = .4;
  params.count              = [];

  % cluster negative EEG detections
  tmpdet            = eegdet; 
  otherpeak         = setdiff(peaktypes,peaktype);
  otherpeak         = otherpeak{1};
  [tmpdet.(otherpeak)] = deal([]); 
  params.peaktype   = peaktype; 
  args              = mmil_parms2args(params);
  [csteegdet,count] = SO_cluster_detections(tmpdet,args{:});

  % cluster all MEG detections using negative EEG peaks
  params.count          = count;
  params.peaktype       = 'allpeak';
  params.ClusterWindow  = .3;
  args                  = mmil_parms2args(params);
  [cstmegdet]           = SO_cluster_detections(megdet,args{:});

  % remove channels with fewer detections than some fixed number
  MinDetections = 100;
  det = csteegdet; badeeg = find(cellfun(@length,{det.(peaktype)}) < MinDetections);
  det = cstmegdet; badmeg = find(cellfun(@length,{det.negpeak}) + cellfun(@length,{det.pospeak})<MinDetections);

  badchans = [badchans badmeg];
  badchans = [badchans badeeg+nmeg];

  % epoch each channel separately
  neegchantrl = cellfun(@length,{csteegdet.(peaktype)});
  neegtrl     = max(neegchantrl);
  nmegchantrl = cellfun(@length,{cstmegdet.negpeak}) + cellfun(@length,{cstmegdet.pospeak});
  nmegtrl     = max(nmegchantrl);
  nchantrial  = [nmegchantrl neegchantrl];
  ntrial      = max(neegtrl,nmegtrl);
  ntime       = length(-prestim:1/grad.sfreq:poststim) - 1;

  dat = zeros(nchan,ntime,ntrial);
  for k = 1:nchan
    if ismember(k,badchans), continue; end
    if k <= nmeg
      det   = cstmegdet(k);
      tmp   = ts_data_selection(grad,'channels',k,'verbose',0);
      samp  = cat(2,det.negpeak,det.pospeak);
    else                      % EEG channel
      det   = csteegdet(k-nmeg);
      tmp   = ts_data_selection(eeg,'channels',k-nmeg,'verbose',0);
      samp  = det.(peaktype);
    end
    tmp     = ts_epoch_data(tmp,'samples',samp,'prestim',prestim,'poststim',poststim);
    dat(k,:,1:nchantrial(k)) = tmp.epochs.data;
  end
  T       = tmp.epochs.time;
  tmp     = eeg.sensor_info;  
  tmpsens = grad.sensor_info; 
  
  
  tmpsens = cat(2,tmpsens,tmp);

  if c == 1
    epoch_data = grad;
    epoch_data.epochs.data = dat;
    epoch_data.epochs.time = T;
    epoch_data.sensor_info = tmpsens;
    epoch_data.num_sensors = length(tmpsens);
    % epoch_data: t-domain, sensor-specific SOs in a zero-padded matrix
  else
    epoch_data.epochs(c)              = epoch_data.epochs(1);
    epoch_data.epochs(c).data         = dat;
  end
  epoch_data.epochs(c).num_trials     = size(dat,3);
  epoch_data.epochs(c).num_chan_trials= nchantrial;
  epoch_data.epochs(c).meg_detections = cstmegdet;
  epoch_data.epochs(c).eeg_detections = csteegdet;
  epoch_data.epochs(c).count          = count;
  epoch_data.epochs(c).params         = params;
  epoch_data.epochs(c).name           = ['eeg_' peaktype];
  epoch_data.epochs(c).event_code     = c;
  clear dat tmp tmpsens samp det
  
end
% remove bad channels
badchans    = unique(badchans);
badmeg      = badchans(badchans<=nmeg);
badeeg      = badchans(badchans>nmeg) - nmeg;
epoch_data  = ts_data_selection(epoch_data,'badchans',badchans);
for c = 1:length(epoch_data.epochs)
  epoch_data.epochs(c).meg_detections(badmeg) = [];
  epoch_data.epochs(c).eeg_detections(badeeg) = [];
  epoch_data.epochs(c).num_chan_trials(badchans)  = [];
end
save('s2_epoch_data.mat','epoch_data','params','-v7.3');
% save: grad, eeg, megdet, eegdet, epoch_data

% Condition 1 = MEG & EEG epoched wrt EEG negative peaks
% Condition 2 = MEG & EEG epoched wrt EEG positive peaks

clear grad eeg

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Randomly select a subset of each channel's epochs so that each channel
% has the same number of trials
seldata = rmsubfield(epoch_data,'epochs.data');
for c   = 1:length(seldata.epochs)
  trials= seldata.epochs(c).num_chan_trials;
  nsel  = min(trials); % select the same number of trials from each channel (= minimum # in any channel)
  dat   = zeros(seldata.num_sensors,ntime,nsel);
  for k = 1:seldata.num_sensors
    trl = trials(k);
    sel = randperm(trl);        % random list of integers from 1 to the number of trials in this channel
    sel = sort(sel(1:nsel));    % select the first nsel random numbers
    dat(k,:,:)    = epoch_data.epochs(c).data(k,:,sel);
    randsel{c}{k} = sel;
  end
  seldata.epochs(c).selected_trials = randsel{c};
  seldata.epochs(c).data            = dat;
  seldata.epochs(c).num_trials      = nsel;
  clear dat
end
% seldata.epochs = rmfield(seldata.epochs,{'meg_detections','eeg_detections','count','params','num_chan_trials','selected_trials');
ts_ezplot(seldata,'layout',params.layout,'showlabels','yes');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDED 21-Mar-2011:
% Based on visual inspection of the data, select a time interval with slow
% oscillations for the time-frequency analysis:
toilim  = [1000 1300]; % ex) [1000 1300] = 300sec = 5min from 1000sec to 1300sec
seldata = ts_data_selection(seldata,'toilim',toilim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% wavelet analysis on eeg, grad1, grad2, & norm
tic
foi     = 4:2:100; % Hz, frequencies of interest
sf      = 3;       % Hz, spectral resolution
tfdata  = ts_freqanalysis_fieldtrip(ts_data_selection(seldata,'condition',1),...
                                          'foi'         ,foi,...
                                          'sf'          ,sf,...
                                          'trials_flag' ,1,...
                                          'save_flag'   ,0,...
                                          'verbose'     ,0);                                       
tmpdata = ts_freqanalysis_fieldtrip(ts_data_selection(seldata,'condition',2),...
                                          'foi'         ,foi,...
                                          'sf'          ,sf,...
                                          'trials_flag' ,1,...
                                          'save_flag'   ,0,...
                                          'verbose'     ,0);                                       
toc

tfdata.timefreq(2) = tmpdata.timefreq; 
clear tmpdata seldata

% save results
% save('s2_data.mat','grad','eeg','-v7.3');
% save('s2_detections.mat','megdet','eegdet','params');
% save('s2_epoch_data.mat','epoch_data','params','-v7.3');
% save('s2_timefreq_data.mat','timefreq_data','params','-v7.3');

% % after Out of Memory error:
% timefreq_data = tfdata; %clear tfdata; 
% save('s2_timefreq_data_pospeak','timefreq_data','-v7.3');
% % Then:
% load s2_timefreq_data_negpeak; tfdata(1) = data_out; clear data_out;
% % load s2_timefreq_data_pospeak; 
% tfdata(2) = timefreq_data; clear timefreq_data;
% tfdata(1).timefreq(2) = tfdata(2).timefreq(1);
% tfdata = tfdata(1);

timefreq_data = ts_data_selection(tfdata,'condition',1);
save('s2_timefreq_data_negpeak.mat','timefreq_data','-v7.3');
timefreq_data = ts_data_selection(tfdata,'condition',2);
save('s2_timefreq_data_pospeak.mat','timefreq_data','-v7.3');
clear timefreq_data

% z-score average gamma power for one event
evcode        = 2;                                                      % condition number (1=negpeak, 2=pospeak)
gamma_foi     = [4 100]; % [30 100], [20 100],[30 80]                   % gamma band, Hz
blcwindow     = [-.7 -.4];                                              % baseline window, sec, used for calculating z-scores
timefreq_data = ts_data_selection(tfdata,'events',evcode);              % select one condition
f             = timefreq_data.timefreq.frequencies;
gamma_fix     = find(f>=gamma_foi(1) & f<=gamma_foi(2));                % indices to gamma frequencies
timefreq_data = ts_trials2avg(timefreq_data,'fieldname','timefreq');    % average power (over trials)
timefreq_data = ts_zscore(timefreq_data,'blcwindow',blcwindow,'baselinetype','zscore'); % calculate z-score

% plot eeg
zlim = [-10 10]; % color limits
ts_ezplot(ts_data_selection(timefreq_data,'channels',[204:258 263:271]),'showlabels','yes','zlim',zlim);
% plot meg (gradiometers)
zlim = [-5 5]; % color limits
ts_ezplot(ts_data_selection(timefreq_data,'chantype','grad'),'showlabels','yes','layout',params.layout,'zlim',zlim);

% average over gamma frequencies
tmpdata = rmfield(timefreq_data,'timefreq');
tmpdata.epochs = rmfield(timefreq_data.timefreq,{'cmplx','power','frequencies'});
tmpdata.epochs.data = mean(timefreq_data.timefreq.power(:,:,gamma_fix),3);
ts_ezplot(ts_data_selection(tmpdata,'channels',[204:258 263:271]),'showlabels','yes');
ts_ezplot(ts_data_selection(tmpdata,'chantype','grad'),'showlabels','yes','layout',params.layout);

% % artifact rejection based on average over gamma frequencies
% [grad_outdata,badchannels,badtrials] = ts_visual_reject(ts_data_selection(tmpdata,'chantype','grad'),'method','summary','metric','absmax');
% ts_ezplot(grad_outdata,'layout',params.layout,'showlabels','yes');
% 
% % remove bad trials from TF trial spectra and recalc z-score average
% trials = setdiff(1:timefreq_data.timefreq.num_trials,badtrials{1});
% timefreq_data = ts_data_selection(tfdata,'events',evcode,'trials',trials,'chantype','grad');
% timefreq_data = ts_data_selection(timefreq_data,'badchans',badchannels);
% timefreq_data = ts_trials2avg(timefreq_data,'fieldname','timefreq');
% timefreq_data = ts_zscore(timefreq_data,'blcwindow',blcwindow,'baselinetype','zscore');
% ts_ezplot(timefreq_data,'layout',params.layout,'showlabels','yes','zlim',zlim);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear tfdata timefreq_data

% calculate norm
planar    = rmsubfield(epoch_data,'epochs.data');
sens      = planar.sensor_info;
nchan     = length(sens);
% Find grad1-grad2 pairs (do not include sensors at a location unless both
% grads are present
grad1_ind = strmatch('grad1',{sens.typestring});
grad2_ind = strmatch('grad2',{sens.typestring});
badg1 = find(~ismember(grad1_ind+1,grad2_ind));
badg2 = find(~ismember(grad2_ind-1,grad1_ind));
grad1_ind(badg1) = [];
grad2_ind(badg2) = [];
nloc        = length(grad1_ind);
planar_time = planar.epochs(1).time;
for c       = 1:length(planar.epochs)
  tmpdat    = zeros(nloc,size(planar_time,2), planar.epochs(c).num_trials);
  nchantrl  = epoch_data.epochs(c).num_chan_trials;
  for kk = 1:nloc      
    idx          = [grad1_ind(kk) grad2_ind(kk)];
    nn           = min(nchantrl(grad1_ind(kk)),nchantrl(grad2_ind(kk)));
    Bx           = epoch_data.epochs(c).data(idx(1),:,1:nn);
    By           = epoch_data.epochs(c).data(idx(2),:,1:nn);
    tmpdat(kk,:,1:nn)      = (Bx.^2 + By.^2).^(1/2);
    tmp_nchantrial(kk)     = nn;
    tmp_detections{kk}     = planar.epochs(c).meg_detections(idx);
    tmpsens(kk)            =  sens(grad1_ind(kk));
    tmpsens(kk).label      = sprintf('(%s,%s)',sens(grad1_ind(kk)).label,sens(grad2_ind(kk)).label);
    tmpsens(kk).typestring = '2-norm';
    tmpsens(kk).lognum     = min([sens(grad1_ind(kk)).lognum sens(grad2_ind(kk)).lognum]) - 1;
  %   planar.epochs.meg_detections
  end
  planar.epochs(c).num_chan_trials    = tmp_nchantrial;
  planar.epochs(c).meg_detections = tmp_detections;
  planar.epochs(c).data           = tmpdat;
  clear tmpdat Bx By tmp_nchantrial tmp_detections
end
planar.sensor_info = tmpsens;
planar.num_sensors = length(tmpsens);
planar.epochs      = rmfield(planar.epochs,'eeg_detections');
clear tmpsens

% Randomly select a subset of each channel's epochs so that each channel
% has the same number of trials
seldata = rmsubfield(planar,'epochs.data');
for c   = 1:length(seldata.epochs)
  trials= seldata.epochs(c).num_chan_trials;
  nsel  = min(trials); % select the same number of trials from each channel (= minimum # in any channel)
  dat   = zeros(seldata.num_sensors,ntime,nsel);
  for k = 1:seldata.num_sensors
    trl = trials(k);
    sel = randperm(trl);        % random list of integers from 1 to the number of trials in this channel
    sel = sort(sel(1:nsel));    % select the first nsel random numbers
    dat(k,:,:)    = planar.epochs(c).data(k,:,sel);
    randsel{c}{k} = sel;
  end
  seldata.epochs(c).selected_trials = randsel{c};
  seldata.epochs(c).data            = dat;
  seldata.epochs(c).num_trials      = nsel;
  clear dat
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ADDED 21-Mar-2011:
% Based on visual inspection of the data, select a time interval with slow
% oscillations for the time-frequency analysis:
toilim  = [1000 1300]; % ex) [1000 1300] = 300sec = 5min from 1000sec to 1300sec
seldata = ts_data_selection(seldata,'toilim',toilim);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

planar = seldata;
clear seldata

% time-frequency analysis
% wavelet analysis on eeg, grad1, grad2, & norm
tic
foi     = 4:3:100; % Hz, frequencies of interest
sf      = 5;       % Hz, spectral resolution
tfnorm  = ts_freqanalysis_fieldtrip(ts_data_selection(planar,'condition',1),...
                                          'foi'         ,foi,...
                                          'sf'          ,sf,...
                                          'trials_flag' ,1,...
                                          'save_flag'   ,0,...
                                          'verbose'     ,0);                                       
tmpdata = ts_freqanalysis_fieldtrip(ts_data_selection(planar,'condition',2),...
                                          'foi'         ,foi,...
                                          'sf'          ,sf,...
                                          'trials_flag' ,1,...
                                          'save_flag'   ,0,...
                                          'verbose'     ,0); 
toc

tfnorm.timefreq(2) = tmpdata.timefreq; 
clear tmpdata

timefreq_data = ts_data_selection(tfnorm,'condition',1);
save('s2_norm_timefreq_data_negpeak.mat','timefreq_data','-v7.3');
timefreq_data = ts_data_selection(tfnorm,'condition',2);
save('s2_norm_timefreq_data_pospeak.mat','timefreq_data','-v7.3');
clear timefreq_data

% z-score average gamma power for one event
evcode        = 2;                                                      % condition number (1=negpeak, 2=pospeak)
gamma_foi     = [4 100]; % [30 100], [20 100],[30 80]                   % gamma band, Hz
blcwindow     = [-.7 -.4];                                              % baseline window, sec, used for calculating z-scores
timefreq_data = ts_data_selection(tfnorm,'events',evcode);              % select one condition
f             = timefreq_data.timefreq.frequencies;
gamma_fix     = find(f>=gamma_foi(1) & f<=gamma_foi(2));                % indices to gamma frequencies
timefreq_data = ts_trials2avg(timefreq_data,'fieldname','timefreq');    % average power (over trials)
timefreq_data = ts_zscore(timefreq_data,'blcwindow',blcwindow,'baselinetype','zscore'); % calculate z-score

% plot eeg
zlim = [-10 10]; % color limits
ts_ezplot(ts_data_selection(timefreq_data,'channels',[204:258 263:271]),'showlabels','yes','zlim',zlim);
% plot meg (gradiometers)
zlim = [-5 5]; % color limits
ts_ezplot(ts_data_selection(timefreq_data,'chantype','grad'),'showlabels','yes','layout',params.layout,'zlim',zlim);

% average over gamma frequencies
tmpdata = rmfield(timefreq_data,'timefreq');
tmpdata.epochs = rmfield(timefreq_data.timefreq,{'cmplx','power','frequencies'});
tmpdata.epochs.data = mean(timefreq_data.timefreq.power(:,:,gamma_fix),3);
ts_ezplot(ts_data_selection(tmpdata,'channels',[204:258 263:271]),'showlabels','yes');
ts_ezplot(ts_data_selection(tmpdata,'chantype','grad'),'showlabels','yes','layout',params.layout);
clear tmpdata
