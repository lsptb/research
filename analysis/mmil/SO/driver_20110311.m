% params.chantype = 'eeg';      %'grad';%'grad';%'eeg';      % chantype must be 'eeg' or 'grad'
% params.peaktype = 'pospeak';  %'allpeak';%'pospeak';%'negpeak';  % peaktype to use for aggregate count (pospeak,negpeak,bothpeaks/allpeaks)
% params.cluster_count = [];    %'eeg_negpeak';%'eeg_pospeak';
% params.toilim   = [];
% params.cluster_ClusterWindow  = .4; % cluster window size, sec; .4 eeg, .2 grad
% params.MinChansPerCluster     = 30; % 30 M/EEG, 20 iEEG
% params.layout = [];                 % = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I. Detect slow oscillations in EEG and MEG

subjnum                   = 2;
cluster_thresh            = 3;
IntervalSelection_thresh  = 3;

params                    = SO_params(subjnum);
params.cluster_thresh     = cluster_thresh;
params.IntervalSelection_CountThreshold = IntervalSelection_thresh;

% EEG
params.chantype               = 'eeg';
params.cluster_ClusterWindow  = .4;
params.MinChansPerCluster     = 20;

params.peaktype               = 'negpeak';
SO_analyze_20110311(subjnum,params);

params.peaktype               = 'pospeak';
SO_analyze_20110311(subjnum,params);

% MEG
params.chantype               = 'grad';
params.cluster_ClusterWindow  = .3;
params.MinChansPerCluster     = 30;

% Reference: EEG (negative peak)
params.cluster_count          = 'eeg_negpeak';
params.peaktype               = 'allpeak';
SO_analyze_20110311(subjnum,params);

% Reference: EEG (positive peak)
params.cluster_count          = 'eeg_pospeak';
params.peaktype               = 'allpeak';
SO_analyze_20110311(subjnum,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% II. Epoch the EEG & MEG data based on SO detected in EEG (use sufficient padding for wavelet analysis)

prestim   = 1; % sec
poststim  = 1;
blcwindow = [-1 -.75];

cd(params.SubjDir);

% 1. load EEG & MEG data and detection clusters
load matfiles/proc_grad_epoch_data_ICA.mat;                 megdata   = epoch_data; clear epoch_data
load matfiles/proc_eeg_epoch_data_ICA.mat;                  eegdata   = epoch_data; clear epoch_data
load matfiles/SO_clustered_eeg_negpeak_detections_sws.mat;  eegnegdet = detections;
load matfiles/SO_clusters_eeg_negpeak_noflip.mat;           eegnegcst = clusters_Near_pks(2);
load matfiles/SO_clustered_eeg_pospeak_detections_sws.mat;  eegposdet = detections;
load matfiles/SO_clusters_eeg_pospeak_noflip.mat;           eegposcst = clusters_Near_pks(1);

clear eegepoch megepoch planarepoch
for i = 1:2
  if i==1, eegcst = eegnegcst; eegdet = eegnegdet; end
  if i==2, eegcst = eegposcst; eegdet = eegposdet; end

  % 2. epoch EEG data
  Fs       = eegcst.sfreq;
  tstart   = eegcst.abs_toilim(1);
  
  % epoch using peaks from detection histogram
  trigtime = [eegcst.epochs.HistTime];
  trigsamp = round((trigtime-tstart)*Fs);

  % % epoch using earliest detection in each oscillation
  % trigtime = [eegcst.epochs.RefTime];
  % trigsamp = round((trigtime-tstart)*Fs);
  
  % % epoch using detections in a reference channel
  % refchan  = 'EEG 001'; refind = strmatch(refchan,{eegdet.label});
  % trigsamp = eegdet(refind).negpeak;
  
  eegepoch(i) = ts_epoch_data(eegdata,'samples',trigsamp,'prestim',prestim,'poststim',poststim);

  % 3. calculate norm
  planar    = megdata;
  sens      = planar.sensor_info;
  nchan     = length(sens);
  grad1_ind = strmatch('grad1',{sens.typestring});
  grad2_ind = strmatch('grad2',{sens.typestring});
  nloc      = length(grad1_ind);
  planar_time = planar.epochs.time;
  tmpdat    = zeros(nloc,size(planar_time,2), planar.epochs.num_trials);
  for kk = 1:nloc      
    aa           = planar.epochs.data;
    Bx           = aa(grad1_ind(kk),:);
    By           = aa(grad2_ind(kk),:);
    tmpdat(kk,:) = (Bx.^2 + By.^2).^(1/2);
    tmpsens(kk)            =  sens(grad1_ind(kk));
    tmpsens(kk).label      = sprintf('(%s,%s)',sens(grad1_ind(kk)).label,sens(grad2_ind(kk)).label);
    tmpsens(kk).typestring = '2-norm';
    tmpsens(kk).lognum     = min([sens(grad1_ind(kk)).lognum sens(grad2_ind(kk)).lognum]) - 1;
  end

  planar.sensor_info = tmpsens;
  planar.num_sensors = length(tmpsens);
  planar.epochs.data = tmpdat;
  clear tmpdat tmpsens aa Bx By 

  % 4. epoch grad1+grad2 and norm
  megepoch(i)    = ts_epoch_data(megdata,'samples',trigsamp,'prestim',prestim,'poststim',poststim);
  planarepoch(i) = ts_epoch_data(planar ,'samples',trigsamp,'prestim',prestim,'poststim',poststim);
  planarepoch(i) = ts_preproc(planarepoch(i),'baseline_flag',1,'blcwindow',blcwindow);

  megepoch(i).sensor_info = cat(2,megepoch(i).sensor_info,eegepoch(i).sensor_info);
  megepoch(i).num_sensors = length(megepoch(i).sensor_info);
  megepoch(i).epochs.data = cat(1,megepoch(i).epochs.data,eegepoch(i).epochs.data);
end
clear eegepoch

% Combine conditions
%   Condition 1 = MEG (ref: eeg negpeak)
%   Condition 2 = MEG (ref: eeg pospeak)
%   Condition 3 = EEG (negpeak)
%   Condition 4 = EEG (pospeak)

megepoch    = ts_combine_data(megepoch);
planarepoch = ts_combine_data(planarepoch);

% megavg    = ts_trials2avg(megepoch);
% planaravg = ts_trials2avg(planarepoch);
% 
% ts_ezplot(eegepoch);
% ts_ezplot(megavg);
% ts_ezplot(planaravg);

% artifact rejection

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% III. Perform wavelet analysis on SO epochs
clear megdata eegdata planar megavg planaravg

% wavelet analysis
foi = 5:5:100; % Hz, frequencies of interest; ex. alpha band (8-12Hz)
sf  = 7;       % Hz, spectral resolution
tic
tfplanar = ts_freqanalysis_fieldtrip(planarepoch,...
                                          'foi'         ,foi,...
                                          'sf'          ,sf,...
                                          'trials_flag' ,1,...
                                          'save_flag'   ,0,...
                                          'verbose'     ,0);        
toc                                        
tfmeg    = ts_freqanalysis_fieldtrip(megepoch,...
                                          'foi'         ,foi,...
                                          'sf'          ,sf,...
                                          'trials_flag' ,1,...
                                          'save_flag'   ,0,...
                                          'verbose'     ,0);        
toc                                        

badtrials = [69 70];

ts_ezplot(tfplanar,'showlabels','yes','blc',1,'baselinetype','zscore','blcwindow',blcwindow);
ts_ezplot(tfmeg,'layout',params.layout,'showlabels','yes','blc',1,'baselinetype','zscore','blcwindow',blcwindow);

% t = tfplanar.timefreq.time;
% s = tfplanar.sensor_info;
% tmpdat = ts_matrix2epoch(squeeze(mean(tfplanar.timefreq.power,3)),'time',t,'sens',s);
% [data reject_data] = ts_reject(tmpdat,'reject_manual_flag',1);

zlim = [-5 5];

tfplanar = ts_data_selection(tfplanar,'trials',setdiff(1:tfplanar.timefreq.num_trials,badtrials));
ts_ezplot(tfplanar,'showlabels','yes','blc',1,'baselinetype','zscore','blcwindow',blcwindow,'zlim',zlim);

tfmeg = ts_data_selection(tfmeg,'trials',setdiff(1:tfmeg.timefreq.num_trials,badtrials));
ts_ezplot(tfmeg,'showlabels','yes','blc',1,'baselinetype','zscore','blcwindow',blcwindow,'zlim',zlim);


% visualization

% average spectral power over frequency band

% TF-domain artifact rejection

% stats

% plv

% average plv over frequency band

% average plv over roi

% define plv-roi conditions => stats








