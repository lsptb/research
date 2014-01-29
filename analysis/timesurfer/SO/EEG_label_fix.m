%% Try to map EEG indices in raw to fixed data by visual inspection
% 06-Oct-2010 by Jason Sherfey

% Test data
% !cd /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/debug/eeg_label_fix
% !cp /space/mdeh1/10/halgdev/projects/jsherfey/sleep/sleep_fixed/s1/sleep_1_2_raw.fif .
% !cp /space/mdeh1/10/halgdev/projects/jsherfey/sleep/sleep_fixed/s1/sleep_s1_2_filt_raw.fif .

cd /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/debug/eeg_label_fix

rawdata = ts_MNE_loadfif('sleep_1_2_raw.fif');
fixdata = ts_MNE_loadfif('sleep_s1_2_filt_raw.fif');

rawdata = ts_data_selection(rawdata,'chantype','eeg');
fixdata = ts_data_selection(fixdata,'chantype','eeg');

rawsens = rawdata.sensor_info;
fixsens = fixdata.sensor_info;

[sel2,raw2fix] = match_str({fixsens.label},{rawsens.label}); 
  % isequal({fixsens.label},{rawsens(raw2fix).label})
  % isequal(raw2fix,[61 62 1:60]')

rawloc = {rawsens.loc};
fixloc = {fixsens.loc};
rawlocfix = rawloc(raw2fix);

rawdata.epochs.time = rawdata.epochs.time - min(rawdata.epochs.time);
visualizer(rawdata);

raw2fix             = [61 62 1:60];
rawdata.sensor_info = rawdata.sensor_info(raw2fix);
rawdata.num_sensors = length(rawdata.sensor_info);
rawdata.epochs.data = rawdata.epochs.data(raw2fix,:);

% /space/mdeh1/10/halgdev/projects/jsherfey/sleep/sleep_all/sleep_subj3_1/sleep_s3_3_raw.fif


%% Fix files that already exist
raw2fix = [61 62 1:60];
file    = 'proc_eeg_epoch_data_ICA.mat';
load(file); % epoch_data
epoch_data.sensor_info = epoch_data.sensor_info(raw2fix);
epoch_data.num_sensors = length(epoch_data.sensor_info);
epoch_data.epochs.data = epoch_data.epochs.data(raw2fix,:);
save(file,'epoch_data','-v7.3');

raw2fix = [61 62 1:60];
file    = 'SO_eeg_detections.mat';
load(file); % events, detections, params
events      = events(raw2fix);
detections  = detections(raw2fix);
save(file,'events','detections','params');

