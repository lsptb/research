%% Purpose: prepare bipolar data from continuous or epoched data

cd /space/mdkm1/2/kmdev/projects/jsherfey/sleep/KC/eeg

%% MG30 Bipolar channel pairs

bipolarlist = { ...
              'FP1','F7'; ...
              'F7','T3'; ...
              'T3','T5'; ...
              'T5','O1'; ...
              'FP2','F8'; ...
              'F8','T4'; ...
              'T4','T6'; ...
              'T6','O2'; ...
              'FP1','F3'; ...
              'F3','C3'; ...
              'C3','P3'; ...
              'P3','O1'; ...
              'FP2','F4'; ...
              'F4','C4'; ...
              'C4','P4'; ...
              'P4','O2'; ...
              'FZ','CZ'; ...
              'CZ','PZ'; ...
              'LPF1','LPF2'; ...
              'LPF2','LPF3'; ...
              'LPF3','LPF4'; ...
              'LPF4','LPF5'; ...
              'LPF5','LPF6'; ...
              'LSF1','LSF2'; ...
              'LSF2','LSF3'; ...
              'LSF3','LSF4'; ...
              'LSF4','LSF5'; ...
              'LSF5','LSF6'; ...
              'LCIN1','LCIN2'; ...
              'LCIN2','LCIN3'; ...
              'LCIN3','LCIN4'; ...
              'LCIN4','LCIN5'; ...
              'LCIN5','LCIN6'; ...
              'LAT1','LAT2'; ...
              'LAT2','LAT3'; ...
              'LAT3','LAT4'; ...
              'LAT4','LAT5'; ...
              'LAT5','LAT6'; ...
              'LAT6','LAT7'; ...
              'LAT7','LAT8'; ...
              'LPT1','LPT2'; ...
              'LPT2','LPT3'; ...
              'LPT3','LPT4'; ...
              'LPT4','LPT5'; ...
              'LPT5','LPT6'; ...
              'LPT6','LPT7'; ...
              'LPT7','LPT8'; ...
              'RPF1','RPF2'; ...
              'RPF2','RPF3'; ...
              'RPF3','RPF4'; ...
              'RPF4','RPF5'; ...
              'RPF5','RPF6'; ...
              'RSF1','RSF2'; ...
              'RSF2','RSF3'; ...
              'RSF3','RSF4'; ...
              'RSF4','RSF5'; ...
              'RSF5','RSF6'; ...
              'RCIN1','RCIN2'; ...
              'RCIN2','RCIN3'; ...
              'RCIN3','RCIN4'; ...
              'RCIN4','RCIN5'; ...
              'RCIN5','RCIN6'; ...
              'RAT1','RAT2'; ...
              'RAT2','RAT3'; ...
              'RAT3','RAT4'; ...
              'RAT4','RAT5'; ...
              'RAT5','RAT6'; ...
              'RAT6','RAT7'; ...
              'RAT7','RAT8'; ...
              'RPT1','RPT2'; ...
              'RPT2','RPT3'; ...
              'RPT3','RPT4'; ...
              'RPT4','RPT5'; ...
              'RPT5','RPT6'; ...
              'RPT6','RPT7'; ...
              'RPT7','RPT8'; ...
              };

%% Continuous data           
% Load continuous data and set up epoch_data structure
file = '/space/md3/5/halgdev/projects/Spindles/Intracranial_Data/MG30/MG30_Sleep1.mat';
load(file,'data'); data.epochs = data.cont; data = rmfield(data,'cont');

% Prepare bipolar data
sens    = {data.sensor_info.label};
elec1   = cellfun(@(x)strmatch(x,sens,'exact'),bipolarlist(:,1));
elec2   = cellfun(@(x)strmatch(x,sens,'exact'),bipolarlist(:,2));
bipolar = data.epochs.data(elec1,:,:) - data.epochs.data(elec2,:,:);
bipolar = ts_matrix2epoch(bipolar,'time',data.epochs.time,'continuous',1);
label   = cellfun(@(x,y)[x '-' y],bipolarlist(:,1),bipolarlist(:,2),'uniformoutput',false);
[bipolar.sensor_info.label] = deal(label{:});

visualizer(data);     % view original continuous data
visualizer(bipolar);  % view bipolar data

%% Epoch data
% load epoch_data
file  = 'MG30_KC_Choosen_from_Scalp_long.eeg';
data  = ts_iEEG_eeg2epoch(file);

% Prepare bipolar data
sens    = {data.sensor_info.label};
elec1   = cellfun(@(x)strmatch(x,sens,'exact'),bipolarlist(:,1));
elec2   = cellfun(@(x)strmatch(x,sens,'exact'),bipolarlist(:,2));
bipolar = data.epochs.data(elec1,:,:) - data.epochs.data(elec2,:,:);
bipolar = ts_matrix2epoch(bipolar,'time',data.epochs.time);
label   = cellfun(@(x,y)[x '-' y],bipolarlist(:,1),bipolarlist(:,2),'uniformoutput',false);
[bipolar.sensor_info.label] = deal(label{:});

% Average data
bp  = ts_data_selection(bipolar,'toilim',[-.5 1.45]);
bp  = ts_preproc(bp,'blc','yes','bpfilter','yes','bpfreq',[.1 30],'bandpass_detrend_flag',0);
avg = ts_trials2avg(bp);
visualizer(avg);  % view average of bipolar epochs

