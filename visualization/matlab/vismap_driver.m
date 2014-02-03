cd /space/mdeh3/9/halgdev/projects/jsherfey/code/visualizer

addpath /space/mdeh3/9/halgdev/projects/jsherfey/sleep/shared
addpath /space/mdeh3/9/halgdev/projects/jsherfey/sleep/meeg/tmp
rootdir = '/space/mdeh3/9/halgdev/projects/jsherfey/sleep';
datapath  = [rootdir '/ieeg/CH17/matfiles/central'];

subjnum = 17;
params  = SO_params(subjnum,[rootdir '/ieeg']);

combolayout = '/space/mdeh3/9/halgdev/projects/jsherfey/sleep/so_paper1/iEEG/CH17/CH17_central_allchans_RefLabel.lay';

toilim = [8500 10000];

% Load iEEG data (continuous epoch_data structure, ntrial=1)
load([datapath '/proc_eeg_epoch_data_ICA.mat'],'epoch_data');
toilim = [8500 10000];

% I = 1000:5000;
% data = epoch_data;
% data.epochs.time = data.epochs.time(I);
% data.epochs.data = data.epochs.data(:,I);
data = ts_data_selection(epoch_data,'toilim',toilim); 
clear epoch_data
for ch = 1:data.num_sensors
  data.epochs.data(ch,:) = ts_freq_filt(data.epochs.data(ch,:)',data.sfreq,[.1 30],[0 0],'bandpass')';
end
dsfact = 4;
data.epochs.time = data.epochs.time(1:dsfact:end);
data.epochs.data = data.epochs.data(:,1:dsfact:end);
data.sfreq = data.sfreq / dsfact;
vismap(data,combolayout);

%% reorder channels for visualizer
chans = ...
  [81 89 97 105 113 121 65 73 ...
   82 90 98 106 114 122 66 74 ...
   83 91 99 107 115 123 67 75 ...
   84 92 100 108 116 124 68 76 ...
   85 93 101 109 117 125 69 77 ...
   86 94 102 110 118 126 70 78 ...
   87 95 103 111 119 127 71 79 ...
   88 96 104 112 120 128 72 80 ...
   64 56 48 40 32 24 16 8 ...
   63 55 47 39 31 23 15 7 ...
   62 54 46 38 30 22 14 6 ...
   61 53 45 37 29 21 13 5 ...
   60 52 44 36 28 20 12 4 ...
   59 51 43 35 27 19 11 3 ...
   58 50 42 34 26 18 10 2 ...
   57 49 41 33 25 17 9 1];

data2 = data;
data2.sensor_info = data2.sensor_info(chans);
data2.epochs.data = data2.epochs.data(chans,:);
visualizer(data2);





