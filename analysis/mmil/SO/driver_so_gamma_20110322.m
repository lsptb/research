% 1. TFR continuous EEG and norm(grad1,grad2)
%   foi = .2:.2:3; % 15 freqs
%   foi = 2:2:100; % 50 freqs
%   Note: loop over raw FIF files => load then TF analysis for EEG then norm
%         => one timefreq_data per FIF file per chantype (EEG, norm)
%         (run on cluster)
% 2. Calc mean TFR after epoching wrt refchan and individual channels
% 3. Calc SPLV from continuous TFR b/w all channel pairs in EEG and norm
%   => save one SPLV file for each TF file
%   (repeat for broadband TFR and low-freq TFR)
% 4. SCA on continuous SPLV (low and broad freqs)
%   => extract sync clusters from continuous SCA measures
%   (chan_strength > threshold around SO in EEG => sum over cluster ~ f(t))
%   Note: consider both cluster size and expansion/contraction rates & 
%         potential relationships with detection cluster size and distribution
% 
% Relate mean TFR from (2) and metric from (4) to SO propagation, distribution, 
% and gamma power.

% Memory management:
% - load all files for detecting SO
% - save one events file per raw file
% - process files separately with their separate events files


%% Part I: slow oscillation detection & clustering
% (run on ip101, ip84, ip36, ip56?)

% load EEG

% EEG detections
%   => clear EEG
%   => save EEG detections

% cluster EEG with EEG detections
%   => clear EEG detections
%   => save EEG clusters

% load MEG (grad)

% MEG detections
%   => clear MEG
%   => save MEG detections

% cluster MEG detections with EEG clusters
%   => clear MEG detections
%   => save MEG clusters (and clear)

%% Part II: epoching and splitting data into manageable chunks
% Note: save all epoch_data to single precision w/ Fs > Nyquist

% 1. MEG
% load MEG
% load MEG clusters

% epoch & save cond1: grad1 (100 trials/file)
% epoch & save cond1: grad2 (100 trials/file)

% calculate norm(grad1,grad2) for each 100 trials
%   => save norm epochs (100 trials/file)
%   => clear MEG, MEG clusters, epochs: grad1, grad2, norm

% repeat for cond2

% 2. EEG
% load EEG
% load EEG clusters

% epoch & save cond1 (100 trials/file)
% epoch & save cond2 (100 trials/file)
%   => clear EEG, EEG clusters, epochs: eeg 

%% Part III: time-frequency analysis
% Note save all timefreq_data to single precision w/ Fs = 100Hz
% (run on mmilcluster)

% Note: at this point we have several epoch_data files containing 100
% trials each for two conditions and both MEG (grad1, grad2, norm) and EEG.
% Next we want to perform a wavelet analysis on each one and save out
% corresponding timefreq_data structures. The timefreq_data structures will
% be used for the gamma analysis.
% ---------------------------------------------------------------------
% << one cluster job for each condition, typestring, trial grouping >>
% ---------------------------------------------------------------------

% Load a single epoch_data structure (100 trials, 1 typestring)

% Remove bad channels

% Perform wavelet analysis  (may need to loop over trials)
%   => clear epoch_data

% Downsample timefreq_data and convert to single precision
%   => save timefreq_data (complex)
%   => clear timefreq_data

% ---------------------------------------------------------------
% Each category is divided into several files with ntrials <= 100
% note: use ts_combine_data to concatenate them later
% EEG pospeak: EEG
% EEG pospeak: grad1
% EEG pospeak: grad2
% EEG pospeak: norm
% EEG negpeak: EEG
% EEG negpeak: grad1
% EEG negpeak: grad2
% EEG negpeak: norm


%% Part III: gamma analysis
% ---------------------------------------------------------------
% HOW TO SHOW GAMMA MODULATION?????
% ---------------------------------------------------------------

% rotate and stack spectral power => visualizer
% (what reference? i.e., what is t=0?)

% z-score for each trial

% freqcorr for each trial

% z-score for average spectral power
% (what reference? i.e., what is t=0?)

% freqcorr then average spectral power
% (what reference? i.e., what is t=0?)

% average z-score or freqcorr over frequency band
% (what reference? i.e., what is t=0?)



%% Part IV: ...









