function params = SO_params(SubjID)
if nargin == 0
  SubjID = input('Subject ID:');
elseif nargin == 1 && iscell(SubjID)
  SubjID = SubjID{1};
elseif nargin > 1
  error('Parameters can be retrieved for only one subject at a time.');
end

params = []; 

% control flags (note: will always create or load post-ICA data)
params.SO_detections_flag           = 1; % detect SO
params.SO_cluster_detections_flag   = 1; % cluster detections
params.SO_cluster_corrcoef1_flag    = 1; % calc R(delay,distance) per cluster
params.ShowPlots                    = 1;

if (params.SO_detections_flag<params.SO_cluster_detections_flag) || ...
    (params.SO_cluster_detections_flag==0 && any([params.SO_cluster_corrcoef1_flag params.SO_cluster_corrcoef2_flag params.SO_cluster_corrcoef3_flag params.SO_cluster_corrcoef4_flag]))
  error('The first 3 steps must be completed in order & before calculating correlation coefficients.');
end
rootdir         = '/home/halgdev/projects/sleep/MEG/SO';
% rootdir         = '/space/emc2/1/halgdev/projects/sleep/Intracranial_SO';

% Generic parameters
params.chantype = 'eeg';%'grad';%'grad';%'eeg';      % chantype must be 'eeg' or 'grad'
params.peaktype = 'pospeak';%'allpeak';%'pospeak';%'negpeak';  % peaktype to use for aggregate count (pospeak,negpeak,bothpeaks/allpeaks)
params.cluster_count = [];%'eeg_negpeak';%'eeg_pospeak';%'eeg_pospeak';
CountThresh1    = 5;    % 5 for eeg
CountThresh2    = 5;% 5/10   % 10 for grad+allpeak with cluster_count=[]; else 5
params.toilim   = [];
params.cluster_ClusterWindow  = .4; % .4/.2           % cluster window size, sec; .4 eeg, .2 grad
params.layout = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
params.MinChansPerCluster                     = 30;         % 30 M/EEG,20 iEEG

% IMPORTANT!!
% The pre-ICA preprocessing will persist through all stages of the analysis
params.ICA_bpfreq     = [.01 100]; % upper limit should be < (Fs after dsfact / 2)
params.ICA_notchflag  = 0;
params.ICA_blcflag    = 'yes';
params.ICA_blcwindow  = [-inf inf]; % demean the complete time series after concatenation
params.ICA_maxsteps   = 20; % way more than sufficient to ensure convergence in well-conditioned data
  % note that the preprocessing is necessary to make the data well-conditioned
params.ICA_ntrial     = 50; % # sec to display from the total concatenated recording

% SO detection
params.bpfilter          = 1;
params.blc               = 'yes';
params.decimate          = 0;
params.smooth            = 0;
params.hilbertpeaks      = 1;
params.derivpeaks        = 1;
params.onlysinglepeaks   = 0;            % 1-skip cycles w/ multiple peaks; 0-select peak with largest amplitude.
params.zerocross         = 1;
params.monotonic         = 0;
params.gtmedian          = 1;
params.return_zerocross  = 0;
params.bpfreq            = [.01 4];      % Hz
params.blcwindow         = [];           % sec, [] = entire time series; [begin end]
params.decimate_factor   = [];
params.smooth_window     = 0;            % sec, window used for moving average smoothing
params.zero2zero_limits  = [.25 1];      % sec, time bw zero crossings must be within these limits
params.mindist           = 2*min(params.zero2zero_limits);
% sec, (0=skip)(minimum distance bw detections of same polarity)
params.debug             = 0;
params.debug_toilim      = [1210 1230];

% Eliminate isolated detections
params.peakpairs_flag    = 0;
params.peakpairs_tau     = 1;

params.cluster_method             = 'histogram';
params.cluster_thresh             = CountThresh1;%5;%'meanstd3';  % threshold for defining histogram peaks
params.cluster_StepSize           = .01;          % step size in sliding aggregate count
  % 1 / max(Fc) where Fc = cutoffs for the pre-ICA BPF
params.cluster_IntegrationWindow  = .05;         % size of sliding window in aggregate count, sec
  % Rationale: the min time b/w zero-crossings around a detected peak is
  % 250ms (zero2zero_limits) => min period of wave with both peaks detected
  % is 500ms (& fastest detected wave is 2Hz). Since the ClusterWindow is
  % 300ms, detections within 150ms of a peak will be involved in the
  % cluster to which it belongs. That implies there could be a 50ms overlap
  % in clusters defined for the pos and neg peaks of the same wave
  % according to the histogram method of clustering. In order to make sure
  % a trough can be detected between the local histogram peaks
  % (for the pos and neg clusters), the IntegrationWindow should be <<
  % 50ms. Setting it to 50/2 = 25ms will allow for enough overlap between
  % steps for the curve to be smooth and, at the same time, enough distance
  % to prevent smoothing over the closest detected peaks.
params.cluster_MinSeparation      = .05;         % combine peaks closer than this, sec
  % this smooths over ripples near histogram peaks without smoothing over
  % distinct peaks in the slow oscillation. See rationale for the
  % integration window.
  % similar to Murphy, Massimini, Volgushev, and others
params.cluster_peaktype           = params.peaktype;%'both';       % peaktype to use for aggregate count (pospeak,negpeak,bothpeaks/allpeaks)

% Interval selection
% params.MinChansPerCluster                     = 30;         % 20
  % same as Murphy?
params.AASM_EpochLength                       = 30;         % sec
params.IntervalSelection_CountThreshold       = CountThresh2;          % threshold crossing level
  % 1. Stage N3 (AASM scoring): >=20% of a 30sec epoch must contain large
  % amplitude slow waves (.5-2Hz) = (2sec-.5sec).
  % => 1.25sec mean / 20% of 30sec = 6sec => ceil(6/1.25) = 5 waves in 30sec
  % => multiply by 2 since pos & neg peaks could both produce hist peaks
  % => 10 clusters in 30sec
  % NOTE: if zero-to-zero approximates (wave period)/2 then
  % zero2zero_limits = [.25 1] suggests we have selected waves .5-2Hz which
  % complies with the AASM standard.
params.IntervalSelection_CombineIfLessThan    = 45;         % sec
  % combine if 1 epoch drops below thresh between two above
params.IntervalSelection_MinCrossingTime      = 60;         % sec
  % require at least two epochs in a row to exceed thresh
  
% flip matrix
params.flipmatrix_FlipWindow                  = params.cluster_ClusterWindow;%.2;         % sec, make this the same as ClusterWindow
params.flipmatrix_EpochPad                    = .5;         % sec, make this large enough to include the reflocked average peak
params.flipmatrix_MinChanPercentVar           = 'meanstd';  % used to mark bad channels
    % is there a better way to set this param?

% consistency analysis
params.consistency_RiThreshold                = .05; % 'meanstd'; % most inconsistent fraction of channels
params.corrcoef_MinChansPerCluster            = 20;
% params.consistency_RejectAllInconsistent_flag = 0;
% other
params.corrcoef_alpha = .05;  % significance level for the corr coef (dist vs delay) analysis

params.SulciFactor = 3; % scaling factor for converting angular distances to euclidean distances between electrodes

% Subject-specific parameters
switch SubjID
  case {'s1',1}
    % head partially outside dewar; very little frontal coverage.
    params.SubjID     = 's1';
%     params.badlabels  = {'C1','Cz','C6','CPz','CP4'};
%     params.badlabels  = {'MEG 0811','MEG 0813','MEG 0823','MEG 0822','MEG 2132','MEG 2433','MEG 2423','MEG 2422','MEG 2413','MEG 2441','MEG 2443','EEG 033','EEG 034','EEG 037','EEG 045','EEG 047'};
%     params.badlabels  = {'MEG 0811','MEG 0813','MEG 0823','MEG 0822','MEG 2132','MEG 2433','MEG 2423','MEG 2422','MEG 2413','MEG 2441','MEG 2443','EEG 029','EEG 030','EEG 033','EEG 034','EEG 036','EEG 037','EEG 044','EEG 045','EEG 047','EOG 061','EOG 062','EOG 063','EMG 064','EEG 065','EEG 073','EEG 074'};
    params.badlabels  = {'MEG 0123','MEG 0733','MEG 0743','MEG 0813','MEG 0823','MEG 0822','MEG 1022','MEG 1912','MEG 2033','MEG 2132','MEG 2133','MEG 2443','EEG 030','EEG 033','EEG 034','EEG 036','EEG 037','EEG 044','EEG 045','EEG 047','EOG 061','EOG 062','EOG 063','EMG 064','EEG 065','EEG 073','EEG 074'};
    params.eeg_index  = [61 62 1:60];
    params.datafiles  = {...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_1_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_2_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_3_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_4_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_5_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_6_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_DC_s1_7_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_DC_s1_8_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_9_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_10_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_11_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_12_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_13_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles/sleep_s1_14_raw.mat' ...
                        };
    params.matfile_index = 1:14; % 1:12
    params.ICA_dsfact    = 3;
        % EEG: lrate = 1.5e-05, wchange = 54.2 at step 20; removed IC #15 
        params.BrainVolume   = 1450;%1260; % cm^3
  case {'s2',2}
    params.SubjID     = 's2';
    params.badlabels  = {'MEG 2112','EEG 018','EEG 029','EEG 033','EEG 037','EEG 042'};
    params.eeg_index  = [61 62 1:60];
    params.datafiles  = {...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_1_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_2_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_3_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_4_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_5_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_6_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_7_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_8_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_9_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_10_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_11_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_DC_s2_12_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_DC_s2_13_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_14_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_15_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_16_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s2/matfiles/sleep_s2_17_raw.mat' ...
                        };
    params.matfile_index = 1:11;
    params.ICA_dsfact    = 3;
    % ICA was VERY good at removing EKG from grads for s2.
    params.BrainVolume   = 1450;%1260; % cm^3
  case {'s3',3}
    error('This subject was removed from the analysis because head was outside of the dewar.');
    params.BrainVolume   = 1450;%1130; % cm^3
    params.eeg_index     = [61 62 1:60];
  case {'s4',4}
    params.SubjID     = 's4';
    params.badlabels  = {'EEG 009','EEG 033','EEG 045','EEG 055'};
    params.eeg_index  = [61 62 1:60];
    params.datafiles  = {...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_1_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_2_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_3_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_4_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_5_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_6_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_7_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_8_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_9_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_10_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_11_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_12_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_13_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_14_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_15_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_16_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_17_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_18_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_19_raw.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s4/matfiles/sleep_s4_20_raw.mat' ...
                        };
    params.matfile_index = 1:2; % 1:19
    params.ICA_dsfact    = 3;
      % ICA showed several components that may have had EKG but also signal.
    % None were removed. The interval displayed was ok. There was a problem
    % saving the post-ICA data because of it's size. I copied
    % proc_epoch_data_1.mat to proc_epoch_data_ICA.mat b/c they should be
    % the same (since no ICs were removed).
    params.BrainVolume   = 1450;%1260; % cm^3
  case {'s5',5}
    % SL_0 is the initial file which was not named with any # at beginning
    params.SubjID = 's5';
    params.badlabels  = {'MEG 2043','EEG 014','EEG 028','EEG 048','EEG 060'};
    params.eeg_index  = [];
    params.datafiles  = {...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s5/matfiles/SL_0_nd01_060726.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s5/matfiles/SL_1_nd01_060726.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s5/matfiles/SL_2_nd01_060726.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s5/matfiles/SL_3_nd01_060726.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s5/matfiles/SL_4_nd01_060726.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s5/matfiles/SL_5_nd01_060726.mat' ...
                        };
    params.matfile_index = 1:2; % 1:6
    params.ICA_dsfact    = 5;
    % ICA showed several components that may have had EKG but also signal.
    % None were removed. The interval displayed was ok.
    params.BrainVolume   = 1450;%1260; % cm^3
  case {'s6',6}
    % Note: experiment was stopped in the middle of SL_3 until subject was sleepy
    % => exclude SL_3
    params.SubjID = 's6';
    params.badlabels  = {'MEG 2221','MEG 0532','MEG 2043','MEG 2243','EEG 60'};
    params.eeg_index  = [];
    params.datafiles  = {...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_0_ma01_060729.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_1_ma01_060729.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_3_ma01_060729.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_4_ma01_060729.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_5_ma01_060729.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_6_ma01_060729.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s6/matfiles/SL_7_ma01_060729.mat' ...
                        };
    params.matfile_index = 1:2; %[1 2 4 5 6 7];
    params.ICA_dsfact    = 5;
      % ICA was VERY good at removing EKG from grads for s6.
      params.BrainVolume   = 1450;%1260; % cm^3
  case {'s7',7}
    error('This subject was removed from the analysis because never entered SWS');
%     % Note: Subject did not go to deep sleep and last few recordings are awake signals
%     % => exclude 6 & 7
%     params.SubjID = 's7';
%     params.badlabels  = {'MEG 2412','MEG 2413','MEG 2043','EEG 003','EEG 004',...
%       'EEG 014','EEG 016','EEG 028','EEG 030','EEG 048','EEG 060'};
%     params.datafiles  = {...
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_1_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_2_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_3_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_4_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_5_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_6_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_7_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/awake_tb01_060801.fif' ...  
%     '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/emptyroom_tb01_060801.fif'};
%     params.matfile_index = 1:5;
%     params.ICA_dsfact    = 5;
%       % ICA was OK at removing EKG from grads in s7. The IC (#11) had what
%       % looked like noise superimposed on the EKG (hopefully not neural
%       % activity) and the EKG spike amplitudes varied slightly over time.
%       % => failed to find any N3 intervals; signal dominated by STRONG
%       % ~.17Hz (T~6sec) artifact present in all channels.
%       % Maybe: repeat ICA and try to remove the very low-freq artifact.
%       % - I used the rate of EKG peaks in the IC showing it to determine
%       % the what stretch along the x-axis corresponded to ~6sec and found
%       % the IC showing a large wave with that period. Two ICs showed EKG
%       % (15 and 52) and one showed a clear 6sec period (4). I removed all
%       % three.
%       % Visual inspection of the data showed no SWS which is consistent with 
%       % the technician report. The large .17Hz wave still dominated.
%       % CONCLUSION: Exclude this subject from the SO analysis.
    params.BrainVolume   = 1450;%1260; % cm^3
    params.eeg_index     = [];
  case {'s8',8}
    params.SubjID = 's8';
    params.badlabels  = {'MEG 0532','MEG 2043','MEG 2221','MEG 2243','EEG 035','EEG 060'};
    params.eeg_index  = [];
    params.datafiles  = {...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s8/matfiles/SL_1_nb01_060808.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s8/matfiles/SL_2_nb01_060808.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s8/matfiles/SL_3_nb01_060808.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s8/matfiles/SL_4_nb01_060808.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s8/matfiles/SL_5_nb01_060808.mat' ...
                        '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s8/matfiles/SL_6_nb01_060808.mat' ...
                        };
    params.matfile_index = 1:6;
    params.ICA_dsfact    = 5;
        % no ICs were removed because the time interval displayed was not useful. 
        % TODO: Need to rerun with different ntrial (try .005) and w/o displaying edge effects.
    params.BrainVolume   = 1450;%1260; % cm^3
  case {'ch17','CH17',17}
    % intracranial: 2 grids + seeg/strips?
    rootdir       = '/space/emc2/1/halgdev/projects/sleep/Intracranial_SO';
    params.SubjID = 'CH17';
    params.badlabels = [];
    params.datafiles = { ...
            '/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 4010953________ 2005AUG30 101426 CF s.mat' ...
						'/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 4010953________ 2005AUG30 101641 CF s.mat' ...
						'/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 4010953________ 2005AUG30 101922 CF s.mat' ...
						'/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 4010953________ 2005AUG30 101950 CF s.mat' ...
						'/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 4010953________ 2005AUG30 102256 CF s.mat' ...
						'/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 4010953________ 2005AUG30 102559 CF s.mat' ...
            '/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 401 09 53______ 2005AUG29 134122 C1 s.mat' ...
						'/space/emc2/1/halgdev/projects/sleep/Intracranial_SO/CH17/matfiles/Clau 401 09 53______ 2005AUG29 135854 C1 s.mat' ...
            };    
    params.matfile_index = 7:8; % 1:6 (portable), 7:8 (central)
    params.ICA_dsfact    = 1;
    params.BrainVolume   = 1450;%1260; % cm^3
    params.iEEG_flag     = 1;
    % central layouts (2 grids)
    params.layout        = {'/space/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg/CH17/CH17_central_grid1.lay',...
                            '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/ieeg/CH17/CH17_central_grid2.lay'};
  otherwise
    error('Subject ID not recognized.');
end
params.SubjDir  = [rootdir '/' params.SubjID];

%% References & params from lit
% 1. Stage N3 (AASM scoring): >=20% of a 30sec epoch must contain large
% amplitude slow waves (.5-2Hz) = (2sec-.5sec).
% => 1.25sec mean / 20% of 30sec = 6sec => ceil(6/1.25) = 5 waves in 30sec
% => multiply by 2 since pos & neg peaks could both produce hist peaks
% => 10 waves in 30sec

%% Brain volume
% http://hypertextbook.com/facts/2001/ViktoriyaShchupak.shtml

