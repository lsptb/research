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
params.SO_detections_flag          = 1;
params.SO_cluster_detections_flag  = 1;
params.SO_cluster_corrcoef1_flag   = 1; % no flip
params.SO_cluster_corrcoef2_flag   = 1; % flip matrix
params.SO_cluster_corrcoef3_flag   = 1; % hist window
params.SO_cluster_corrcoef4_flag   = 1; % consistency

if (params.SO_detections_flag<params.SO_cluster_detections_flag) || ...
    (params.SO_cluster_detections_flag==0 && any([params.SO_cluster_corrcoef1_flag params.SO_cluster_corrcoef2_flag params.SO_cluster_corrcoef3_flag params.SO_cluster_corrcoef4_flag]))
  error('The first 3 steps must be completed in order & before calculating correlation coefficients.');
end

% Subject-specific parameters
switch SubjID
  case {'s1',1}
    params.SubjID     = 's1';
    params.badlabels  = {'C1','Cz','C6','CPz','CP4'};
    params.datafiles  = {...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_2_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_3_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_4_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_fixed/s1/sleep_1_5_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_6_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_DC_s1_7_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_DC_s1_8_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_9_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_1/sleep_s1_10_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_11_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_12_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_13_raw.fif' ...
        '/space/emc2/2/halgdev/data/MEG_MGH/sleep_subj1_2/sleep_s1_14_raw.fif'};
    params.matfile_index = 1:5;
    params.ICA_dsfact    = 3;
  case {'s2',2}
    params.SubjID     = 's2';
    params.badlabels  = {'MEG 2112','EEG 018','EEG 029','EEG 033','EEG 037','EEG 042'};
    params.datafiles  = {...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_1_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_2_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_3_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_4_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_5_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_6_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_7_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_8_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_1/sleep_s2_9_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_10_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_11_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_DC_s2_12_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_DC_s2_13_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_14_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_15_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_16_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj2_2/sleep_s2_17_raw.fif'};      
    params.matfile_index = 1:9;
    params.ICA_dsfact    = 3;
  case {'s3',3}
    error('This subject was removed from the analysis because head was outside of the dewar.');
  case {'s4',4}
    params.SubjID     = 's4';
    params.badlabels  = {'EEG 009','EEG 033','EEG 045','EEG 055'};
    params.datafiles  = {...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_1_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_2_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_3_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_4_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_5_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_6_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_7_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_8_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_9_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_10_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_11_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_12_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_13_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_14_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_15_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_16_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_17_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_18_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_19_raw.fif' ...
        '/home/halgdev/data/MEG_MGH/sleep_subj4_2/sleep_s4_20_raw.fif'};
%         '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_sleep_noise.fif' ...
%         '/home/halgdev/data/MEG_MGH/sleep_subj4_1/sleep_s4_sleep_noise_filt.fif' ...
    params.matfile_index = 1:12; 
      % WHETHER THESE ARE THE BEST FILES SHOULD BE VERIFIED!
    params.ICA_dsfact    = 3;
  case {'s5',5}
    params.SubjID = 's5';
    params.badlabels  = {'MEG 2043','EEG 014','EEG 028','EEG 048','EEG 060'};
    params.datafiles  = {...
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_0_nd01_060726.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_1_nd01_060726.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_2_nd01_060726.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_3_nd01_060726.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_4_nd01_060726.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/SL_5_nd01_060726.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/wake_nd01_060726.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_nd01_060726/emptyroom_nd01_060726.fif'};      
    params.matfile_index = 1:6;
    params.ICA_dsfact    = 5;
  case {'s6',6}
    params.SubjID = 's6';
    params.badlabels  = {'MEG 2221','MEG 0532','MEG 2043','MEG 2243','EEG 60'};
    params.datafiles  = {...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_0_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_1_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_3_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_4_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_5_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_6_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/SL_7_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/awake_ma01_060729.fif' ...
    '/home/halgdev/data/MEG_UCSD/SL_ma01_060729/emptyroom_ma01_060729.fif'};
    params.matfile_index = 1:7;
    params.ICA_dsfact    = 5;
  case {'s7',7}
    params.SubjID = 's7';
    params.badlabels  = {'MEG 2412','MEG 2413','MEG 2043','EEG 003','EEG 004',...
      'EEG 014','EEG 016','EEG 028','EEG 030','EEG 048','EEG 060'};
    params.datafiles  = {...
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_1_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_2_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_3_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_4_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_5_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_6_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/SL_7_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/awake_tb01_060801.fif' ...  
    '/home/halgdev/data/MEG_UCSD/SL_tb01_060801/emptyroom_tb01_060801.fif'};
    params.matfile_index = 1:7;
    params.ICA_dsfact    = 5;
  case {'s8',8}
    params.SubjID = 's8';
    params.badlabels  = {'MEG 0532','MEG 2043','MEG 2221','MEG 2243','EEG 035','EEG 060'};
    params.datafiles  = {...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_1_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_2_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_3_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_4_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_5_nb01_060808.fif' ...
        '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/SL_6_nb01_060808.fif'};
    %     '/home/halgdev/data/MEG_UCSD/SL_nb01_060808/emptyroom_nb01_060808.fif' ...
    params.matfile_index = 1:6;
    params.ICA_dsfact    = 5;
  otherwise
    error('Subject ID not recognized.');
end

% Generic parameters
rootdir         = '/space/emc2/1/halgdev/projects/sleep/MEG/SO';
params.SubjDir  = [rootdir '/' params.SubjID];
params.chantype = {'grad1','grad2'};      % chantype must be a cell array of strings
params.toilim   = [];%[1200 1500];

% IMPORTANT!!
% The pre-ICA preprocessing will persist through all stages of the analysis
% params.ICA_dsfact     = 3;
params.ICA_bpfreq     = [.1 100]; % upper limit should be < (Fs after dsfact / 2)
params.ICA_notchflag  = 0;
params.ICA_blcflag    = 'yes';
params.ICA_blcwindow  = [-inf inf]; % demean the complete time series
params.ICA_maxsteps   = 20; % way more than sufficient to ensure convergence in well-conditioned data
  % note that the preprocessing is necessary to make the data well-conditioned
params.ICA_ntrial     = .01; % this should be a small fraction of the total concatenated recording

% SO detection
params.bpfilter          = 1;
params.blc               = 'yes';
params.decimate          = 0;
params.smooth            = 0;
params.hilbertpeaks      = 1;
params.derivpeaks        = 1;
params.onlysinglepeaks   = 0;            % 1-skip cycles w/ multiple peaks; 0-select peak with largest amplitude.
params.zerocross         = 1;
params.monotonic         = 1;
params.gtmedian          = 1;
params.return_zerocross  = 0;
params.bpfreq            = [.1 4];      % Hz
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
params.cluster_thresh             = 'meanstd';    % threshold for defining histogram peaks
params.cluster_StepSize           = .01;          % step size in sliding aggregate count
  % 1 / max(Fc) where Fc = cutoffs for the pre-ICA BPF
params.cluster_IntegrationWindow  = .025;         % size of sliding window in aggregate count, sec
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
params.cluster_MinSeparation      = .025;         % combine peaks closer than this, sec
  % this smooths over ripples near histogram peaks without smoothing over
  % distinct peaks in the slow oscillation. See rationale for the
  % integration window.
params.cluster_ClusterWindow      = .3;           % cluster window size, sec
  % similar to Murphy, Massimini, Volgushev, and others
params.cluster_peaktype           = 'both';       % peaktype to use for aggregate count

% Interval selection
params.MinChansPerCluster                     = .25;        % fraction of nchan
  % same as Murphy?
params.AASM_EpochLength                       = 30;         % sec
params.IntervalSelection_CountThreshold       = 10;         % threshold crossing level
  % 1. Stage N3 (AASM scoring): >=20% of a 30sec epoch must contain large
  % amplitude slow waves (.5-2Hz) = (2sec-.5sec).
  % => 1.25sec mean / 20% of 30sec = 6sec => ceil(6/1.25) = 5 waves in 30sec
  % => multiply by 2 since pos & neg peaks could both produce hist peaks
  % => 10 waves in 30sec
params.IntervalSelection_CombineIfLessThan    = 45;         % sec
  % combine if 1 epoch drops below thresh between two above
params.IntervalSelection_MinCrossingTime      = 60;         % sec
  % require at least two epochs in a row to exceed thresh
  
% flip matrix
params.flipmatrix_FlipWindow                  = .3;         % sec, make this the same as ClusterWindow
params.flipmatrix_EpochPad                    = .5;         % sec, make this large enough to include the reflocked average peak
params.flipmatrix_MinChanPercentVar           = 'meanstd';  % used to mark bad channels
    % is there a better way to set this param?
  
% histogram method
% 

% consistency analysis
params.consistency_RiThreshold                = 'meanstd';
    % is there a better way to set this param?
params.consistency_RejectAllInconsistent_flag = 1;

% other
params.corrcoef_alpha = .05;  % significance level for the corr coef (dist vs delay) analysis
params.layout         = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';

%% References & params from lit
% 1. Stage N3 (AASM scoring): >=20% of a 30sec epoch must contain large
% amplitude slow waves (.5-2Hz) = (2sec-.5sec).
% => 1.25sec mean / 20% of 30sec = 6sec => ceil(6/1.25) = 5 waves in 30sec
% => multiply by 2 since pos & neg peaks could both produce hist peaks
% => 10 waves in 30sec
% 2. 


