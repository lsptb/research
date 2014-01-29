parms = [];
parms.SensorEventFile         = '';
parms.SensorEventPhaseFile    = '';
parms.FlipFile                = '';
parms.CorrSensorEventFile     = '';
parms.ClusterSensorEventFile  = '';
 
parms.rootoutdir        = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg';
parms.writelog          = 1;
parms.overwrite         = 0;

% Data
usedatainmemory         = 0;
parms.clearall          = ~usedatainmemory;
parms.loadflag          = ~usedatainmemory;
parms.toiflag           = ~usedatainmemory;
parms.ICAflag           = ~usedatainmemory;
parms.preprocflag       = 0;%~usedatainmemory;
parms.SubjectID         = 's2';
parms.matfile_index     = 1:9;%17;
parms.toilim            = 'all';
parms.RSSflag           = 0;

% SO detection
parms.detectionflag     = 1;
parms.bpfilter          = 1;
parms.blc               = 'yes';
parms.decimate          = 0;
parms.smooth            = 0;
parms.hilbertpeaks      = 1;
parms.derivpeaks        = 1;
parms.onlysinglepeaks   = 0;              % 1-skip cycles w/ multiple peaks; 0-select peak with largest amplitude.
parms.zerocross         = 1;
parms.monotonic         = 1;
parms.gtmedian          = 1;
parms.return_zerocross  = 0;

parms.bpfreq            = [.01 4];%       % Hz
parms.blcwindow         = [];             % sec, [] = entire time series; [begin end]
parms.decimate_factor   = [];
parms.smooth_window     = 0;            % sec, window used for moving average smoothing
parms.zero2zero_limits  = [.25 1];        % sec, time bw zero crossings must be within these limits
parms.mindist           = 2*min(parms.zero2zero_limits);%.2;             % sec, (0=skip)(minimum distance bw detections of same polarity)

parms.debug             = 0;
parms.debug_toilim      = [1210 1230];

parms.peakpairs_flag    = 0; % only process pos/neg peak pairs
parms.peakpairs_tau     = 1; % seconds, find pairs with |tk1-tk2|<tau
parms.post_toilim       = [800 2400];%[12 91] * 30; % 30sec hypnogram epochs 12-91

% Epoching, averaging, & flipping
parms.calc_flip_matrix  = 0; % this takes a long time if above FlipFile=''
parms.refavgflag        = 1;
parms.peaktype          = {'negpeak','pospeak'};%{};%{'negpeak'};%{'negpeak','pospeak'};  % which ref peaks to use
parms.Ref               = [];%{};%{'MEG0143','MEG0213','MEG1323','MEG1423','MEG0723','MEG2122'};%,'MEG0723'}{'MEG0213','MEG0212','MEG0143','MEG0142'};  % index, cellarr of labels, 'all','max' (max # of detections)
parms.detrend           = 0;                      % whether to detrend before filtering before epoching
% parms.blc               = 'no';                   % whether to blc after epoching (using all times)
parms.EpochPad          = 1;                      % sec, epoch padding
parms.FlipWindow        = .2;                     % sec, look at peaks in tk+/-FlipWindow when determining peak for comparing polarities
parms.CoRefAvgWindow    = .2;                     % sec, chan k must have detection within tj+/-tau for the j-th ref detection to be used
% note: filtering or not (and fc) is inherited from SO detection params
parms.plotflag          = 0;
parms.saveflag          = 1;
parms.closeflag         = 0;
parms.autoscale         = 0;
parms.allaxes           = [];

% Processing flipped data
parms.corrdetectionflag = 0;
parms.clusterflag       = 0;
parms.ClusterWindow     = .2;                     % sec, cluster detections within window tk+/-ClusterWindow
parms.PSD_foilim        = [20 60];
parms.PSD_foi           = 10:100;         % if non-empty, PSD is calculated at foi
parms.PSD_offset        = .075;              % sec, center of psd window wrt refpeak
parms.PSD_pad           = .075;           % sec
  % PSD_window = refpeak + PSD_offset + [-PSD_pad:PSD_pad]
parms.Coh_flag          = 0;
parms.Coh_foilim        = [20 60];
parms.Coh_NFFTfactor    = 1/20; % NFFT = round(Fs*Coh_NFFTfactor)
parms.Coh_Noverlap      = [];
parms.Coh_time          = 'psd';          % 'epoch','cluster','psd'

% TF analysis
parms.timefreqflag      = 0;
parms.dsfact            = 4;
parms.tftoilim          = [1080 1250]; % [1050 1300];
parms.foi               = [2:4:54 66:4:100];
parms.sf                = 6;
parms.fband             = [10 100];  % TF wave

