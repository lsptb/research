parms = [];
parms.SensorEventFile         = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8/s8_SO_init_peaks_filt0.01-4Hz_toi600-1350_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat';
parms.FlipFile                = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8/s8_SO_flip_matrix_filt0.01-4Hz_toi600-1350_grad1grad2_refall_smooth0.1sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_nodoublepks_13-Jun-2010.mat';
parms.CorrSensorEventFile     = '';
parms.ClusterSensorEventFile  = '';%fakefile';
 
parms.rootoutdir        = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg';
parms.writelog          = 0;
parms.overwrite         = 0;

% Data
usedatainmemory         = 0;
parms.clearall          = ~usedatainmemory;
parms.loadflag          = ~usedatainmemory;
parms.toiflag           = ~usedatainmemory;
parms.preprocflag       = ~usedatainmemory;
parms.SubjectID         = 's8';
parms.matfile_index     = 1:3;
parms.toilim            = [600 1350];
parms.RSSflag           = 0;

% SO detection
parms.detectionflag     = 1;
parms.bpfilter          = 1;
parms.blc               = 'yes';
parms.decimate          = 0;
parms.smooth            = 1;
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
parms.smooth_window     = .05;            % sec, window used for moving average smoothing
parms.zero2zero_limits  = [.25 1];        % sec, time bw zero crossings must be within these limits
parms.mindist           = 2*min(parms.zero2zero_limits);%.2;             % sec, (0=skip)(minimum distance bw detections of same polarity)

parms.debug             = 0;
parms.debug_toilim      = [1210 1230];

parms.peakpairs_flag    = 1; % only process pos/neg peak pairs
parms.peakpairs_tau     = 1; % seconds, find pairs with |tk1-tk2|<tau

% Epoching, averaging, & flipping
parms.calc_flip_matrix  = 1; % this takes a long time if above FlipFile=''
parms.refavgflag        = 1;
% parms.noflipflag        = 1;
% parms.flipflag          = 1;
% parms.chantype          = {'grad1'};%{'grad1','grad2'};      % which chans to epoch/avg
parms.peaktype          = {'negpeak'};%{'negpeak','pospeak'};  % which ref peaks to use
parms.Ref               = {'MEG0143','MEG0213','MEG1323','MEG1423','MEG0723','MEG2122'};%,'MEG0723'}{'MEG0213','MEG0212','MEG0143','MEG0142'};  % index, cellarr of labels, 'all','max' (max # of detections)
parms.detrend           = 0;                      % whether to detrend before filtering before epoching
% parms.blc               = 'no';                   % whether to blc after epoching (using all times)
parms.EpochPad          = 1;                      % sec, epoch padding
parms.FlipWindow        = .2;                     % sec, look at peaks in tk+/-FlipWindow when determining peak for comparing polarities
parms.CoRefAvgWindow    = .2;                     % sec, chan k must have detection within tj+/-tau for the j-th ref detection to be used
% note: filtering or not (and fc) is inherited from SO detection params
parms.plotflag          = 1;
parms.saveflag          = 0;
parms.closeflag         = 0;
parms.autoscale         = 0;
parms.allaxes           = [];

% Processing flipped data
parms.corrdetectionflag = 1;
parms.clusterflag       = 1;
parms.ClusterWindow     = .2;                     % sec, cluster detections within window tk+/-ClusterWindow
parms.PSD_foilim        = [20 60];
parms.PSD_foi           = 10:100;          % if non-empty, PSD is calculated at foi
parms.PSD_offset        = 'epochpad';     % 'epochpad',array [ms], or []
parms.PSD_fraction      = 150/2000;             % fraction of the epoch over which to calc FFT
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

