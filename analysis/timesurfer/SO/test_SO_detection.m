fid = 1;
load /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8/matfiles/SL_2_nb01_060808_grad.mat % data
dat  = ts_data_selection(data,'chantype','grad1');

parms = [];
parms.bpfilter          = 1;
parms.decimate          = 1;
parms.smooth            = 1;
parms.hilbertpeaks      = 1;
parms.derivpeaks        = 0;
parms.zerocross         = 1;
parms.monotonic         = 1;
parms.gtmedian          = 1;
parms.return_zerocross  = 0;

parms.bpfreq            = [.3 5];   % Hz
parms.decimate_factor   = 4;
parms.smooth_window     = .05;      % sec
parms.zero2zero_limits  = [.25 1];  % sec

parms.toilim = [600 1350];
tic
args  = mmil_parms2args(parms);
peaks = SO_detection(dat,args{:});
toc

tmpdat   = ts_data_selection(dat,'toilim',parms.toilim);
t        = tmpdat.epochs.time;
 events  = [];
 for k   = 1:length(peaks)
   events(k).label = dat.sensor_info(k).label;
   events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
   events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
 end
 

peaks_dec  = peaks;
events_dec = events;

parms.decimate = 0;
tic
args  = mmil_parms2args(parms);
peaks = SO_detection(dat,args{:});
toc

tmpdat   = ts_data_selection(dat,'toilim',parms.toilim);
t        = tmpdat.epochs.time;
 events  = [];
 for k   = 1:length(peaks)
   events(k).label = dat.sensor_info(k).label;
   events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
   events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
 end

 