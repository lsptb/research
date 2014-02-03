% Requires: MNE (fif) or ntools (NSpike), Brian's functions, onestream, MMPS
cd /home/jsherfey/svn/dev/plot_contrasts

% parameters
chanlabel = {'MEG 0113'};%,'MEG 0122','MEG 0132'};
datafile  = '/space/mdkm1/10/kmdev/projects/AM/AM_raw/AM_27j_raw.fif';
prefix    = 'proc_contrasts';
trigchan  = 'STI 014';
foi       = 'MEGlong';
sf        = [];
events    = [1 5 7];
contrasts = {[1 5],[1 7],[5 7]};
blcwindow = [-.25 0];

data   = ts_process_data(datafile,'trigchan',trigchan,'chanlabels',chanlabel,'bpfilter','yes','bpfreq',[.1 55],'blc',1,'blcwindow',blcwindow);
tfavg  = ts_freqanalysis(data,'foi',foi,'sf',sf,'events',events,'trials_flag',1,'save_flag',1,'prefix',prefix,'overwrite',1);
tfwave = ts_freqband_average(tfavg,'blc',1,'baselinetype','zscore','blcwindow',blcwindow);
[tfwave,tfrej] = ts_reject(tfwave,'reject_auto_flag',1,'reject_grad',5,'prescale_grad',1);
tfwave = ts_data_selection(tfwave,'reject_data',tfrej);

for k = 1:length(contrasts)
  erp_stats(k)    = ts_statistics_wrapper(data,'events',contrasts{k});
  tfwave_stats(k) = ts_statistics_wrapper(tfwave,'events',contrasts{k});
end

alldata{1}  = ts_trials2avg(data);
allstats{1} = erp_stats;
alldata{2}  = ts_trials2avg(tfwave);
allstats{2} = tfwave_stats;

ts_plotcontrasts(alldata,allstats,'contrasts',contrasts,'channels',1);

