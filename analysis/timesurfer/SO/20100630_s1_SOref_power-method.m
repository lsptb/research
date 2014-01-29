% MEG 0133, 0143, 0142, 0213
% t = 2100-2200 sec
tic
toilim = [2100 2200];
labels = {data.sensor_info.label}; %{'MEG 0133','MEG 0143','MEG 0142','MEG 0213'};
seldat = ts_data_selection(data,'chanlabel',labels,'toilim',toilim);
nn     = length(labels);

solim  = [.1 1.5];
foi    = .1:.1:3;
sf     = 3;

seldat = ts_preproc(seldat,'blc','yes','bpfilter','yes','bpfreq',[.01 4]);
tfdata = ts_freqanalysis_fieldtrip(seldat,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0);
% visualizer(seldat,tfdata)
zdata  = ts_zscore(tfdata,'baselinetype','zscore','blcwindow','all');
% visualizer(seldat,zdata)

tfmean = tfdata; tfmean.sensor_info(nn).label='mean power';
tfmean.timefreq.power(nn,:,:) = mean(tfdata.timefreq.power,1); % repeat with sum
% visualizer(seldat,tfmean)

fix    = find(foi>=solim(1)&foi<=solim(2));
tmpdat = seldat; tmpdat.sensor_info(nn).label='mean power';
tmpdat.epochs.data(nn,:) = mean(squeeze(tfmean.timefreq.power(nn,fix,:)),2);
% visualizer(tmpdat)

visualizer(tmpdat,tfmean)
toc



% dat = tmpdat.epochs.data(nn,:);
% StdDevWindow = min(50,floor(length(dat)/10));
% StdDevWindow = StdDevWindow + ~mod(StdDevWindow,2);
% npoints      = length(dat) - StdDevWindow + 1;
% offset       = floor(StdDevWindow/2);
% tmp          = zeros(size(dat));
% for l = 1:npoints
%   tmp(l+offset) = std(dat(l:l+StdDevWindow-1));
% end
% dat = tmp; clear tmp
% dat = smooth(dat,50,'lowess');
% tmpdat.epochs.data(nn,:) = dat;
