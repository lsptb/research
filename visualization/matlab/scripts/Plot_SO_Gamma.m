clear all

SubjID    = 's1';
FileIndex = 2;

gamma_freq = [20 80];
blcwindow  = [-.5 .5]; % sec, mu/std window for z-score calculation
% [-.5 .5], [-.8 -.5]

reject_eeg    = 350;
reject_grad   = 3000;
reject_gamma  = 2.5; % z-score wrt (all time & trials) per (chan,freq)

TF_zscore_clim      = [-5 5];
TF_zscore_blcwindow = [-.5 .5];

addpath(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/%s/gamma/scripts',SubjID));
addpath(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/%s/gamma/functions',SubjID));
cd(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/%s/gamma/matfiles',SubjID));

run(sprintf('%s_SO_params',SubjID));
[fpath,fname] = fileparts(params.sws_info(FileIndex).datafile);
prefix = strrep(fname,'_raw','');

%% Artifact rejection
load([prefix '_sws_SO_avg_data.mat'],'avg_data');
load([prefix '_sws_SO_epoch_data.mat'],'epoch_data','params');
load([prefix '_sws_SO_tf_epoch_data_005-100Hz.mat'],'tf_epoch_data');

t     = tf_epoch_data.timefreq(1).time;
f     = tf_epoch_data.timefreq(1).frequencies;
nchan = tf_epoch_data.num_sensors;
ncond = length(tf_epoch_data.timefreq);
ntime = length(t);
nfreq = length(f);
Tix   = find(t>=blcwindow(1)&t<=blcwindow(2));  % for plotting, not rejection
Fix   = find(f>=gamma_freq(1)&f<=gamma_freq(2));

% Gamma z-score calculation for rejection
%   - use one mean & stdev for all times and trials at a given freq & chan
gamma_epoch_data  = rmsubfield(epoch_data,'epochs.data');
gamma_avg_data    = rmsubfield(avg_data,'averages.data');
[gamma_epoch_data.epochs.time] = deal(t);
[gamma_avg_data.averages.time] = deal(t);
for c = 1:ncond
  fprintf('cond %g of %g\n',c,ncond);
  max_ntrial = tf_epoch_data.timefreq(c).num_trials;
  ntrials    = tf_epoch_data.timefreq(c).num_trials_per_chan;
  gamma_mat  = nan(nchan,ntime,max_ntrial,'single');
  gamma_avg  = zeros(nchan,ntime,'single');
  for ch = 1:nchan
    fprintf('chan %g of %g\n',ch,nchan);
    trl = 1:ntrials(ch);
    dat = double(tf_epoch_data.timefreq(c).power(ch,:,Fix,trl));
    clear mu sd
    % calc one mean & stdev per freq for this chan
    for freq = 1:length(Fix)
      tmp          = dat(1,:,freq,:);
      mu(1,1,freq) = mean(tmp(:));
      sd(1,1,freq) = std(tmp(:));
    end
    mu  = repmat(mu,[1 ntime 1 ntrials(ch)]);
    sd  = repmat(sd,[1 ntime 1 ntrials(ch)]);    
    Z   = (dat - mu) ./ sd;
    Z   = mean(Z,3);
    Z   = permute(Z,[1 2 4 3]);
    gamma_mat(ch,:,trl) = single(Z);
    gamma_avg(ch,:)     = single(mean(Z,3));
    clear tmp mu sd Z dat
  end
  gamma_epoch_data.epochs(c).data = gamma_mat;
  gamma_avg_data.averages(c).data = gamma_avg;
  clear gamma_mat gamma_avg
end

% add fake trial_info
if ~isfield(epoch_data.epochs,'trial_info')
  for c = 1:ncond
    tmp = ts_matrix2data(rand(2,10000),'continuous',1); 
    tmp = ts_epoch_data(tmp,'samples',500+[1:epoch_data.epochs(c).num_trials]);
    trl = tmp.epochs.trial_info;
    epoch_data.epochs(c).trial_info = trl;
    gamma_epoch_data.epochs(c).trial_info = trl;
  end
end

% perform rejection on each channel separately
reject_data = []; goodtrials = [];
for c = 1:ncond
  fprintf('cond %g of %g\n',c,ncond);
  ntrials = tf_epoch_data.timefreq(c).num_trials_per_chan;
  for ch  = 1:nchan
    fprintf('chan %g of %g\n',ch,nchan);
    tmp   = ts_data_selection(epoch_data,'conditions',c,'channels',ch,'trials',1:ntrials(ch),'verbose',0);
    [jnk,rej1] = ts_reject(tmp,'prescale_eeg',10^6,'reject_eeg',reject_eeg,'prescale_grad',10^13,'reject_grad',reject_grad,'reject_auto_flag',1,'verbose',0);
    tmp   = ts_data_selection(gamma_epoch_data,'conditions',c,'channels',ch,'trials',1:ntrials(ch),'verbose',0);
    [jnk,rej2] = ts_reject(tmp,'prescale_eeg',1,'reject_eeg',reject_gamma,'prescale_grad',1,'reject_grad',reject_gamma,'reject_auto_flag',1,'verbose',0);
    reject_data(ch).badtrials{c} = unique([rej1.badtrials{1} rej2.badtrials{1}]);
    goodtrials{c,ch} = setdiff(1:ntrials(ch),reject_data(ch).badtrials{c});
  end
end

% reject_data.badchans    = [];
% reject_data.badtrials   = badtrials;
% reject_data.event_code  = [1 2];

% averaging without artifacts (epoch_data, tf_epoch_data, gamma_epoch_data)
load([prefix '_sws_SO_tf_avg_data_005-100Hz.mat'],'tf_avg_data');
avg_data_rej        = rmsubfield(avg_data,'averages.data');
gamma_avg_data_rej  = rmsubfield(avg_data,'averages.data');
tf_avg_data_rej     = rmsubfield(tf_avg_data,'timefreq.power');
[gamma_avg_data_rej.averages.time] = deal(gamma_epoch_data.epochs(1).time);
for c = 1:ncond
  fprintf('cond %g of %g\n',c,ncond);
  ntrials                             = tf_epoch_data.timefreq(c).num_trials_per_chan;
  avg_data_rej.averages(c).data       = zeros(nchan,numel(avg_data_rej.averages(c).time),'single');
  gamma_avg_data_rej.averages(c).data = zeros(nchan,ntime,'single');
  tf_avg_data_rej.timefreq(c).power   = zeros(nchan,ntime,nfreq,'single');
  for ch  = 1:nchan
    fprintf('chan %g of %g\n',ch,nchan);
    avg_data_rej.averages(c).data(ch,:)       = mean(epoch_data.epochs(c).data(ch,:,goodtrials{c,ch}),3);
    gamma_avg_data_rej.averages(c).data(ch,:) = mean(gamma_epoch_data.epochs(c).data(ch,:,goodtrials{c,ch}),3);
    tf_avg_data_rej.timefreq(c).power(ch,:,:) = mean(tf_epoch_data.timefreq(c).power(ch,:,:,goodtrials{c,ch}),4);    
  end
end

% TIME DOMAIN  
title_string = sprintf('%s (%g-%gsec)',prefix,params.sws_info(FileIndex).toilim);

cond_labels = {'negEEG','posEEG'};
ts_ezplot(avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',title_string,'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',title_string,'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',title_string,'autoscale',0,'chantype','grad2','layout',meg_layout);

ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad2','layout',meg_layout);

ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad2','layout',meg_layout);

% TIME-FREQUENCY DOMAIN
zlim = TF_zscore_clim;
blim = TF_zscore_blcwindow;
cond_labels = {'negEEG','posEEG'};
ts_ezplot(tf_avg_data_rej,'events',[1 2],'title',title_string,'showlabels','yes','autoscale',0,'chantype','eeg','blc','yes','baselinetype','zscore','blcwindow',blim,'zlim',zlim,'layout',eeg_layout,'cond_labels',cond_labels);
ts_ezplot(tf_avg_data_rej,'events',[1 2],'title',title_string,'showlabels','yes','autoscale',0,'chantype','grad1','blc','yes','baselinetype','zscore','blcwindow',blim,'zlim',zlim,'layout',meg_layout,'cond_labels',cond_labels);
ts_ezplot(tf_avg_data_rej,'events',[1 2],'title',title_string,'showlabels','yes','autoscale',0,'chantype','grad2','blc','yes','baselinetype','zscore','blcwindow',blim,'zlim',zlim,'layout',meg_layout,'cond_labels',cond_labels);

% tmpdata = avg_data_rej;
% tmpdata.averages([3 4]) = gamma_avg_data_rej.averages;
% tmpdata.averages(3).event_code = 3;
% tmpdata.averages(4).event_code = 4;
% for k=1:4,tmpdata.averages(k).data = tmpdata.averages(k).data / max(abs(tmpdata.averages(k).data(:))); end
% ts_ezplot(tmpdata,'showlabels','yes','cond_labels',{'negEEG (ERP)','posEEG (ERP)','negEEG (gamma)','posEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','eeg','layout',eeg_layout);
% 
% tmpdata = avg_data_rej;
% tmpdata.averages([3 4]) = gamma_avg_data_rej.averages;
% tmpdata.averages(3).event_code = 3;
% tmpdata.averages(4).event_code = 4;
% for k=1:4
%   for ch = 1:tmpdata.num_sensors
%     tmpdata.averages(k).data(ch,:) = tmpdata.averages(k).data(ch,:) / max(abs(tmpdata.averages(k).data(ch,:))); 
%   end
% end
% ts_ezplot(tmpdata,'showlabels','yes','cond_labels',{'negEEG (ERP)','posEEG (ERP)','negEEG (gamma)','posEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','eeg','layout',eeg_layout);

tmpdata = avg_data_rej;
tmpdata.averages = [tmpdata.averages gamma_avg_data_rej.averages];
tmpdata.averages(3).event_code = 3;
tmpdata.averages(4).event_code = 4;
for k=1:4
  for ch = 1:tmpdata.num_sensors
    tmpdata.averages(k).data(ch,:) = tmpdata.averages(k).data(ch,:) / max(abs(tmpdata.averages(k).data(ch,:))); 
  end
end
tmpdata = ts_preproc(tmpdata,'blc','yes');
ts_ezplot(tmpdata,'showlabels','yes','events',[1 2 3],'linewidth',2,'graphcolor','kbr','cond_labels',{'negEEG (ERP)','posEEG (ERP)','negEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','eeg','layout',eeg_layout);
ts_ezplot(tmpdata,'showlabels','yes','events',[1 2 4],'linewidth',2,'graphcolor','kbr','cond_labels',{'negEEG (ERP)','posEEG (ERP)','posEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','eeg','layout',eeg_layout);

ts_ezplot(tmpdata,'showlabels','yes','events',[1 2 3],'linewidth',1.5,'graphcolor','kbr','cond_labels',{'negEEG (ERP)','posEEG (ERP)','negEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','grad1','layout',meg_layout);
ts_ezplot(tmpdata,'showlabels','yes','events',[1 2 4],'linewidth',1.5,'graphcolor','kbr','cond_labels',{'negEEG (ERP)','posEEG (ERP)','posEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','grad1','layout',meg_layout);

ts_ezplot(tmpdata,'showlabels','yes','events',[1 2 3],'linewidth',1.5,'graphcolor','kbr','cond_labels',{'negEEG (ERP)','posEEG (ERP)','negEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','grad2','layout',meg_layout);
ts_ezplot(tmpdata,'showlabels','yes','events',[1 2 4],'linewidth',1.5,'graphcolor','kbr','cond_labels',{'negEEG (ERP)','posEEG (ERP)','posEEG (gamma)'},'title',[title_string ' (ERP & gamma)'],'autoscale',1,'chantype','grad2','layout',meg_layout);



%% gamma calculation for plotting

gamma_epoch_data  = rmsubfield(epoch_data,'epochs.data');
gamma_avg_data    = rmsubfield(avg_data,'averages.data');
[gamma_epoch_data.epochs.time] = deal(t);
[gamma_avg_data.averages.time] = deal(t);
for c = 1:ncond
  fprintf('cond %g of %g\n',c,ncond);
  max_ntrial = tf_epoch_data.timefreq(c).num_trials;
  ntrials    = tf_epoch_data.timefreq(c).num_trials_per_chan;
  gamma_mat  = nan(nchan,ntime,max_ntrial,'single');
  gamma_avg  = zeros(nchan,ntime,'single');
  for ch = 1:nchan
    trl = 1:ntrials(ch);
    tmp = double(tf_epoch_data.timefreq(c).power(ch,:,Fix,trl));
    mu  = mean(tmp(1,Tix,:,:),2);
    sd  = std (tmp(1,Tix,:,:),0,2);
    Z   = (tmp - repmat(mu,[1 ntime 1 1])) ./ repmat(sd,[1 ntime 1 1]);
    Z   = mean(Z,3);
    Z   = permute(Z,[1 2 4 3]);
    gamma_mat(ch,:,trl) = single(Z);
    gamma_avg(ch,:)     = single(mean(Z,3));
    clear tmp mu sd Z
  end
  gamma_epoch_data.epochs(c).data = gamma_mat;
  gamma_avg_data.averages(c).data = gamma_avg;
  clear gamma_mat gamma_avg
end

% averaging without artifacts (epoch_data, tf_epoch_data, gamma_epoch_data)
load([prefix '_sws_SO_tf_avg_data_005-100Hz.mat'],'tf_avg_data');
gamma_avg_data_rej                    = rmsubfield(avg_data,'averages.data');
[gamma_avg_data_rej.averages.time]    = deal(gamma_epoch_data.epochs(1).time);
for c = 1:ncond
  ntrials                             = tf_epoch_data.timefreq(c).num_trials_per_chan;
  gamma_avg_data_rej.averages(c).data = zeros(nchan,ntime,'single');
  for ch  = 1:nchan
    gamma_avg_data_rej.averages(c).data(ch,:) = mean(gamma_epoch_data.epochs(c).data(ch,:,goodtrials{c,ch}),3);
  end
end

ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad2','layout',meg_layout);

ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad2','layout',meg_layout);


%% yet another gamma calculation for plotting

tf_avg_data_rej_zscore = ts_zscore(ts_data_selection(tf_avg_data_rej,'events',[1 2]),'baselinetype','zscore','blcwindow',blim);
c=1; gamma_avg_data_rej.averages(c).data = mean(tf_avg_data_rej_zscore.timefreq(c).power,3);
c=2; gamma_avg_data_rej.averages(c).data = mean(tf_avg_data_rej_zscore.timefreq(c).power,3);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(gamma_avg_data_rej,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad2','layout',meg_layout);

tf_avg_data_zscore = ts_zscore(ts_data_selection(tf_avg_data,'events',[1 2]),'baselinetype','zscore','blcwindow',blim);
c=1; gamma_avg_data.averages(c).data = mean(tf_avg_data_zscore.timefreq(c).power,3);
c=2; gamma_avg_data.averages(c).data = mean(tf_avg_data_zscore.timefreq(c).power,3);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','eeg','layout',eeg_layout);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad1','layout',meg_layout);
ts_ezplot(gamma_avg_data,'showlabels','yes','cond_labels',cond_labels,'title',[title_string ' (gamma)'],'autoscale',0,'chantype','grad2','layout',meg_layout);

