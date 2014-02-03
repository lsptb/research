lay = {'/space/mdkm1/2/kmdev/projects/jsherfey/sleep/KC/CH17_central_grid1.lay',...
       '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/KC/CH17_central_grid2.lay'};

% Load data
file  = 'CH17_evoked_KC_portable_telemetry_BS55-65IIR48dB.eeg';
trls  = [2 3 5 6 8:10 12:14 25:34 36:39 41];

file  = 'CH17_evoked_KC_central_telemetry_BS55-65IIR48dB.eeg';
trls  = [44 46 48:51 53:57 59:65];

data  = ts_iEEG_eeg2epoch(file);
labels={data.sensor_info.label};
dat   = ts_data_selection(data,'trials',trls);
T     = dat.epochs.time;
catdat = dat;
catdat.epochs.data = reshape(catdat.epochs.data,[dat.num_sensors length(T)*dat.epochs.num_trials]);
catdat.epochs.time = (1:size(catdat.epochs.data,2))/data.sfreq;
dat   = ts_preproc(dat,'bpfilter','yes','bpfreq',[.1 4],'bandpass_detrend_flag',0,'blc','yes');
chans = 65:128; lbl = labels(chans);
% dat  = ts_data_selection(dat,'channels',1:64);
dat  = ts_data_selection(dat,'channels',chans);
figure
a = nearest(T,.2);
b = nearest(T,1);
clear ind
for k = 1:dat.epochs.num_trials
  mat = dat.epochs.data(:,:,k);
  for j = 1:size(mat,1)
    ind(j) = find(mat(j,a:b)==min(mat(j,a:b)),1);
  end
  [ind,I] = sort(ind);
  mat     = mat(I,:);
  subplot(2,1,1),imagesc(T,1:dat.num_sensors,mat); 
  title(sprintf('trial %g',trls(k)));
  set(gca,'xlim',[-2 5]); 
  ylim = get(gca,'ylim');
  ytick = linspace(ylim(1),ylim(2),dat.num_sensors);
  set(gca,'YTick',ytick);
  set(gca,'YTickLabel',lbl(I));  
  subplot(2,1,2),plot(T,mat(1:40,:));axis tight; set(gca,'xlim',[-.5 3]); 
  tmp = catdat;
  tmp.epochs.data = tmp.epochs.data(I,:);
  tmp.sensor_info = tmp.sensor_info(I);
%   visualizer(tmp);
  pause
%   close(findobj('tag','plots'));
end

dat   = ts_data_selection(data,'trials',trls);

% select one trial and grid
gridnum = 1;
xlim    = [.3 .8]; nr =5; nc = 5; bpfreq=[.01 4];
trl     = 1;
% plotlim = [-1 3];
f1 = figure('color','w','position',[0 150 850 800]);
% f2 = figure('color','w','position',[750 150 950 800]);
f3 = figure('color','w');
for trl = 1:dat.epochs.num_trials
  tmp = ts_data_selection(dat,'toilim',xlim,'trials',trl,'channels',64*(gridnum-1) + 1:64);
  t = tmp.epochs.time;
  % filter
  tmp = ts_preproc(tmp,'bpfilter','yes','bpfreq',bpfreq,'bandpass_detrend_flag',0,'blc','yes');
  mat = tmp.epochs.data;
  % keep channels with abs(values) > mean(abs all chans) + 2*std(abs all chans)
  vec = mat(:);
  vec = vec(vec<0);
  mu  = mean(vec);
  sd  = std(vec);
%   % ts_ezplot(tmp,'showlabels','yes');
%   badchans = find(any(mat<mu-sd,2)==0);
  mat(badchans,:) = 0;
%   % tmp.epochs.data = mat;
%   % ts_ezplot(tmp,'showlabels','yes');
  % find peaks
%   tk   = zeros(tmp.num_sensors,1);
%   for k = 1:tmp.num_sensors
%     if ismember(k,badchans), continue; end
%     level    = mu - sd;
%     [ind,t0] = crossing(mat(k,:),t,level);
% %   %   plot(t,mat(k,:),'-',t(ind),mat(k,ind),'*'); hline(level,'k'); hline(-level,'k');
% %   %   title(num2str(k))
% %   %   pause
%     if length(ind) < 2
%       tk(k) = find(min(mat(k,:))==mat(k,:),1);
%     else
%       tk(k) = find(min(mat(k,ind(1):ind(2)))==mat(k,:),1);
%     end
%   end
%   del = ts_data_selection(tmp,'badchans',badchans);
%   del.epochs.time = 0;
%   del.epochs.data = t(tk(tk~=0))';
%   tmp.epochs.data(badchans,:) = 0;
  figure(f1); clf
  ts_ezplot(tmp,'showlabels','yes','title',sprintf('trial %g',trls(trl)),'newfig',0);
%   tmp = ts_data_selection(tmp,'badchans',badchans);
%   figure(f2); clf
%   topo(del,'layout',lay{1},'iEEG_flag',1,'showlabels','yes','fig',f2); 
%   set(gca,'clim',xlim); colorbar
%   title(sprintf('delay map (trial %g)',trls(trl)));
%   figure(f3); clf
%   chans = 1:del.num_sensors;
%   pts = cellfun(@(x,y)mat(x,y),num2cell(chans),num2cell(tk(tk~=0)'));
%   plot(t,mat,'-',t(tk(tk~=0)),pts,'*'); axis tight; title('KC overlay');
  topo(tmp,'xlim',xlim,'nrows',nr,'ncols',nc,'layout',lay{gridnum},'iEEG_flag',1,'showlabels','no','fig',f3); colorbar
  set(gcf,'name',sprintf('amplitude (trial %g)',trls(trl)));
  pause
end


xlim = [0 1]; nr = 5; nc = 5; 
thistrl = 1;
ieegdat = ts_data_selection(dat,'channels',1:64,'toilim',xlim,'trials',thistrl);
topo(ieegdat,'xlim',xlim,'nrows',nr,'ncols',nc,'layout',lay{1},'iEEG_flag',1) % ,'showlabels','yes'
% ts_ezplot(ieegdat,'showlabels','yes');
ieegdat = ts_data_selection(dat,'channels',65:128,'toilim',xlim);
figure; topo(ieegdat,'xlim',xlim,'nrows',nr,'ncols',nc,'layout',lay{2},'iEEG_flag',1,'showlabels','yes')
% ts_ezplot(ieegdat,'showlabels','yes');

% Filter before averaging
proc  = ts_preproc(dat,'bpfilter','yes','bpfreq',[.1 30],'bandpass_detrend_flag',0,'blc','yes');
% Average over multiple K-Complexes
avg   = ts_trials2avg(proc);
visualizer(avg);


% Display all K-Complexes per channel
out = ts_visual_reject(dat,'method','channel');
% Display all channels per K-Complex
out = ts_visual_reject(ts_data_selection(dat,'channels',[1:124]),'method','trial');

% Time-frequency analysis of one K-Complex
KCnumber = 1;   % the K-Complex epoch number
tfr = ts_freqanalysis_fieldtrip(ts_data_selection(dat,'trials',KCnumber),'foi',[.5 1:30],'sf',1,'trials_flag',1,'save_flag',0);
% Display spectral power z-score with the single K-Complex waveform
visualizer(ts_data_selection(dat,'trials',KCnumber),ts_zscore(tfr));
