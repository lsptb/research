function plot_ECoG_Events(data,varargin)
% xlim:     interval for sorting based on peaks
% tlim:     interval to display
% bpfreq:   filter corner freqs
% nshow:    number of channels to overlay
% toprows:  number of topoplot rows
% topcols:  number of topoplot columns
% avgovertime: whether to average over time or plot single time point (yes or no)
% layout:   grid layout
%           note: the first layout should correspond to the clusters structure
params = mmil_args2parms( varargin, ...
                         { 'xlim',[],[],...
                           'tlim',[],[],...
                           'prestim',1,[],...
                           'poststim',2,[],...
                           'bpfreq',[],[],...
                           'nshow',30,[],...
                           'toprows',4,[],...
                           'topcols',4,[],...
                           'avgovertime','yes',[],...
                           'layout',[],[],...
                           'markers',[],[],...
                           'clusters',[],[],...
                           'trials',[],[],...
                           'zlim',[],[],...
                           },false);
                           
if isempty(params.layout) || (isempty(params.markers) && isempty(params.clusters))
  fprintf('Both a grid layout and times for epoching must be specified.\n');
  return;
end
if ~iscell(params.layout),  params.layout = {params.layout}; end
if isstruct(params.clusters)
  cst = params.clusters;
  if isempty(params.trials), params.trials = 1:length(cst.epochs); end
  params.markers = [cst.epochs(params.trials).RefTime];
%   params.markers = [cst.epochs(params.trials).HistTime];
else
  cst = [];
end
if isempty(params.trials), params.trials = 1:length(params.markers); end

for l = 1:length(params.layout)
  [chNum,X,Y,Width,Height,Lbl,Rem] = textread(params.layout{l},'%f %f %f %f %f %q %q');
  for i=1:length(Lbl)
    if ~isempty(Rem{i})
      % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
      Lbl{i} = [Lbl{i} ' ' Rem{i}];
    end
  end
  [sel1,sel2] = match_str({data.sensor_info.label},Lbl);
  g{l} = sel1;
end

tlim = params.tlim; if isempty(tlim), tlim = [-params.prestim params.poststim]; end
xlim = params.xlim; if isempty(xlim), xlim = [-params.prestim params.poststim]; end
nr = params.toprows; nc = params.topcols; lay = params.layout;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels  = {data.sensor_info.label};
T       = data.epochs.time;
Fs      = data.sfreq;
prepad  = round(params.prestim*Fs);
postpad = round(params.poststim*Fs);
ntrial  = length(params.markers);

% data   = ts_iEEG_eeg2epoch(file);
% dat    = ts_data_selection(data,'trials',trls);
% dat    = ts_preproc(dat,'bpfilter','yes','bpfreq',bpfreq,'bandpass_detrend_flag',0,'blc','yes');
% labels ={dat.sensor_info.label};
% T      = dat.epochs.time;
% if ~isempty(params.bpfreq)
%   data   = ts_preproc(data,'bpfilter','yes','bpfreq',params.bpfreq,'bandpass_detrend_flag',0,'bandpass_baseline_flag',1);
% end

sz = get(0,'ScreenSize'); w = sz(3); h = sz(4); 
L  = min(2*w/5,h/2);
if 1
  f1 = figure(1); if isempty(get(f1,'CurrentAxes')), set(gcf,'color','w','position',[0 h/5 w/5 4*h/5]); end
  f6 = figure(6); if isempty(get(f6,'CurrentAxes')), set(gcf,'color','w','position',[w/5 L L L]); end
  if length(params.layout)>1, f3 = figure(3); if isempty(get(f3,'CurrentAxes')), set(gcf,'color','w','position',[w/5+L,w/5+L L L]); end; end
  f4 = figure(4); if isempty(get(f4,'CurrentAxes')), set(gcf,'color','w','position',[w/5 0 L L]); end
  if length(params.layout)>1, f5 = figure(5); if isempty(get(f5,'CurrentAxes')), set(gcf,'color','w','position',[w/5+L 0 L L]); end; end
  if ~isempty(cst), f2 = figure(2); if isempty(get(f2,'CurrentAxes')), set(gcf,'color','w','position',[0 0 w/5 h/5]); end; end
end

for k = 1:ntrial
  trl = params.trials(k);
  tk  = params.markers(k);
  tmp = ts_data_selection(data,'toilim',tk + [-params.prestim params.poststim]);
  tmp.epochs.time = tmp.epochs.time - tk;
  t   = tmp.epochs.time;
  a   = nearest(t,xlim(1));
  b   = nearest(t,xlim(2));
  tmp.epochs.data = blc(tmp.epochs.data);
  if ~isempty(params.bpfreq)
    tmp.epochs.data = ts_freq_filt(tmp.epochs.data',Fs,params.bpfreq,[0 0],'bandpass')';
  else
    tmp.epochs.data = ts_freq_filt(tmp.epochs.data',Fs,60,5,'notch')';
  end
  
  % sort channels based on neg peak
  mat = tmp.epochs.data;
%   clear ind
%   for j = 1:size(mat,1)
%     ind(j) = find(mat(j,:)==min(mat(j,a:b)),1);
%   end
%   delays  = t(ind);
%   delays  = delays - min(delays);
%   [ind,I] = sort(ind);
%   mat     = mat(I,:);
  
  % plot sorted data
  figure(f1); clf
  subplot(2,1,1),imagesc(t,1:tmp.num_sensors,mat); vline(0,'k');
  title(sprintf('trial %g (t = %5.5g sec)',trl,tk));
  set(gca,'xlim',tlim);
  ylim = get(gca,'ylim');
%   ytick = linspace(ylim(1),ylim(2),tmp.num_sensors);
%   set(gca,'YTick',ytick);
%   set(gca,'YTickLabel',labels(I));  
  subplot(2,1,2),plot(t,mat(1:min(params.nshow,size(mat,1)),:));axis tight; 
  set(gca,'xlim',tlim); vline(0,'k');
  
  %  delay maps
  if ~isempty(cst)
    del = ts_data_selection(tmp,'chanlabels',cst.epochs(trl).InvolvedChans);
    [sel1,sel2] = match_str({del.sensor_info.label},cst.epochs(trl).InvolvedChans);
    delays      = cst.epochs(trl).Delays(sel2);
    del.epochs.data = delays';
    del.epochs.time = 0;
    zlim = [min(delays(:)) 1.5*max(delays(:))];
    
    figure(f2); clf
    subplot(2,1,1);
    topo(del,'zlim',zlim,'layout',lay{1},'iEEG_flag',1,'showlabels','yes','fig',0); 
    colorbar; axis square; title('Delay map, involved grid 1 only')
    subplot(2,1,2);
    plot(cst.epochs(trl).Delays,cst.epochs(trl).Distance*cst.deg2meters,'.');
    axis([0 .3 0 max(cellfun(@max,{cst.epochs.Distance}))*cst.deg2meters]); lsline;
    title(sprintf('R = %3.3g (p = %3.3g, N = %g)',cst.epochs(trl).R,cst.epochs(trl).p,cst.epochs(trl).N));
    ylabel('distance'); xlabel('delay (sec)');
  end
  
  % multiplot
  figure(f6); clf
  zlim = [min(mat(:)) max(mat(:))];
  l    = max(abs(zlim));
  zlim = [-l l];
  for layoutnum = 1:length(g)
    if layoutnum>1,figure(f3); clf; end
    ts_ezplot(ts_data_selection(tmp,'channels',g{layoutnum}),'showlabels','yes','title',sprintf('grid 1 (trial %g), REF = %s',trl,params.markers(trl)),'newfig',0,'layout',lay{layoutnum});%,'vline',0); % cst.epochs(trl).RefChan)
  end
%   ts_ezplot(ts_data_selection(tmp,'channels',g{1}),'showlabels','yes','title',sprintf('grid 1 (trial %g), REF = %s',trl,params.markers(trl),'newfig',0,'layout',lay{1});%,'vline',0); % cst.epochs(trl).RefChan)
%   figure(f3); clf
%   ts_ezplot(ts_data_selection(tmp,'channels',g{2}),'showlabels','yes','title',sprintf('grid 2 (trial %g)',trl),'newfig',0,'layout',lay{2});%,'vline',0);
  
  % topoplot
  figure(f4); clf
  if strcmp(params.avgovertime,'no')
    zlim = 1.2*max(abs(tmp.epochs.data(:)));
    zlim = [-zlim zlim];
  else
    zlim = 'absmax';
  end
  if ~isempty(params.zlim), zlim = params.zlim; end
  for layoutnum = 1:length(g)
    if layoutnum>1,figure(f5); clf; end
    tmptmp=ts_data_selection(tmp,'channels',g{layoutnum});
    topo(tmptmp,'avgovertime',params.avgovertime,'zlim',zlim,'nrows',nr,'ncols',nc,'layout',lay{layoutnum},'iEEG_flag',1,'showlabels','no','fig',0); 
    set(gcf,'name',sprintf('amplitude (trial %g)',trl)); colorbar
  end  
%   tmptmp=ts_data_selection(tmp,'channels',g{1}); topo(tmptmp,'avgovertime',params.avgovertime,'zlim',zlim,'nrows',nr,'ncols',nc,'layout',lay{1},'iEEG_flag',1,'showlabels','no','fig',0); 
% %   tmptmp=ts_data_selection(tmp,'channels',g{1}); topo(tmptmp,'highlight',strmatch(cst.epochs(trl).RefChan,{tmptmp.sensor_info.label},'exact'),'avgovertime',params.avgovertime,'zlim',zlim,'nrows',nr,'ncols',nc,'layout',lay{1},'iEEG_flag',1,'showlabels','no','fig',0); 
% 
%   set(gcf,'name',sprintf('grid 1, amplitude (trial %g)',trl)); colorbar
%   figure(f5); clf
%   tmptmp=ts_data_selection(tmp,'channels',g{2}); topo(tmptmp,'avgovertime',params.avgovertime,'zlim',zlim,'nrows',nr,'ncols',nc,'layout',lay{2},'iEEG_flag',1,'showlabels','no','fig',0); 
%   set(gcf,'name',sprintf('grid 2, amplitude (trial %g)',trl)); colorbar
  if ~isempty(cst), figure(f2); end
  pause
end
