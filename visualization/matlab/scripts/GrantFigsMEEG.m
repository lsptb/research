% Fig X1AB - topoplot detection density (EEG, MEG)
cd /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles
% load SO_clusters_grad_allpeak_noflip
load SO_clusters_grad_allpeak_noflip_ref-eeg_negpeak.mat
% load lastrun_version1/SO_clusters_noflip.mat
megcst = clusters_Near_all(1);
% % load lastrun/SO_clustered_detections.mat
% % megdet = detections;
% % Ndet   = cellfun(@length,{megdet.cluster_number});
% % megbad = find(Ndet==0);
% megcst = clusters_Near_pks(1);
% megcst = clusters_Near_pks(2);
load SO_clusters_eeg_pospeak_noflip
eegcstpos = clusters_Near_pks(1);
% % load SO_clustered_eeg_pospeak_detections.mat
% % eegdetpos = detections;
load SO_clusters_eeg_negpeak_noflip
eegcstneg = clusters_Near_pks(2);
% % load SO_clustered_eeg_negpeak_detections.mat
% % eegdetneg = detections;
load proc_eeg_epoch_data_ICA
eegdata = epoch_data; clear epoch_data
load proc_grad_epoch_data_ICA
megdata = epoch_data; clear epoch_data

%% Detection topoplot
% clustered detections
% figure
clusters=megcst;
% for j = 2:3
%   if j==1,     clusters = megcst;    chantype = 'meg';
%   elseif j==2, clusters = eegcstpos; chantype = 'eegpos';
%   elseif j==3, clusters = eegcstneg; chantype = 'eegneg';
%   end
  sens      = clusters.sensor_info;
  nchan     = length(clusters.sensor_info);
  labels    = {clusters.sensor_info.label};
  origins   = arrayfun(@(x)x.RefChan,clusters.epochs,'uniformoutput',false);
  RefTime   = arrayfun(@(x)x.RefTime,clusters.epochs);
  ntrial    = length(origins);
  FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clusters.epochs);
  MaxDelay  = arrayfun(@(x)max(x.Delays),clusters.epochs);
  AllInvChn   = {clusters.epochs.InvolvedChans};
  AllInvChn   = [AllInvChn{:}];
  OriginCount = zeros(nchan,1);
  InvolvCount = zeros(nchan,1);
  tmp     = {clusters.epochs.InvolvedChans}; 
  invchan = cellfun(@(x)[x{:}],tmp,'uniformoutput',false);
  for k   = 1:nchan
    label = labels{k};
    OriginCount(k) = numel(find(ismember(origins,label)));
    InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
    ix = find(~cellfun(@isempty,regexp(invchan,label)));
  end
  badchans    = find(InvolvCount==0);
    %   if strcmp(chantype,'meg')
    g1 = strmatch('grad1',{sens.typestring}); [badg1,jnk] = match_str({sens(g1).label},{sens(badchans).label}); % badg1=[];
    g2 = strmatch('grad2',{sens.typestring}); [badg2,jnk] = match_str({sens(g2).label},{sens(badchans).label}); % badg2=[];
    Z           = InvolvCount ./ ntrial;  
    Z(isnan(Z)) = 0;
    bin = Z>0;
    avg = [Z(g1) + Z(g2)] ./ (bin(g1) + bin(g2));
    avg(isnan(avg)) = 0;
    tmp = ts_matrix2avg([avg avg],'sens',clusters.sensor_info(g1));
    tmp.averages.data = avg;
    tmp.averages.time = 0;
    badchan = find(avg==0)
%%     figure; ts_ezplot(tmp,'badchans',badchan,'highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
%     title('grad detections'); set(gca,'clim',[0 .5]); axis square
%%  % all detections
  load SO_grad_detections
  Z = (cellfun(@length,{detections.pospeak}) + cellfun(@length,{detections.negpeak})) / 2;
  Z(isnan(Z)) = 0;
  bin = Z>0;
  avg = [Z(g1) + Z(g2)] ./ (bin(g1) + bin(g2));
  avg = avg';
    avg(isnan(avg)) = 0;
    tmp = ts_matrix2avg([avg avg],'sens',clusters.sensor_info(g1));
    tmp.averages.data = avg;
    tmp.averages.time = 0;
    params = SO_params(1);
    badchan = ismember({detections.label},params.badlabels);
    badchan = find(badchan(g1) + badchan(g2) > 0);
%     badchan = find(avg==0)
    zlim = [];
%%    figure; ts_ezplot(tmp,'badchans',badchan,'highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
%%    
  %% Delay maps
 %   PeakTime = 462; % 1833.9;%
  PeakTime = 946.4;%1917.5;%2600.7;%2167.6; % %
  load SO_grad_detections
  ix = nearest([megcst.epochs.HistTime],PeakTime);
  t0 = megcst.epochs(ix).HistTime; PeakTime = t0;
  tmp = ts_data_selection(megdata,'toilim',t0 + [-1 1],'chanlabel',megcst.epochs(ix).InvolvedChans);
  tmp = ts_preproc(tmp,'bpfilter','yes','bpfreq',[.1 4],'bandpass_baseline_flag',1,'bandpass_detrend_flag',0,'verbose',0);
  tmp.epochs.data = megcst.epochs(ix).Delays';
  tmp.epochs.time = 0;
    
    
    [sel1,sel2] = match_str({megdata.sensor_info.label},megcst.epochs(ix).InvolvedChans);
    Z = zeros(megdata.num_sensors,1);
    Z(sel1) = megcst.epochs(ix).Delays(sel2);
    Z(ismember({megdata.sensor_info.label},params.badlabels)) = 0;
    bin = Z>0;
    del = [Z(g1) + Z(g2)] ./ (bin(g1) + bin(g2));
    del(isnan(del)) = 0;
    badchan = find(del==0 | isnan(del));
    tmp = ts_matrix2avg([del del],'sens',megdata.sensor_info(strmatch('grad1',{megdata.sensor_info.typestring})));
    tmp.averages.data = del;
    tmp.averages.time = 0;
    lay='/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
    tmp = ts_data_selection(tmp,'badchans',find([tmp.averages.data]==0));
%   figure; topo(tmp,'chantype','meg','layout',lay,'highlight',[],'showlabels',0,'xlim',[],'nrows',1,'ncols',1,'electrodes','on','style','straight')
    badlabels = [];%{'MEG 1442','MEG 1423','MEG 1533','MEG 0143'};
    figure; topo(tmp,'badlabels',badlabels,'showlabels',0,'electrodes','on','chantype','meg','layout',lay,'highlight',[],'xlim',[],'nrows',1,'ncols',1,'style','straight')
    axis square; colorbar; caxis([0 .25]); title(sprintf('Grad delay map at %gsec',PeakTime));

    
%     [delsort,ind]=sort(tmp2.averages.data);
%     {tmp2.sensor_info(ind).label}'
  
%     typ = {megdata.sensor_info(sel1).typestring};
%     typ = typ(sel2);
%     g1 = strmatch('grad1',{megdata.sensor_info.label}; [badg1,jnk] = match_str(megcst.epochs(ix).InvolvedChans,params.badlabels); % badg1=[];
%     g2 = strmatch('grad2',{megdata.sensor_info.label}; %[badg2,jnk] = match_str({sens(g2).label},{sens(badchans).label}); % badg2=[];
%     Z           = InvolvCount ./ ntrial;  
%     Z(isnan(Z)) = 0;
%     bin = Z>0;
%     avg = [Z(g1) + Z(g2)] ./ (bin(g1) + bin(g2));
%     avg(isnan(avg)) = 0;
%     tmp = ts_matrix2avg([avg avg],'sens',clusters.sensor_info(g1));
  
  
  
%   [badchan,sel] = match_str(megcst.epochs(ix).InvolvedChans,{params.badlabels{:} 'MEG 1542'})
% %  figure; ts_ezplot(tmp,'badchans',badchan,'highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); 
%   lay=[];%'/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
%   figure; topo(tmp,'badlabels',{params.badlabels{:} 'MEG 1542'},'chantype','meg','layout',lay,'highlight',[],'showlabels',1,'xlim',[],'nrows',1,'ncols',1,'electrodes','on','style','straight')
%   figure; topo(tmp,'badlabels',{params.badlabels{:}},'chantype','meg','layout',lay,'highlight',[],'showlabels',1,'xlim',[],'nrows',1,'ncols',1,'electrodes','on','style','straight')
%   zlim = [];
%   badchan = [];
%   dat = ts_data_selection(megdata,'chantype','grad1');
%   dat.epochs.time = 0;
%   [sel1,sel2] = match_str({megdata.sensor_info.label},megcst.epochs(ix).InvolvedChans)
%   Z = zeros(megdata.num_sensors,1);
%   Z(sel1) = megcst.epochs(ix).Delays(sel2);
%   bin = Z>0;
%   del = [Z(g1) + Z(g2)] ./ (bin(g1) + bin(g2));
%   del(isnan(del)) = 0;
%   badchan = find(del==0 | isnan(del));
%   tmp = ts_matrix2avg([del del],'sens',clusters.sensor_info(g1));
%   tmp.averages.data = del;
%   tmp.averages.time = 0;
% %   badchan = [];
%   figure; ts_ezplot(tmp,'badchans',badchan,'highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
%   axis square; colorbar; caxis([0 .15]); title(sprintf('Grad delay map at %gsec',PeakTime));

  figure; plot(megcst.epochs(ix).Delays,megcst.epochs(ix).Distance,'.'); 
  title(sprintf('grad at %gs (R=%g,p=%g)',megcst.epochs(ix).HistTime,megcst.epochs(ix).R,megcst.epochs(ix).p));
  axis([0 .4 0 180]); xlabel('delay (s)'); ylabel('angular distance (deg)');
  
  load SO_eeg_detections
  badchan = find(ismember({detections.label},params.badlabels));
  ix = nearest([eegcstneg.epochs.HistTime],PeakTime);
  t0 = eegcstneg.epochs(ix).HistTime;
  tmp = ts_data_selection(eegdata,'toilim',t0 + [-1 1],'chanlabel',eegcstneg.epochs(ix).InvolvedChans);
  tmp = ts_preproc(tmp,'bpfilter','yes','bpfreq',[.1 4],'bandpass_baseline_flag',1,'bandpass_detrend_flag',0,'verbose',0);
  tmp.epochs.data = eegcstneg.epochs(ix).Delays';
  tmp.epochs.time = 0;

  figure; topo(tmp,'badlabels',[],'chantype','eeg','layout',[],'highlight',[],'showlabels',0,'xlim',[],'nrows',1,'ncols',1,'electrodes','on','style','straight')
  axis square; colorbar; caxis([0 .25]); title(sprintf('EEG delay map at %gsec',PeakTime));
  
  
  figure; plot(eegcstneg.epochs(ix).Delays,eegcstneg.epochs(ix).Distance,'.');
  title(sprintf('eeg at %gs (R=%g,p=%g)',eegcstneg.epochs(ix).HistTime,eegcstneg.epochs(ix).R,eegcstneg.epochs(ix).p));
  axis([0 .4 0 180]); xlabel('delay (s)'); ylabel('angular distance (deg)');
  
%   tmp = ts_data_selection(eegdata,'chanlabel',eegcstneg.epochs(ix).InvolvedChans);
%   [X,Y] = get_cartesian_coords(tmp,[]);
%   D     = dist([X Y]');
%   D     = D(1,:);
%   [sel1,sel2] = match_str(eegcstneg.epochs(ix).InvolvedChans,{tmp.sensor_info.label})
%   figure; plot(eegcstneg.epochs(ix).Delays,D,'.'); title('eeg'); axis([0 .4 0 .4]);
%   
  
%   Z = z;
%   Z(isnan(Z)) = 0;
%   bin = Z>0;
%   avg = [Z(g1) + Z(g2)] ./ (bin(g1) + bin(g2));
%   avg = avg';
%     avg(isnan(avg)) = 0;
%     tmp = ts_matrix2avg([avg avg],'sens',clusters.sensor_info(g1));
%     tmp.averages.data = avg;
%     tmp.averages.time = 0;
%     params = SO_params(1);
%     badchan = ismember({detections.label},params.badlabels);
%     badchan = find(badchan(g1) + badchan(g2) > 0);
% %     badchan = find(avg==0)
%     zlim = [];
%     badchan = [];
%     figure; ts_ezplot(tmp,'badchans',badchan,'highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
  
  
%% Histograms (num chans per cluster / num clusters per chan)

% Fig X1AB - topoplot detection density (EEG, MEG)
cd /space/mdkm1/2/kmdev/projects/jsherfey/sleep/meeg/s1/matfiles
% load SO_clusters_grad_allpeak_noflip
% load SO_clusters_grad_allpeak_noflip_ref-eeg_pospeak.mat
load SO_clusters_grad_allpeak_noflip_ref-eeg_negpeak.mat
megcst = clusters_Near_all(1);
load SO_clusters_eeg_pospeak_noflip
eegcstpos = clusters_Near_pks(1);
% % load SO_clustered_eeg_pospeak_detections.mat
% % eegdetpos = detections;
load SO_clusters_eeg_negpeak_noflip
eegcstneg = clusters_Near_pks(2);
% % load SO_clustered_eeg_negpeak_detections.mat
% % eegdetneg = detections;

% Positive EEG
clusters    = eegcstpos;
nchan       = length(clusters.sensor_info);
labels      = {clusters.sensor_info.label};
origins     = arrayfun(@(x)x.RefChan,clusters.epochs,'uniformoutput',false);
AllInvChn   = {clusters.epochs.InvolvedChans};
AllInvChn   = [AllInvChn{:}];
OriginCount = zeros(nchan,1);
InvolvCount = zeros(nchan,1);
tmp     = {clusters.epochs.InvolvedChans}; 
invchan = cellfun(@(x)[x{:}],tmp,'uniformoutput',false);
for k   = 1:nchan
  label = labels{k};
  OriginCount(k) = numel(find(ismember(origins,label)));
  InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
end
ClusterSize = arrayfun(@(x)length(x.InvolvedChans),clusters.epochs);

figure
subplot(2,2,1),try hist(InvolvCount(InvolvCount~=0),25); end; title('pos EEG'); xlabel('# clusters');
subplot(2,2,3),try hist(ClusterSize,25); end; xlabel('cluster size')

% Negative EEG
clusters    = eegcstneg;
nchan       = length(clusters.sensor_info);
labels      = {clusters.sensor_info.label};
origins     = arrayfun(@(x)x.RefChan,clusters.epochs,'uniformoutput',false);
AllInvChn   = {clusters.epochs.InvolvedChans};
AllInvChn   = [AllInvChn{:}];
OriginCount = zeros(nchan,1);
InvolvCount = zeros(nchan,1);
tmp     = {clusters.epochs.InvolvedChans}; 
invchan = cellfun(@(x)[x{:}],tmp,'uniformoutput',false);
for k   = 1:nchan
  label = labels{k};
  OriginCount(k) = numel(find(ismember(origins,label)));
  InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
end
ClusterSize = arrayfun(@(x)length(x.InvolvedChans),clusters.epochs);
subplot(2,2,2),try hist(InvolvCount(InvolvCount~=0),25); end; title('neg EEG'); xlabel('# clusters');
subplot(2,2,4),try hist(ClusterSize,25); end; xlabel('cluster size')

%  EEG



