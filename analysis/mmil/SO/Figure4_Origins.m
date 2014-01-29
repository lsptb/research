%% Group analysis (18-Aug-2010)
clear all
subj    = [1 2 4 5 6 8]; nsubj = length(subj);
alphas  = [.05 1];
singlepeaks = 0;
tic
%% Origins
clear Nidi Nidi_mu Nidi_std Nidi_min Nidi_max Ninv_mu Ninv_std Ninv_min Ninv_max SizeDelayvar_R SizeDelayvar_p
clear Speed_mu Speed_std Speed_min Speed_max detrate detcnt Nalldet Ncstdet Ncondet
clear originmax1a originmax1b originmax2a originmax2b originmu1a originmu1b originmu2a originmu2b
clear ALLDET CSTDET CONDET GrandDelays1 GrandDelays2 ALLCSTDET
[GrandDelays1{1:204}] = deal([]);
[GrandDelays2{1:204}] = deal([]);
Nalldet = zeros(1,204);
Ncstdet = zeros(1,204);
Ncondet = zeros(1,204);
allvardelay = []; allNinv = [];
for s   = 1:nsubj
  SubjID = sprintf('s%g',subj(s));
  cd(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/%s/matfiles',SubjID));
  if singlepeaks
    % load proc_epoch_data_ICA; data = epoch_data; clear epoch_data
    load SO_detections_singlepeaks;                       alldet  = detections;
    load SO_clustered_detections_singlepeaks;             cstdet  = detections;
    load SO_clustered_consistent_detections_singlepeaks;  condet  = detections;
    load SO_clusters_noflip_singlepeaks;
    clusters_Near_all = clusters_Near_all(1);
    clusters_Nmax_all = clusters_Nmax_all(1);
    load SO_clusters_consistency_singlepeaks;
    clusters_Cear_all = clusters_Cear_all(1);
    clusters_Cmax_all = clusters_Cmax_all(1);
  else
    % load proc_epoch_data_ICA; data = epoch_data; clear epoch_data
    load SO_detections;                       alldet  = detections;
    load SO_clustered_detections;             cstdet  = detections;
    load SO_clustered_consistent_detections;  condet  = detections;
    load SO_clusters_noflip;
    clusters_Near_all = clusters_Near_all(1);
    clusters_Nmax_all = clusters_Nmax_all(1);
    load SO_clusters_consistency;
    clusters_Cear_all = clusters_Cear_all(1);
    clusters_Cmax_all = clusters_Cmax_all(1);
%     load SO_clustered_flipmatrix_detections;  matdet  = detections;
%     clusters_Mear_all = clusters_Mear_all(1);
%     clusters_Mmax_all = clusters_Mmax_all(1);
  end
  clusters    = clusters_Nmax_all;  % clusters_Cmax_all % clusters_Mmax_all
  clustersear = clusters_Near_all;  % clusters_Cear_all % clusters_Mear_all
  det         = cstdet;             % condet % matdet
  
  % IDI
  Nidi        = diff([clusters_Nmax_all.epochs.RefTime]);
  Nidi_mu(s)  = mean(Nidi);
  Nidi_std(s) = std(Nidi);
  Nidi_min(s) = min(Nidi);
  Nidi_max(s) = max(Nidi);
  
  % # involved channels (cluster size)
  Ncstsiz     = cellfun(@length,{clusters_Nmax_all.epochs.InvolvedChans});
  Ninv_mu(s)  = mean(Ncstsiz);
  Ninv_std(s) = std(Ncstsiz);
  Ninv_min(s) = min(Ncstsiz);
  Ninv_max(s) = max(Ncstsiz);

  % correlation bw var(delay) and cluster size
  vardelay    = arrayfun(@(x)var(x.Delays),clusters_Nmax_all.epochs);
  [tmpR,tmpP] = corrcoef([Ncstsiz' vardelay']);
  SizeDelayvar_R(s) = tmpR(1,2);
  SizeDelayvar_p(s) = tmpP(1,2);
  allNinv           = [allNinv Ncstsiz];
  allvardelay       = [allvardelay vardelay];
  
  % Speed
  RR           = [clusters.epochs.R];
  pp           = [clusters.epochs.p];
  tmp          = abs([clusters.epochs.Speed]);
  tmp          = tmp(pp<.05 & RR>.5);
  tmp          = tmp(~isnan(tmp));
  tmp          = tmp(~isinf(tmp));
  Speed_mu(s)  = mean(tmp);
  Speed_std(s) = std(tmp);
  Speed_min(s) = min(tmp);
  Speed_max(s) = max(tmp);
%   1,toc
  % # detections per minute (use cstdet)
  T     = 0;
  nn    = 0;
  tref  = [clusters_Nmax_all.epochs.RefTime];
  for j = 1:length(IntervalT0)
    t0  = IntervalT0(j);
    tf  = IntervalTf(j);
    nn  = length(find(tref>=t0&tref<=tf)) + nn;
    T   = T + (tf-t0);
  end
  detrate(s) = round(60*nn/T/2); % cycles per minute
%   2,toc
  % Total # of detections in session (clusterID from clustered detections)
  detcnt(s)  = round(length(clusters_Nmax_all.epochs)/2); 
  
  % I. Origin
  % consistent or noflip, all trials, origin is max R within 50ms of earliest detection
  labels    = {det.label};
  nchan     = length(labels);
  Ncondet   = cellfun(@length,{det.cluster_number});
  badchans  = find(Ncondet==0);
  sens      = clusters.sensor_info;

  % 1. Define as sensor w/ largest R within 50ms of earliest detection
  origins   = arrayfun(@(x)x.RefChan,clusters.epochs,'uniformoutput',false);
  RefTime   = arrayfun(@(x)x.RefTime,clusters.epochs);
  RefDel    = arrayfun(@(x)x.RefTimeFromFirst,clusters.epochs);
  ntrial    = length(origins);
  FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clusters.epochs);
  MaxDelay  = arrayfun(@(x)max(x.Delays),clusters.epochs);
  Rvals     = [clusters.epochs.R];%arrayfun(@(x)max(x.R),clusters_Cmax_all.epochs);
  pvals     = [clusters.epochs.p];%arrayfun(@(x)max(x.p),clusters_Cmax_all.epochs);
  Nvals     = [clusters.epochs.N];%arrayfun(@(x)max(x.N),clusters_Cmax_all.epochs);
%   3,toc
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
%     ix = find(cellfun(@(x)ismember(label,x),{clusters.epochs.InvolvedChans}));
    ix = find(~cellfun(@isempty,regexp(invchan,label)));
%     if s==1, GrandDelays1{k} = []; end
    GrandDelays1{k} = [GrandDelays1{k} [clusters.epochs(ix).Delays]];
  end
%   4,toc
  if s==1
    GrandOriginCount1 = OriginCount;
    GrandInvolvCount1 = InvolvCount;
    Grandntrial1      = ntrial;
  else
    GrandOriginCount1 = GrandOriginCount1 + OriginCount;
    GrandInvolvCount1 = GrandInvolvCount1 + InvolvCount;
    Grandntrial1      = Grandntrial1 + ntrial;
  end
  tmp             = OriginCount ./ InvolvCount; tmp(isnan(tmp)) = 0;
  originmax1a(s)  = max(tmp); originmu1a(s) = mean(tmp(tmp~=0));
  tmp             = OriginCount ./ ntrial; tmp(isnan(tmp)) = 0;
  originmax1b(s)  = max(tmp); originmu1b(s) = mean(tmp(tmp~=0));
  
  % 2. Define as earliest detection
  origins   = arrayfun(@(x)x.RefChan,clustersear.epochs,'uniformoutput',false);
  RefTime   = arrayfun(@(x)x.RefTime,clustersear.epochs);
  ntrial    = length(origins);
  FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clustersear.epochs);
  MaxDelay  = arrayfun(@(x)max(x.Delays),clustersear.epochs);
  Rvals     = [clustersear.epochs.R];%arrayfun(@(x)max(x.R),clustersear.epochs);
  pvals     = [clustersear.epochs.p];%arrayfun(@(x)max(x.p),clustersear.epochs);
  Nvals     = [clustersear.epochs.N];%arrayfun(@(x)max(x.N),clustersear.epochs);

  AllInvChn   = {clustersear.epochs.InvolvedChans};
  AllInvChn   = [AllInvChn{:}];
  OriginCount = zeros(nchan,1);
  InvolvCount = zeros(nchan,1);
  tmp     = {clustersear.epochs.InvolvedChans}; 
  invchan = cellfun(@(x)[x{:}],tmp,'uniformoutput',false);
  for k   = 1:nchan
    label = labels{k};
    OriginCount(k) = numel(find(ismember(origins,label)));
    InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
%     ix = find(cellfun(@(x)ismember(label,x),{clustersear.epochs.InvolvedChans}));
    ix = find(~cellfun(@isempty,regexp(invchan,label)));
%     if s==1, GrandDelays2{k} = []; end
    GrandDelays2{k} = [GrandDelays2{k} [clustersear.epochs(ix).Delays]];    
  end
%   5,toc
  if s==1
    GrandOriginCount2 = OriginCount;
    GrandInvolvCount2 = InvolvCount;
    Grandntrial2      = ntrial;
  else
    GrandOriginCount2 = GrandOriginCount2 + OriginCount;
    GrandInvolvCount2 = GrandInvolvCount2 + InvolvCount;
    Grandntrial2      = Grandntrial2 + ntrial;
  end
  tmp             = OriginCount ./ InvolvCount; tmp(isnan(tmp)) = 0;
  originmax2a(s)  = max(tmp); originmu2a(s) = mean(tmp(tmp~=0));
  tmp             = OriginCount ./ ntrial; tmp(isnan(tmp)) = 0;
  originmax2b(s)  = max(tmp); originmu2b(s) = mean(tmp(tmp~=0));

  % 3. detection density (alldet, cstdet, condet)
  Nalldet = Nalldet + cellfun(@length,{alldet.negpeak}) + cellfun(@length,{alldet.pospeak});
  Ncstdet = Ncstdet + cellfun(@length,{cstdet.negpeak}) + cellfun(@length,{cstdet.pospeak});
  Ncondet = Ncondet + cellfun(@length,{condet.negpeak}) + cellfun(@length,{condet.pospeak});
  if s==1
    ALLDET  = alldet;
    CSTDET  = cstdet;
    CONDET  = condet;
  else
    tmp     = cellfun(@(x,y)[x y],{ALLDET.pospeak},{alldet.pospeak},'uniformoutput',false);
    [ALLDET.pospeak] = deal(tmp{:});
    tmp     = cellfun(@(x,y)[x y],{ALLDET.negpeak},{alldet.negpeak},'uniformoutput',false);
    [ALLDET.negpeak] = deal(tmp{:});
    tmp     = cellfun(@(x,y)[x y],{CSTDET.pospeak},{cstdet.pospeak},'uniformoutput',false);
    [CSTDET.pospeak] = deal(tmp{:});
    tmp     = cellfun(@(x,y)[x y],{CSTDET.negpeak},{cstdet.negpeak},'uniformoutput',false);
    [CSTDET.negpeak] = deal(tmp{:});
    tmp     = cellfun(@(x,y)[x y],{CONDET.pospeak},{condet.pospeak},'uniformoutput',false);
    [CONDET.pospeak] = deal(tmp{:});
    tmp     = cellfun(@(x,y)[x y],{CONDET.negpeak},{condet.negpeak},'uniformoutput',false);
    [CONDET.negpeak] = deal(tmp{:});
    clear tmp
  end  
%   6,toc
  [detections,count] = SO_cluster_detections(alldet,'method',params.cluster_method,'thresh',params.cluster_thresh,...
    'StepSize',params.cluster_StepSize,'IntegrationWindow',params.cluster_IntegrationWindow,'MinSeparation',params.cluster_MinSeparation,...
    'ClusterWindow',params.cluster_ClusterWindow);
  if s==1
    ALLCSTDET = detections;
  else
    tmp = cellfun(@(x,y)[x y],{ALLCSTDET.negpeak},{detections.negpeak},'uniformoutput',false);
    [ALLCSTDET.negpeak] = deal(tmp{:});
    tmp = cellfun(@(x,y)[x y],{ALLCSTDET.pospeak},{detections.pospeak},'uniformoutput',false);
    [ALLCSTDET.pospeak] = deal(tmp{:});  
    n0  = {ALLCSTDET.cluster_number}; n0 = max([n0{:}]) + 1;
    tmp = cellfun(@(x,y)[x y+n0],{ALLCSTDET.cluster_number},{detections.cluster_number},'uniformoutput',false);
    [ALLCSTDET.cluster_number] = deal(tmp{:});
    tmp = cellfun(@(x,y)[x y+n0],{ALLCSTDET.pospeak_cluster_number},{detections.pospeak_cluster_number},'uniformoutput',false);
    [ALLCSTDET.pospeak_cluster_number] = deal(tmp{:});
    tmp = cellfun(@(x,y)[x y+n0],{ALLCSTDET.negpeak_cluster_number},{detections.negpeak_cluster_number},'uniformoutput',false);
    [ALLCSTDET.negpeak_cluster_number] = deal(tmp{:});
  end
  clear detections tmp
  toc
end

% 1. Max R origin
OriginCount = GrandOriginCount1;
InvolvCount = GrandInvolvCount1;
ntrial      = Grandntrial1;
    
% (# times chan is origin) / (# times chan is in cluster)
X = OriginCount ./ InvolvCount; X(isnan(X)) = 0;
% (# times chan is origin) / (total # of clusters)
Y = OriginCount ./ ntrial;      Y(isnan(Y)) = 0;
Z = InvolvCount ./ ntrial;      Z(isnan(Z)) = 0;

badg1=[];%g1 = strmatch('grad1',{sens.typestring}); [badg1,jnk] = match_str({sens(g1).label},{sens(badchans).label});
badg2=[];%g2 = strmatch('grad2',{sens.typestring}); [badg2,jnk] = match_str({sens(g2).label},{sens(badchans).label});

figure('name','Origin: Max R')
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 .06];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / InvolvedCount');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 .01];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / TotalNumClusters');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 .25];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('InvolvedCount / TotalNumClusters');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar
set(gcf,'color','w')

% 2. Earliest detection origin
OriginCount = GrandOriginCount2;
InvolvCount = GrandInvolvCount2;
ntrial      = Grandntrial2;
    
% (# times chan is origin) / (# times chan is in cluster)
X = OriginCount ./ InvolvCount; X(isnan(X)) = 0;
% (# times chan is origin) / (total # of clusters)
Y = OriginCount ./ ntrial;      Y(isnan(Y)) = 0;
Z = InvolvCount ./ ntrial;      Z(isnan(Z)) = 0;

badg1=[];%g1 = strmatch('grad1',{sens.typestring}); [badg1,jnk] = match_str({sens(g1).label},{sens(badchans).label});
badg2=[];%g2 = strmatch('grad2',{sens.typestring}); [badg2,jnk] = match_str({sens(g2).label},{sens(badchans).label});

figure('name','Origin: Earliest detection')
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 .04];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / InvolvedCount');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 .008];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / TotalNumClusters');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 .3];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('InvolvedCount / TotalNumClusters');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar
set(gcf,'color','w')

% 3. detection density (alldet, cstdet, condet)
X = Nalldet';
Y = Ncstdet';
Z = Ncondet';

figure('name','Detection Count','color','w')
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 3E4];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('All detections');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 1E4];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Clustered detections');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 3E3];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Consistent clustered detections');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar

figure('name','Detection Count / (total count / nchan)','color','w')
X   = X / (sum(X)/204);
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 1.25];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('All detections');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
Y   = Y / (sum(Y)/204);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 1.25];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Clustered detections');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
Z   = Z / (sum(Z)/204);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 1.25];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Consistent clustered detections');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar

toc

[tmpR,tmpP] = corrcoef([allNinv' allvardelay']);

% % cluster size
% Ccstsiz   = cellfun(@length,{clusters_Cmax_all.epochs.InvolvedChans});
% Ncstsiz   = cellfun(@length,{clusters_Nmax_all.epochs.InvolvedChans});
% [Nc,Xc]   = hist(Ccstsiz,25);
% [Nn,Xn]   = hist(Ncstsiz,25);
% figure; plot(Xc,Nc,'b-',Xn,Nn,'r-'); legend('Consistent','No flip')
% xlabel('cluster size (# involved channels)'); ylabel('count'); axis tight; vline(40,'k'); vline(20,'k');

% Detection density
% [detections,count] = SO_cluster_detections(ALLDET,'method',params.cluster_method,'thresh',params.cluster_thresh,...
%   'StepSize',params.cluster_StepSize,'IntegrationWindow',params.cluster_IntegrationWindow,'MinSeparation',params.cluster_MinSeparation,...
%   'ClusterWindow',params.cluster_ClusterWindow);
detections = ALLCSTDET;

clusterID = arrayfun(@(x)(x.cluster_number),detections,'UniformOutput',false);
tmpID     = [clusterID{:}];
clusterID = unique(tmpID); % sorted cluster index, 1:nclusters
nclusters = length(clusterID);
numcst    = cellfun(@length,{detections.cluster_number});
[numdetN,numdetX]   = hist(numcst,25);
cnum      = {detections.cluster_number};
cnum      = [cnum{:}];
Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));
[cnt,ninv]= hist(Ninvolved,40);
% figure,plot(ninv,cnt,'b-');
% figure; try hist(Ninvolved,40); end

figure('name','Grand detection density','color','w');
X   = [100*numcst / nclusters]';
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 25];
subplot(3,1,1); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Detection density');
subplot(3,1,2); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar
subplot(3,1,3); try hist(Ninvolved,40); end
xlabel('Number of Channels'); ylabel('Number of SO Cycles');
axis tight; %set(gca,'xlim',[0 120]); 

X   = [100*numcst / nclusters]';
N   = cellfun(@length,GrandDelays1)';
g1  = strmatch('grad1',{sens.typestring});
g2  = strmatch('grad2',{sens.typestring});
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 25];
tmp1 = ts_data_selection(tmp,'chantype','grad1');
tmp2 = ts_data_selection(tmp,'chantype','grad2');
tmp  = tmp1;
tmp.averages.data = mean([X(g1) X(g2)],2);
% tmp.averages.data = tmp1.averages.data + tmp2.averages.data;
figure('name','Grand detection density (both grads)','color','w');
subplot(2,1,1); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Detection density (all grads)');
subplot(2,1,2); try hist(Ninvolved,40); end
xlabel('Number of Channels'); ylabel('Number of SO Cycles');
axis tight; %set(gca,'xlim',[0 120])

figure('name','Grand average delay maps','color','w'); zlim = [];
X   = cellfun(@mean,GrandDelays1)';
N   = cellfun(@length,GrandDelays1)';
g1  = strmatch('grad1',{sens.typestring});
g2  = strmatch('grad2',{sens.typestring});
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0;
tmp1 = ts_data_selection(tmp,'chantype','grad1');
tmp2 = ts_data_selection(tmp,'chantype','grad2');
tmp  = tmp1;
tmp.averages.data = (N(g1).*X(g1) + N(g2).*X(g2)) ./ (N(g1)+N(g2));
subplot(4,2,1); ts_ezplot(tmp1,'highlight',[],'electrodes','on','zlim',zlim,'style','both','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); colorbar
title('MaxR Delay Map');
subplot(4,2,3); ts_ezplot(tmp2,'highlight',[],'electrodes','on','zlim',zlim,'style','both','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); colorbar
subplot(4,2,5); ts_ezplot(tmp,'highlight',[],'electrodes','on','zlim',zlim,'style','both','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); colorbar
title('all grads');
subplot(4,2,7); try hist([GrandDelays1{:}],20); end; axis tight; xlabel('delay (s)'); ylabel('count');
X   = cellfun(@mean,GrandDelays2)';
N   = cellfun(@length,GrandDelays2)';
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0;
tmp1 = ts_data_selection(tmp,'chantype','grad1');
tmp2 = ts_data_selection(tmp,'chantype','grad2');
tmp  = tmp1;
tmp.averages.data = (N(g1).*X(g1) + N(g2).*X(g2)) ./ (N(g1)+N(g2));
subplot(4,2,2); ts_ezplot(tmp1,'highlight',[],'electrodes','on','zlim',zlim,'style','both','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); colorbar
title('Earliest Detection Delay Map');
subplot(4,2,4); ts_ezplot(tmp2,'highlight',[],'electrodes','on','zlim',zlim,'style','both','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); colorbar
subplot(4,2,6); ts_ezplot(tmp,'highlight',[],'electrodes','on','zlim',zlim,'style','both','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0); colorbar
title('all grads');
subplot(4,2,8); try hist([GrandDelays2{:}],20); end; axis tight; xlabel('delay (s)')



