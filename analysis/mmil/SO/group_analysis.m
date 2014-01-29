%% Group analysis (18-Aug-2010)
subj    = [1 2 4 5 6 8]; nsubj = length(subj);
alphas  = [.05 1];
singlepeaks = 1;

tic
for k = 1:length(alphas)
  alpha   = alphas(k);
  Rthresh = num2cell(-1:.1:1);

  nr      = length(Rthresh);
  NmaxR1  = zeros(1,nr); nmax1 = 0;
  CmaxR1  = zeros(1,nr); cmax1 = 0;
  MmaxR1  = zeros(1,nr); mmax1 = 0;
  NmaxR2  = zeros(1,nr); nmax2 = 0;
  CmaxR2  = zeros(1,nr); cmax2 = 0;
  MmaxR2  = zeros(1,nr); mmax2 = 0;
  NmaxR   = zeros(1,nr); nmax0 = 0;
  CmaxR   = zeros(1,nr); cmax0 = 0;
  MmaxR   = zeros(1,nr); mmax0 = 0;
  NearR1  = zeros(1,nr); near1 = 0;
  CearR1  = zeros(1,nr); cear1 = 0;
  MearR1  = zeros(1,nr); mear1 = 0;
  NearR2  = zeros(1,nr); near2 = 0;
  CearR2  = zeros(1,nr); cear2 = 0;
  MearR2  = zeros(1,nr); mear2 = 0;
  NearR   = zeros(1,nr); near0 = 0;
  CearR   = zeros(1,nr); cear0 = 0;
  MearR   = zeros(1,nr); mear0 = 0;
  
  for s   = 1:nsubj
    SubjID = sprintf('s%g',subj(s));
    cd(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/%s/matfiles',SubjID));
%     load proc_epoch_data_ICA; data = epoch_data; clear epoch_data
%     load SO_detections;                       alldet  = detections;
%     load SO_clustered_detections;             cstdet  = detections;
%     load SO_clustered_consistent_detections;  condet  = detections;
    if singlepeaks
      load SO_clusters_noflip_singlepeaks;
      clusters_Near_all = clusters_Near_all(1);
      clusters_Nmax_all = clusters_Nmax_all(1);
      load SO_clusters_consistency_singlepeaks;
      clusters_Mear_all = clusters_Cear_all;
      clusters_Mear_pks = clusters_Cear_pks;
      clusters_Mmax_all = clusters_Cmax_all;
      clusters_Mmax_pks = clusters_Cmax_pks;
      clusters_Mear_all = clusters_Mear_all(1);
      clusters_Mmax_all = clusters_Mmax_all(1);   
      load SO_clusters_consistency_singlepeaks;
      clusters_Cear_all = clusters_Cear_all(1);
      clusters_Cmax_all = clusters_Cmax_all(1);
    else
      load SO_clusters_noflip;
      clusters_Near_all = clusters_Near_all(1);
      clusters_Nmax_all = clusters_Nmax_all(1);
      load SO_clusters_consistency;
      clusters_Cear_all = clusters_Cear_all(1);
      clusters_Cmax_all = clusters_Cmax_all(1);
      load SO_clusters_flipmatrix.mat
      clusters_Mear_all = clusters_Mear_all(1);
      clusters_Mmax_all = clusters_Mmax_all(1);   
    end
  
%     labels    = {condet.label};
%     nchan     = length(labels);
%     Ncstdet   = cellfun(@length,{cstdet.cluster_number});
%     Ncondet   = cellfun(@length,{condet.cluster_number});
%     Ncstdetp  = cellfun(@length,{cstdet.pospeak});
%     Ncondetp  = cellfun(@length,{condet.pospeak});
%     Ncstdetn  = cellfun(@length,{cstdet.negpeak});
%     Ncondetn  = cellfun(@length,{condet.negpeak});
  
    ix      = [clusters_Nmax_pks(1).epochs.p] < alpha;
    NmaxR1  = NmaxR1 + cellfun(@(x)length(find([clusters_Nmax_pks(1).epochs(ix).R]>x)),Rthresh);
    nmax1   = nmax1 + length(clusters_Nmax_pks(1).epochs);
    ix      = [clusters_Cmax_pks(1).epochs.p] < alpha;
    CmaxR1  = CmaxR1 + cellfun(@(x)length(find([clusters_Cmax_pks(1).epochs(ix).R]>x)),Rthresh);
    cmax1   = cmax1 + length(clusters_Cmax_pks(1).epochs);
    ix      = [clusters_Mmax_pks(1).epochs.p] < alpha;
    MmaxR1  = MmaxR1 + cellfun(@(x)length(find([clusters_Mmax_pks(1).epochs(ix).R]>x)),Rthresh);
    mmax1   = mmax1 + length(clusters_Mmax_pks(1).epochs);    
    ix      = [clusters_Near_pks(1).epochs.p] < alpha;
    NearR1  = NearR1 + cellfun(@(x)length(find([clusters_Near_pks(1).epochs(ix).R]>x)),Rthresh);
    near1   = near1 + length(clusters_Near_pks(1).epochs);
    ix      = [clusters_Cear_pks(1).epochs.p] < alpha;
    CearR1  = CearR1 + cellfun(@(x)length(find([clusters_Cear_pks(1).epochs(ix).R]>x)),Rthresh);
    cear1   = cear1 + length(clusters_Cear_pks(1).epochs);
    ix      = [clusters_Mear_pks(1).epochs.p] < alpha;
    MearR1  = MearR1 + cellfun(@(x)length(find([clusters_Mear_pks(1).epochs(ix).R]>x)),Rthresh);
    mear1   = mear1 + length(clusters_Mear_pks(1).epochs);

    ix      = [clusters_Nmax_pks(2).epochs.p] < alpha;
    NmaxR2  = NmaxR2 + cellfun(@(x)length(find([clusters_Nmax_pks(2).epochs(ix).R]>x)),Rthresh);
    nmax2   = nmax2 + length(clusters_Nmax_pks(2).epochs);
    ix      = [clusters_Cmax_pks(2).epochs.p] < alpha;
    CmaxR2  = CmaxR2 + cellfun(@(x)length(find([clusters_Cmax_pks(2).epochs(ix).R]>x)),Rthresh);
    cmax2   = cmax2 + length(clusters_Cmax_pks(2).epochs);
    ix      = [clusters_Mmax_pks(2).epochs.p] < alpha;
    MmaxR2  = MmaxR2 + cellfun(@(x)length(find([clusters_Mmax_pks(2).epochs(ix).R]>x)),Rthresh);
    mmax2   = mmax2 + length(clusters_Mmax_pks(2).epochs);
    ix      = [clusters_Near_pks(2).epochs.p] < alpha;
    NearR2  = NearR2 + cellfun(@(x)length(find([clusters_Near_pks(2).epochs(ix).R]>x)),Rthresh);
    near2   = near2 + length(clusters_Near_pks(2).epochs);
    ix      = [clusters_Cear_pks(2).epochs.p] < alpha;
    CearR2  = CearR2 + cellfun(@(x)length(find([clusters_Cear_pks(2).epochs(ix).R]>x)),Rthresh);
    cear2   = cear2 + length(clusters_Cear_pks(2).epochs);
    ix      = [clusters_Mear_pks(2).epochs.p] < alpha;
    MearR2  = MearR2 + cellfun(@(x)length(find([clusters_Mear_pks(2).epochs(ix).R]>x)),Rthresh);
    mear2   = mear2 + length(clusters_Mear_pks(2).epochs);

    ix     = [clusters_Nmax_all.epochs.p] < alpha;
    NmaxR  = NmaxR  + cellfun(@(x)length(find([clusters_Nmax_all.epochs(ix).R]>x)),Rthresh);
    nmax0   = nmax0 + length(clusters_Nmax_all(1).epochs);
    ix     = [clusters_Cmax_all.epochs.p] < alpha;
    CmaxR  = CmaxR  + cellfun(@(x)length(find([clusters_Cmax_all.epochs(ix).R]>x)),Rthresh);
    cmax0   = cmax0 + length(clusters_Cmax_all(1).epochs);
    ix      = [clusters_Mmax_all.epochs.p] < alpha;
    MmaxR   = MmaxR  + cellfun(@(x)length(find([clusters_Mmax_all.epochs(ix).R]>x)),Rthresh);
    mmax0   = mmax0 + length(clusters_Mmax_all(1).epochs);
    ix     = [clusters_Near_all.epochs.p] < alpha;
    NearR  = NearR  + cellfun(@(x)length(find([clusters_Near_all.epochs(ix).R]>x)),Rthresh);
    near0   = near0 + length(clusters_Near_all(1).epochs);
    ix     = [clusters_Cear_all.epochs.p] < alpha;
    CearR  = CearR  + cellfun(@(x)length(find([clusters_Cear_all.epochs(ix).R]>x)),Rthresh);
    cear0   = cear0 + length(clusters_Cear_all(1).epochs);
    ix      = [clusters_Mear_all.epochs.p] < alpha;
    MearR   = MearR  + cellfun(@(x)length(find([clusters_Mear_all.epochs(ix).R]>x)),Rthresh);
    mear0   = mear0 + length(clusters_Mear_all(1).epochs);    
  end
  
  Rthresh = [Rthresh{:}];

  figure('name',sprintf('p<=%g & R>Rthresh',alpha)); lims = [0 1 0 1];
  n1=nmax1; n2=cmax1; n3=mmax1;
  subplot(2,3,1); plot(Rthresh,NmaxR1/n1,'b.-',Rthresh,CmaxR1/n2,'r.-',Rthresh,MmaxR1/n3,'g.-'); axis(lims); vline(0,'k');
  legend('NF','CF','MF'); ylabel('% trials'); title('pos (maxR)');
  n1=nmax2; n2=cmax2; n3=mmax2;
  subplot(2,3,2); plot(Rthresh,NmaxR2/n1,'b.-',Rthresh,CmaxR2/n2,'r.-',Rthresh,MmaxR2/n3,'g.-'); title('neg (maxR)'); axis(lims); vline(0,'k');
  n1=nmax0; n2=cmax0; n3=mmax0;
  subplot(2,3,3); plot(Rthresh,NmaxR/n1,'b.-',Rthresh,CmaxR/n2,'r.-',Rthresh,MmaxR/n3,'g.-'); title('all (maxR)'); axis(lims); vline(0,'k');
  n1=near1; n2=cear1; n3=mear1;
  subplot(2,3,4); plot(Rthresh,NearR1/n1,'b.-',Rthresh,CearR1/n2,'r.-',Rthresh,MearR1/n3,'g.-');  axis(lims); vline(0,'k');
  ylabel('% trials'); xlabel('Rthresh'); title('pos (earliest)');
  n1=near2; n2=cear2; n3=mear2;
  subplot(2,3,5); plot(Rthresh,NearR2/n1,'b.-',Rthresh,CearR2/n2,'r.-',Rthresh,MearR2/n3,'g.-'); title('neg (earliest)'); axis(lims); vline(0,'k');
  n1=near0; n2=cear0; n3=mear0;
  subplot(2,3,6); plot(Rthresh,NearR/n1,'b.-',Rthresh,CearR/n2,'r.-',Rthresh,MearR/n3,'g.-'); title('all (earliest)'); axis(lims); vline(0,'k');
  toc
end


%% Origins
clear Nidi Nidi_mu Nidi_std Nidi_min Nidi_max Ninv_mu Ninv_std Ninv_min Ninv_max SizeDelayvar_R SizeDelayvar_p
clear Speed_mu Speed_std Speed_min Speed_max detrate detcnt Nalldet Ncstdet Ncondet
clear originmax1a originmax1b originmax2a originmax2b originmu1a originmu1b originmu2a originmu2b
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
  end
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
  RR           = [clusters_Cmax_all.epochs.R];
  pp           = [clusters_Cmax_all.epochs.p];
  tmp          = abs([clusters_Cmax_all.epochs.Speed]);
  tmp          = tmp(pp<.05 & RR>.5);
  tmp          = tmp(~isnan(tmp));
  tmp          = tmp(~isinf(tmp));
  Speed_mu(s)  = mean(tmp);
  Speed_std(s) = std(tmp);
  Speed_min(s) = min(tmp);
  Speed_max(s) = max(tmp);
  
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
  
  % Total # of detections in session (clusterID from clustered detections)
  detcnt(s)  = round(length(clusters_Nmax_all.epochs)/2); 
  
  % I. Origin
  % consistent, all trials, origin is max R within 50ms of earliest detection
  labels    = {condet.label};
  nchan     = length(labels);
  Ncondet   = cellfun(@length,{condet.cluster_number});
  badchans  = find(Ncondet==0);
  sens      = clusters_Cmax_all.sensor_info;

  % 1. Define as sensor w/ largest R within 50ms of earliest detection
  origins   = arrayfun(@(x)x.RefChan,clusters_Cmax_all.epochs,'uniformoutput',false);
  RefTime   = arrayfun(@(x)x.RefTime,clusters_Cmax_all.epochs);
  RefDel    = arrayfun(@(x)x.RefTimeFromFirst,clusters_Cmax_all.epochs);
  ntrial    = length(origins);
  FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clusters_Cmax_all.epochs);
  MaxDelay  = arrayfun(@(x)max(x.Delays),clusters_Cmax_all.epochs);
  Rvals     = [clusters_Cmax_all.epochs.R];%arrayfun(@(x)max(x.R),clusters_Cmax_all.epochs);
  pvals     = [clusters_Cmax_all.epochs.p];%arrayfun(@(x)max(x.p),clusters_Cmax_all.epochs);
  Nvals     = [clusters_Cmax_all.epochs.N];%arrayfun(@(x)max(x.N),clusters_Cmax_all.epochs);

  AllInvChn   = {clusters_Cmax_all.epochs.InvolvedChans};
  AllInvChn   = [AllInvChn{:}];
  OriginCount = zeros(nchan,1);
  InvolvCount = zeros(nchan,1);
  for k   = 1:nchan
    label = labels{k};
    OriginCount(k) = numel(find(ismember(origins,label)));
    InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
  end

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
  origins   = arrayfun(@(x)x.RefChan,clusters_Cear_all.epochs,'uniformoutput',false);
  RefTime   = arrayfun(@(x)x.RefTime,clusters_Cear_all.epochs);
  ntrial    = length(origins);
  FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clusters_Cear_all.epochs);
  MaxDelay  = arrayfun(@(x)max(x.Delays),clusters_Cear_all.epochs);
  Rvals     = arrayfun(@(x)max(x.R),clusters_Cear_all.epochs);
  pvals     = arrayfun(@(x)max(x.p),clusters_Cear_all.epochs);
  Nvals     = arrayfun(@(x)max(x.N),clusters_Cear_all.epochs);

  AllInvChn   = {clusters_Cear_all.epochs.InvolvedChans};
  AllInvChn   = [AllInvChn{:}];
  OriginCount = zeros(nchan,1);
  InvolvCount = zeros(nchan,1);
  for k   = 1:nchan
    label = labels{k};
    OriginCount(k) = numel(find(ismember(origins,label)));
    InvolvCount(k) = numel(find(ismember(AllInvChn,label)));
  end
  
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
tmp.averages.time = 0; zlim = [0 .1];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / InvolvedCount');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 .02];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / TotalNumClusters');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 .6];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('InvolvedCount / TotalNumClusters');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar

% 1. Earliest detection origin
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

figure('name','Origin: Max R')
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 .1];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / InvolvedCount');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 .02];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('OriginCount / TotalNumClusters');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 .6];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('InvolvedCount / TotalNumClusters');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar


% 3. detection density (alldet, cstdet, condet)
X = Nalldet';
Y = Ncstdet';
Z = Ncondet';

figure('name','Detection Count')
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 .8E4];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('All detections');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 .9E3];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Clustered detections');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 .9E3];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Consistent clustered detections');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar

figure('name','Detection Count')
X   = 100*X / sum(X);
tmp = ts_matrix2avg([X X],'sens',sens); 
tmp.averages.data = X;
tmp.averages.time = 0; zlim = [0 1];
subplot(2,3,1); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('All detections');
subplot(2,3,4); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
Y   = 100*Y / sum(Y);
tmp = ts_matrix2avg([Y Y],'sens',sens); colorbar
tmp.averages.data = Y;
tmp.averages.time = 0; zlim = [0 1];
subplot(2,3,2); ts_ezplot(tmp,'chantype','grad1','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Clustered detections');
subplot(2,3,5); ts_ezplot(tmp,'chantype','grad2','highlight',[],'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
Z   = 100*Z / sum(Z);
tmp = ts_matrix2avg([Z Z],'sens',sens); colorbar
tmp.averages.data = Z;
tmp.averages.time = 0; zlim = [0 1];
subplot(2,3,3); ts_ezplot(tmp,'chantype','grad1','highlight',badg1,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar; title('Consistent clustered detections');
subplot(2,3,6); ts_ezplot(tmp,'chantype','grad2','highlight',badg2,'electrodes','on','zlim',zlim,'style','straight','topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
colorbar

toc

[tmpR,tmpP] = corrcoef([allNinv' allvardelay']);
