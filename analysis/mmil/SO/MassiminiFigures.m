% From Massimini:
% - Delay maps of neg peak, pos zero-crossing, neg zero-crossing
%   - R(neg peak, pos zero) > R(neg peak, neg zero) (p<.05) in most cycles
% - (amp of following pos wave: neg-to-pos peak) vs (delay of neg peak)
%   - was parabolic in Massimini (traveling wave paper)
% - probability of electrodes being origins
%   - anterior electrodes had higher prob
%   - no individual electrode exceeded 10% prob (range of max: 6.2% to 9%)
%   - dots shown at origins in topo-type plot w/ size proportional to # cycles
%   
% Notes:
% - individual cycles affected 5-164 channels (most affected 30-50 chans)
% - max delay: 40-360ms (115+/-9.9ms)
% - for single cycles: corr(var(delay),(# chans involved)) = .75 (p<.05)
% - linear correlation was significant (p<.01) in 55% of cycles (n=20 chans)
% - constant speed: 1.2 - 7 m/s (2.7+/-.2 m/s)
% - a continuous gradient of latencies was evident in 80% of cycles
% - (# detections in stage 3) = 12.6+/-5.5; (stage 4) = 21+/-7.8
% - (# detections in entire session) = 249-468 (320+/-86)
% - In all subjects, SO inter-detection interval (IDI) = 1.25+/-.13 sec
%   => ~ .8Hz
tic
SubjID = 's1';
cd(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/%s/matfiles',SubjID));
load proc_epoch_data_ICA; data = epoch_data; clear epoch_data
load SO_detections;                       alldet  = detections;
load SO_clustered_detections;             cstdet  = detections;
load SO_clustered_consistent_detections;  condet  = detections;
load SO_clusters_noflip;
clusters_Near_all = clusters_Near_all(1);
clusters_Nmax_all = clusters_Nmax_all(1);
load SO_clusters_consistency;
clusters_Cear_all = clusters_Cear_all(1);
clusters_Cmax_all = clusters_Cmax_all(1);

%% I. Origin
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
% (# times chan is origin) / (# times chan is in cluster)
X = OriginCount ./ InvolvCount; X(isnan(X)) = 0;
% (# times chan is origin) / (total # of clusters)
Y = OriginCount ./ ntrial;      Y(isnan(Y)) = 0;
Z = InvolvCount ./ ntrial;      Z(isnan(Z)) = 0;

g1 = strmatch('grad1',{sens.typestring}); [badg1,jnk] = match_str({sens(g1).label},{sens(badchans).label});
g2 = strmatch('grad2',{sens.typestring}); [badg2,jnk] = match_str({sens(g2).label},{sens(badchans).label});

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
% (# times chan is origin) / (# times chan is in cluster)
X = OriginCount ./ InvolvCount; X(isnan(X)) = 0;
% (# times chan is origin) / (total # of clusters)
Y = OriginCount ./ ntrial;      Y(isnan(Y)) = 0;
Z = InvolvCount ./ ntrial;      Z(isnan(Z)) = 0;

figure('name','Origin: Earliest Detection')
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
Nalldet = cellfun(@length,{alldet.negpeak}) + cellfun(@length,{alldet.pospeak}); X = Nalldet';
Ncstdet = cellfun(@length,{cstdet.negpeak}) + cellfun(@length,{cstdet.pospeak}); Y = Ncstdet';
Ncondet = cellfun(@length,{condet.negpeak}) + cellfun(@length,{condet.pospeak}); Z = Ncondet';

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


%% Histograms
figure('name','Histograms (maxR, all trials)');
% cluster size
Ccstsiz   = cellfun(@length,{clusters_Cmax_all.epochs.InvolvedChans});
Ncstsiz   = cellfun(@length,{clusters_Nmax_all.epochs.InvolvedChans});
[Nc,Xc]   = hist(Ccstsiz,25);
[Nn,Xn]   = hist(Ncstsiz,25);
subplot(2,3,1); plot(Xc,Nc,'b-',Xn,Nn,'r-'); legend('Consistent','No flip')
xlabel('cluster size (# involved channels)'); ylabel('count'); axis tight; vline(40,'k'); vline(20,'k');

% max delay
Cmaxdel   = cellfun(@max,{clusters_Cmax_all.epochs.Delays});
Nmaxdel   = cellfun(@max,{clusters_Nmax_all.epochs.Delays});
[Nc,Xc]   = hist(Cmaxdel,25);
[Nn,Xn]   = hist(Nmaxdel,25);
subplot(2,3,2); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('max delay (s)');

% mean delay
Cmeandel  = cellfun(@mean,{clusters_Cmax_all.epochs.Delays});
Nmeandel  = cellfun(@mean,{clusters_Nmax_all.epochs.Delays});
[Nc,Xc]   = hist(Cmeandel,25);
[Nn,Xn]   = hist(Nmeandel,25);
subplot(2,3,3); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('mean delay (s)');

% inter-detection-interval
Cidi  = diff([clusters_Cmax_all.epochs.RefTime]);
Nidi  = diff([clusters_Nmax_all.epochs.RefTime]);
[Nc,Xc]   = hist(Cidi(Cidi<200),200);
[Nn,Xn]   = hist(Nidi(Nidi<200),200);
subplot(2,3,4); plot(Xc,Nc,'b-',Xn,Nn,'r-');
xlabel('inter-detection interval (s)'); set(gca,'xlim',[0 20]);

% absmax amplitude
Cmaxamp   = cellfun(@(x)max(abs(x)),{clusters_Cmax_all.epochs.DetectionAmps});
Nmaxamp   = cellfun(@(x)max(abs(x)),{clusters_Nmax_all.epochs.DetectionAmps});
[Nc,Xc]   = hist(Cmaxamp,25);
[Nn,Xn]   = hist(Nmaxamp,25);
subplot(2,3,5); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('absmax amplitude (fT/cm2)');

% mean amplitude
Cmeanamp  = cellfun(@mean,{clusters_Cmax_all.epochs.DetectionAmps});
Nmeanamp  = cellfun(@mean,{clusters_Nmax_all.epochs.DetectionAmps});
[Nc,Xc]   = hist(Cmeanamp,25);
[Nn,Xn]   = hist(Nmeanamp,25);
subplot(2,3,6); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('mean amplitude (fT/cm2)');


%% (% pos clustered detections) vs (% neg clustered detections)
nchan = length(alldet);
Nalln = cellfun(@length,{alldet.negpeak});
Nallp = cellfun(@length,{alldet.pospeak});
Ncstn = cellfun(@length,{cstdet.negpeak});
Ncstp = cellfun(@length,{cstdet.pospeak});
fracn = Ncstn ./ Nalln;
fracp = Ncstp ./ Nallp;
figure;
subplot(2,1,1),plot(fracp,fracn,'.'); axis([0 1 0 1]); lsline; axis square
ylabel('neg peaks'); title('% detections in clusters (no flipping)');
Ncstn = cellfun(@length,{condet.negpeak});
Ncstp = cellfun(@length,{condet.pospeak});
fracn = Ncstn ./ Nalln;
fracp = Ncstp ./ Nallp;
subplot(2,1,2),plot(fracp,fracn,'.'); axis([0 1 0 1]); lsline; axis square
xlabel('pos peaks'); ylabel('neg peaks'); title('% detections in clusters (consistency analysis)');


%% delay maps for consecutive cycles


%% epoching high and low R trials
labels    = {condet.label};
nchan     = length(labels);
Ncondet   = cellfun(@length,{condet.cluster_number});
badchans  = find(Ncondet==0);
sens      = clusters_Cmax_all.sensor_info;

alpha   = .05;

Fs      = data.sfreq;
t       = data.epochs.time;
nsamp   = length(t);
maxtime = t(end);
prestim   = 1; % sec
pad       = round(prestim*Fs);

% all trials, p<alpha, R low
Rthresh1  = 0;
Rthresh2  = .35;
origins   = arrayfun(@(x)x.RefChan,clusters_Cmax_all.epochs,'uniformoutput',false);
RefTime   = arrayfun(@(x)x.RefTime,clusters_Cmax_all.epochs);
RefDel    = arrayfun(@(x)x.RefTimeFromFirst,clusters_Cmax_all.epochs);
ntrial    = length(origins);
FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clusters_Cmax_all.epochs);
MaxDelay  = arrayfun(@(x)max(x.Delays),clusters_Cmax_all.epochs);
Rvals     = arrayfun(@(x)max(x.R),clusters_Cmax_all.epochs);
pvals     = arrayfun(@(x)max(x.p),clusters_Cmax_all.epochs);
Nvals     = arrayfun(@(x)max(x.N),clusters_Cmax_all.epochs);
Rix     = find(pvals<alpha & Rvals<Rthresh2 & Rvals>Rthresh1);%find(pvals<alpha & Rvals<Rthresh);
Rvals   = Rvals(Rix);
pvals   = pvals(Rix);
Nvals   = Nvals(Rix);
RefTime = RefTime(Rix);
trig    = [clusters_Cmax_all.epochs(Rix).RefTime];
trig    = trig(trig<maxtime);
samp    = cellfun(@(x)nearest(t,x),num2cell(trig));
keep      = ~((samp-pad)<1 | (samp+pad)>nsamp); % not outbounds
samp      = samp(keep);
ntrl      = length(samp);
begsample = round(samp - pad);
endsample = round(samp + pad);
s0        = num2cell(begsample);
sf        = num2cell(endsample);
T         = cellfun(@(x,y)(data.epochs.time(x:y)),s0,sf,'UniformOutput',false);
T         = [T{:}];
tmp       = cellfun(@(x,y)(data.epochs.data(:,x:y)),s0,sf,'UniformOutput',false);
tmp = cat(3,tmp{:});
tmp = reshape(tmp,[size(tmp,1) size(tmp,2)*size(tmp,3)]);
tmpdatlo  = ts_matrix2epoch(tmp,'continuous',1,'time',[0:size(tmp,2)]/Fs,'sens',sens);
% visualizer(tmpdatlo);

% all trials, p<alpha, R hi
Rthresh1  = .7;
Rthresh2  = 1;
origins   = arrayfun(@(x)x.RefChan,clusters_Cmax_all.epochs,'uniformoutput',false);
RefTime   = arrayfun(@(x)x.RefTime,clusters_Cmax_all.epochs);
RefDel    = arrayfun(@(x)x.RefTimeFromFirst,clusters_Cmax_all.epochs);
ntrial    = length(origins);
FirstDet  = arrayfun(@(x)min(x.DetectionTimes),clusters_Cmax_all.epochs);
MaxDelay  = arrayfun(@(x)max(x.Delays),clusters_Cmax_all.epochs);
Rvals     = arrayfun(@(x)max(x.R),clusters_Cmax_all.epochs);
pvals     = arrayfun(@(x)max(x.p),clusters_Cmax_all.epochs);
Nvals     = arrayfun(@(x)max(x.N),clusters_Cmax_all.epochs);
Rix     = find(pvals<alpha & Rvals<Rthresh2 & Rvals>Rthresh1);
Rvals   = Rvals(Rix);
pvals   = pvals(Rix);
Nvals   = Nvals(Rix);
RefTime = RefTime(Rix);
trig    = [clusters_Cmax_all.epochs(Rix).RefTime];
trig    = trig(trig<maxtime);
samp    = cellfun(@(x)nearest(t,x),num2cell(trig));
keep      = ~((samp-pad)<1 | (samp+pad)>nsamp); % not outbounds
samp      = samp(keep);
ntrl      = length(samp);
begsample = round(samp - pad);
endsample = round(samp + pad);
s0        = num2cell(begsample);
sf        = num2cell(endsample);
T         = cellfun(@(x,y)(data.epochs.time(x:y)),s0,sf,'UniformOutput',false);
T         = [T{:}];
tmp       = cellfun(@(x,y)(data.epochs.data(:,x:y)),s0,sf,'UniformOutput',false);
tmp = cat(3,tmp{:});
tmp = reshape(tmp,[size(tmp,1) size(tmp,2)*size(tmp,3)]);
tmpdathi = ts_matrix2epoch(tmp,'continuous',1,'time',[0:size(tmp,2)]/Fs,'sens',sens);
% visualizer(tmpdathi);

if 0
  figure;
  subplot(1,3,1),try hist(Rvals,100); end; title('R (all)'); set(gca,'xlim',[-1 1]);
  subplot(1,3,2),try hist(pvals,100); end; title('p'); set(gca,'xlim',[0 1]);
  subplot(1,3,3),try hist(Nvals,100); end; title('N'); axis tight
  figure;
  subplot(1,3,1),try hist(Rvals1,100); end; title('R (1)'); set(gca,'xlim',[-1 1]);
  subplot(1,3,2),try hist(pvals1,100); end; title('p'); set(gca,'xlim',[0 1]);
  subplot(1,3,3),try hist(Nvals1,100); end; title('N'); axis tight
  figure;
  subplot(1,3,1),try hist(Rvals2,100); end; title('R (2)'); set(gca,'xlim',[-1 1]);
  subplot(1,3,2),try hist(pvals2,100); end; title('p'); set(gca,'xlim',[0 1]);
  subplot(1,3,3),try hist(Nvals2,100); end; title('N'); axis tight
end

% Rvals1    = arrayfun(@(x)max(x.R),clusters_Cmax_pks(1).epochs,'uniformoutput',false);
% Rvals2    = arrayfun(@(x)max(x.R),clusters_Cmax_pks(2).epochs,'uniformoutput',false);
% rmix      = find(cellfun(@isempty,Rvals1)); clusters_Cmax_pks(1).epochs(rmix) = [];
% rmix      = find(cellfun(@isempty,Rvals2)); clusters_Cmax_pks(2).epochs(rmix) = [];
% Rvals1    = arrayfun(@(x)max(x.R),clusters_Cmax_pks(1).epochs);
% pvals1    = arrayfun(@(x)max(x.p),clusters_Cmax_pks(1).epochs);
% Nvals1    = arrayfun(@(x)max(x.N),clusters_Cmax_pks(1).epochs);
% Rvals2    = arrayfun(@(x)max(x.R),clusters_Cmax_pks(2).epochs);
% pvals2    = arrayfun(@(x)max(x.p),clusters_Cmax_pks(2).epochs);
% Nvals2    = arrayfun(@(x)max(x.N),clusters_Cmax_pks(2).epochs);

% Rix1    = find(pvals1<alpha & Rvals1<Rthresh);
% trig1   = [clusters_Cmax_pks(1).epochs(Rix1).RefTime];
% trig1   = trig(trig1<maxtime);
% samp1   = cellfun(@(x)nearest(t,x),num2cell(trig1));
% keep      = ~((samp1-pad)<1 | (samp1+pad)>nsamp); % not outbounds
% samp1     = samp1(keep);
% ntrl      = length(samp1);
% begsample = round(samp1 - pad);
% endsample = round(samp1 + pad);
% s0        = num2cell(begsample);
% sf        = num2cell(endsample);
% tmp1       = cellfun(@(x,y)(data.epochs.data(:,x:y)),s0,sf,'UniformOutput',false);
% tmp1 = cat(3,tmp1{:});
% tmp1 = reshape(tmp1,[size(tmp1,1) size(tmp1,2)*size(tmp1,3)]);
% tmp1 = ts_matrix2epoch(tmp1,'continuous',1);
% visualizer(tmp1);
% 
% trig2   = [clusters_Cmax_pks(2).epochs(Rix2).RefTime];
% Rix2    = find(pvals2<alpha & Rvals2<Rthresh);
% trig2   = trig(trig2<maxtime);
% samp2   = cellfun(@(x)nearest(t,x),num2cell(trig2));
% keep      = ~((samp2-pad)<1 | (samp2+pad)>nsamp); % not outbounds
% samp2     = samp2(keep);
% ntrl      = length(samp2);
% begsample = round(samp2 - pad);
% endsample = round(samp2 + pad);
% s0        = num2cell(begsample);
% sf        = num2cell(endsample);
% tmp2       = cellfun(@(x,y)(data.epochs.data(:,x:y)),s0,sf,'UniformOutput',false);
% tmp2 = cat(3,tmp2{:});
% tmp2 = reshape(tmp2,[size(tmp2,1) size(tmp2,2)*size(tmp2,3)]);
% tmp2 = ts_matrix2epoch(tmp2,'continuous',1);
% visualizer(tmp2);

%% Detection
detections = cstdet;
params    = SO_params(SubjID);
t         = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
typestr   = ''; % {'','pospeak_','negpeak_'} where ''=>both
cnumfield = sprintf('%scluster_number',typestr);
cindfield = sprintf('%scluster_time_index',typestr);
clusterID = arrayfun(@(x)(x.(cnumfield)),detections,'UniformOutput',false);
tmpID     = [clusterID{:}];
clusterID = unique(tmpID); % sorted cluster index, 1:nclusters
clusterIX = arrayfun(@(x)(x.(cindfield)),detections,'UniformOutput',false);
clusterIX = [clusterIX{:}];
clusterIX = cellfun(@(x)unique(clusterIX(tmpID==x)),num2cell(unique(tmpID)));
nclusters = length(clusterID);
cnum      = {detections.(cnumfield)};
cnum      = [cnum{:}];
Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));
% create cell array listing channels involved in each cluster
cnum          = {detections.(cnumfield)};
tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
tmp           = [tmp{:}]';
tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);
% center times for each cluster
tc            = t(clusterIX);
t0f           = [IntervalT0 IntervalTf];
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  POTENTIAL FIGURE FOR PAPER
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
figure('Name',sprintf('Subject %s: Detection Count',SubjID)); 
nn    = round((max(tc)-min(tc))/params.AASM_EpochLength);
[N,X] = hist(tc,nn);
tmpx  = zeros(1,length(t));
tmpx(clusterIX) = 1;
subplot(4,1,1),plot(t,tmpx,'.-'); axis tight; title('Aggregate slow oscillation detection count')
subplot(4,1,2),try hist(tc,nn); end; axis tight
subplot(4,1,3),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
xlabel('time (sec)'); ylabel(sprintf('count (%gsec bins)',params.AASM_EpochLength));
for k=1:length(t0f), vline(t0f(k),'k'); end
subplot(4,1,4),plot(tc,'.-'); ylabel('cluster time (sec)'); xlabel('cluster number'); axis tight;
clear nn N X tmpx t0f
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Correlation coefficients
labels    = {condet.label};
nchan     = length(labels);
Ncstdet   = cellfun(@length,{cstdet.cluster_number});
Ncondet   = cellfun(@length,{condet.cluster_number});
Ncstdetp  = cellfun(@length,{cstdet.pospeak});
Ncondetp  = cellfun(@length,{condet.pospeak});
Ncstdetn  = cellfun(@length,{cstdet.negpeak});
Ncondetn  = cellfun(@length,{condet.negpeak});


alphas    = [.05 1];
for k = 1:length(alphas)
  alpha   = alphas(k);
  Rthresh = num2cell(-1:.1:1);
  ix      = [clusters_Nmax_pks(1).epochs.p] < alpha;
  NmaxR1  = cellfun(@(x)length(find([clusters_Nmax_pks(1).epochs(ix).R]>x)),Rthresh);
  ix      = [clusters_Cmax_pks(1).epochs.p] < alpha;
  CmaxR1  = cellfun(@(x)length(find([clusters_Cmax_pks(1).epochs(ix).R]>x)),Rthresh);
  ix      = [clusters_Near_pks(1).epochs.p] < alpha;
  NearR1  = cellfun(@(x)length(find([clusters_Near_pks(1).epochs(ix).R]>x)),Rthresh);
  ix      = [clusters_Cear_pks(1).epochs.p] < alpha;
  CearR1  = cellfun(@(x)length(find([clusters_Cear_pks(1).epochs(ix).R]>x)),Rthresh);

  ix      = [clusters_Nmax_pks(2).epochs.p] < alpha;
  NmaxR2  = cellfun(@(x)length(find([clusters_Nmax_pks(2).epochs(ix).R]>x)),Rthresh);
  ix      = [clusters_Cmax_pks(2).epochs.p] < alpha;
  CmaxR2  = cellfun(@(x)length(find([clusters_Cmax_pks(2).epochs(ix).R]>x)),Rthresh);
  ix      = [clusters_Near_pks(2).epochs.p] < alpha;
  NearR2  = cellfun(@(x)length(find([clusters_Near_pks(2).epochs(ix).R]>x)),Rthresh);
  ix      = [clusters_Cear_pks(2).epochs.p] < alpha;
  CearR2  = cellfun(@(x)length(find([clusters_Cear_pks(2).epochs(ix).R]>x)),Rthresh);

  ix     = [clusters_Nmax_all.epochs.p] < alpha;
  NmaxR  = cellfun(@(x)length(find([clusters_Nmax_all.epochs(ix).R]>x)),Rthresh);
  ix     = [clusters_Cmax_all.epochs.p] < alpha;
  CmaxR  = cellfun(@(x)length(find([clusters_Cmax_all.epochs(ix).R]>x)),Rthresh);
  ix     = [clusters_Near_all.epochs.p] < alpha;
  NearR  = cellfun(@(x)length(find([clusters_Near_all.epochs(ix).R]>x)),Rthresh);
  ix     = [clusters_Cear_all.epochs.p] < alpha;
  CearR  = cellfun(@(x)length(find([clusters_Cear_all.epochs(ix).R]>x)),Rthresh);

  Rthresh = [Rthresh{:}];

  figure('name',sprintf('p<=%g & R>Rthresh',alpha)); lims = [0 1 0 1];
  n1=length(clusters_Nmax_pks(1).epochs); n2=length(clusters_Cmax_pks(1).epochs);
  subplot(2,3,1); plot(Rthresh,NmaxR1/n1,'b.-',Rthresh,CmaxR1/n2,'r.-'); axis(lims); vline(0,'k');
  legend('NF','CF'); ylabel('% trials'); title('pos (maxR)');
  n1=length(clusters_Nmax_pks(2).epochs); n2=length(clusters_Cmax_pks(2).epochs);
  subplot(2,3,2); plot(Rthresh,NmaxR2/n1,'b.-',Rthresh,CmaxR2/n2,'r.-'); title('neg (maxR)'); axis(lims); vline(0,'k');
  n1=length(clusters_Nmax_all.epochs);    n2=length(clusters_Cmax_all.epochs);
  subplot(2,3,3); plot(Rthresh,NmaxR/n1,'b.-',Rthresh,CmaxR/n2,'r.-'); title('all (maxR)'); axis(lims); vline(0,'k');

  n1=length(clusters_Near_pks(1).epochs); n2=length(clusters_Cear_pks(1).epochs);
  subplot(2,3,4); plot(Rthresh,NearR1/n1,'b.-',Rthresh,CearR1/n2,'r.-');  axis(lims); vline(0,'k');
  ylabel('% trials'); xlabel('Rthresh'); title('pos (earliest)');
  n1=length(clusters_Near_pks(2).epochs); n2=length(clusters_Cear_pks(2).epochs);
  subplot(2,3,5); plot(Rthresh,NearR2/n1,'b.-',Rthresh,CearR2/n2,'r.-'); title('neg (earliest)'); axis(lims); vline(0,'k');
  n1=length(clusters_Near_all.epochs);    n2=length(clusters_Cear_all.epochs);
  subplot(2,3,6); plot(Rthresh,NearR/n1,'b.-',Rthresh,CearR/n2,'r.-'); title('all (earliest)'); axis(lims); vline(0,'k');
end

toc



