addpath /space/emc2/1/halgdev/projects/sleep/MEG/SO/scripts
datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run1/s1/matfiles';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s1/matfiles';
datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s1/matfiles/proc_epoch_data_ICA.mat';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s2/matfiles/tests/2'; % w=200ms, thresh=meanstd?
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s2/matfiles/tests/3'; % w=400ms, thresh=meanstd
% datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s2/matfiles/proc_epoch_data_ICA.mat';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s4/matfiles/clusterRun2'; % w=200ms, thresh=meanstd?
% datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s4/matfiles/proc_epoch_data_ICA.mat';

% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s1/matfiles'; % w=200ms, thresh=meanstd3
% datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s1/matfiles/proc_epoch_data_ICA.mat';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s2/matfiles'; % w=200ms, thresh=meanstd3
% datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s2/matfiles/proc_epoch_data_ICA.mat';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s5/matfiles'; % w=200ms, thresh=meanstd3
datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s5/matfiles/proc_epoch_data_ICA.mat';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s6/matfiles'; % w=200ms, thresh=meanstd3
% datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s6/matfiles/proc_epoch_data_ICA.mat';
% datapath = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s8/matfiles'; % w=200ms, thresh=meanstd3
% datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s8/matfiles/proc_epoch_data_ICA.mat';

% delay comparisons
% cd /space/emc2/1/halgdev/projects/sleep/MEG/SO/s4/matfiles/clusterRun2
% load SO_detections.mat
% load SO_clustered_detections.mat
% load SO_clustered_flipmatrix_detections.mat
% load SO_clustered_consistent_detections.mat
% load SO_clusters_noflip.mat
% load SO_clusters_histflip.mat
% load SO_clusters_flipmatrix.mat
% load SO_clusters_consistency.mat
% load SO_corrcoef_R-1-1_p10.mat
% load SO_corrcoef_R-1-1_p05.mat
% cd /space/emc2/1/halgdev/projects/sleep/MEG/SO/s2/matfiles/tests/2
% load SO_detections.mat
% load SO_clustered_detections.mat
% load SO_clusters_noflip.mat
% load SO_clustered_consistent_detections.mat
% load SO_clusters_consistency.mat

[a,b]   = fileparts(datapath);
SubjID  = a(end-1:end);

%% inspect flipped data
cd(datapath);
load(datafile); % epoch_data
load SO_clustered_consistent_detections.mat
RiThreshold = mean(Ri);%-std(Ri);
badchans = find(Ri < RiThreshold-std(Ri));
pij(:,badchans) = 0;
flipvec = mode(pij,1)';
flipdat = ts_data_selection(epoch_data,'toilim',[IntervalT0(1) IntervalTf(end)]); clear epoch_data
flipmat = repmat(flipvec,1,length(flipdat.epochs.time));
flipdat.epochs.data = flipdat.epochs.data .* flipmat;
clear flipmat
visualizer(flipdat);

%% 1. Compare involvement of cond1 clusters vs cond2 clusters
tic
cd(datapath);
load SO_detections.mat;           alldet = detections;
load SO_clustered_detections.mat; cstdet = detections;

% Calc: (# clustered detections in sensor i) / (total # of detections in i)
nchan = length(alldet);
Nalln = cellfun(@length,{alldet.negpeak});
Nallp = cellfun(@length,{alldet.pospeak});
Ncstn = cellfun(@length,{cstdet.negpeak});
Ncstp = cellfun(@length,{cstdet.pospeak});
fracn = Ncstn ./ Nalln;
fracp = Ncstp ./ Nallp;

figure('name',SubjID); fig(1).name = SubjID;
subplot(2,1,1),plot(fracp,fracn,'.'); axis([0 1 0 1]); lsline; axis square
ylabel('neg peaks'); title('% detections in clusters (no flipping)');
fig(1).subplot(1).x = fracp; fig(1).subplot(1).y = fracn;
fig(1).subplot(1).title = '% detections in clusters (no flipping)';

load SO_clustered_consistent_detections.mat; cstdet = detections;
Ncstn = cellfun(@length,{cstdet.negpeak});
Ncstp = cellfun(@length,{cstdet.pospeak});
fracn = Ncstn ./ Nalln;
fracp = Ncstp ./ Nallp;

subplot(2,1,2),plot(fracp,fracn,'.'); axis([0 1 0 1]); lsline; axis square
xlabel('pos peaks'); ylabel('neg peaks'); title('% detections in clusters (consistency analysis)');
fig(1).subplot(2).x = fracp; fig(1).subplot(2).y = fracn;
fig(1).subplot(2).title = '% detections in clusters (consistency analysis)';
toc

%% 2. Compare onset std of cond1 clusters vs cond2 clusters
tic
cd(datapath);
load SO_clusters_noflip.mat;      clust1 = clusters1_1;
load SO_clusters_consistency.mat; clust2 = clusters4_1;

% Calc: std of onsets
clusters      = clust1;
sens          = clusters(1).sensor_info;
nchan         = length(sens);
recstd1       = zeros(nchan,2);
delay{1}      = {clusters(1).epochs.Delays};
delay{2}      = {clusters(2).epochs.Delays};
for ch  = 1:nchan
  this  = sens(ch).label;
  for c = 1:length(clusters)
    cstind        = find(arrayfun(@(x)(ismember(this,x.InvolvedChans)),clusters(c).epochs));
    if isempty(cstind), continue; end
    chind         = arrayfun(@(x)find(ismember(x.InvolvedChans,this)),clusters(c).epochs(cstind));
    recstd1(ch,c) = std(cellfun(@(x,y)(x(y)),delay{c}(cstind),num2cell(chind)));
  end
end
clear clusters delay cstind chind
clusters      = clust2;
sens          = clusters(1).sensor_info;
nchan         = length(sens);
recstd2       = zeros(nchan,2);
delay{1}      = {clusters(1).epochs.Delays};
delay{2}      = {clusters(2).epochs.Delays};
for ch  = 1:nchan
  this  = sens(ch).label;
  for c = 1:length(clusters)
    cstind        = find(arrayfun(@(x)(ismember(this,x.InvolvedChans)),clusters(c).epochs));
    if isempty(cstind), continue; end
    chind         = arrayfun(@(x)find(ismember(x.InvolvedChans,this)),clusters(c).epochs(cstind));
    recstd2(ch,c) = std(cellfun(@(x,y)(x(y)),delay{c}(cstind),num2cell(chind)));
  end
end
clear clusters delay cstind chind

figure('name',SubjID); fig(2).name = SubjID;
lim = max([1.5*max(recstd1(:)) 1.5*max(recstd2(:))]);
subplot(2,1,1),plot(recstd1(recstd1(:,1)~=0,1),recstd1(recstd1(:,1)~=0,2),'.'); axis([0 lim 0 lim]); lsline;
ylabel('cond2'); title('SD in clusters (no flipping)');
subplot(2,1,2),plot(recstd2(recstd2(:,1)~=0,1),recstd2(recstd2(:,1)~=0,2),'.'); axis([0 lim 0 lim]); lsline;
ylabel('cond2'); title('SD in clusters (consistency analysis)'); xlabel('cond1, [sec]');
fig(2).subplot(1).x = recstd1(recstd1(:,1)~=0,1); fig(2).subplot(1).y = recstd1(recstd1(:,1)~=0,2); fig(2).subplot(1).title = 'SD in clusters (no flipping)';
fig(2).subplot(2).x = recstd2(recstd2(:,1)~=0,1); fig(2).subplot(2).y = recstd2(recstd2(:,1)~=0,2); fig(2).subplot(2).title = 'SD in clusters (consistency analysis)';
toc

%% 3. Spatial distribution of detection frequency. Topoplot cluster count.
tic
cd(datapath);
% load SO_detections.mat;
% load SO_clustered_detections.mat;
load SO_clustered_consistent_detections.mat
Nalln = cellfun(@length,{detections.negpeak});
Nallp = cellfun(@length,{detections.pospeak});

cd(datapath);
load(datafile); % epoch_data
sens = epoch_data.sensor_info;
clear epoch_data

Nall    = Nalln + Nallp;
zlim    = [0 1];
figure('name',SubjID); fig(3).name = SubjID;
tmp = Nall / max(Nall);   avgdat  = ts_data_selection(ts_matrix2avg([tmp;tmp]','sens',sens),'toi',0);
subplot(2,3,1); ts_ezplot(avgdat,'layout',params.layout,'topoplot',1,'toprows',1,'topcols',1,'zlim',zlim,'title','all trials','newfig',0,'chantype','grad1'); fig(3).subplot(1).topodata = avgdat; fig(3).subplot(1).title = 'all trials';
subplot(2,3,4); ts_ezplot(avgdat,'layout',params.layout,'topoplot',1,'toprows',1,'topcols',1,'zlim',zlim,'title','all trials','newfig',0,'chantype','grad2'); fig(3).subplot(4).topodata = avgdat; fig(3).subplot(4).title = 'all trials';
tmp = Nallp / max(Nallp); avgdat  = ts_data_selection(ts_matrix2avg([tmp;tmp]','sens',sens),'toi',0); 
subplot(2,3,2); ts_ezplot(avgdat,'layout',params.layout,'topoplot',1,'toprows',1,'topcols',1,'zlim',zlim,'title','pos trials','newfig',0,'chantype','grad1'); fig(3).subplot(2).topodata = avgdat; fig(3).subplot(2).title = 'pos trials';
subplot(2,3,5); ts_ezplot(avgdat,'layout',params.layout,'topoplot',1,'toprows',1,'topcols',1,'zlim',zlim,'title','pos trials','newfig',0,'chantype','grad2'); fig(3).subplot(5).topodata = avgdat; fig(3).subplot(5).title = 'pos trials';
tmp = Nalln / max(Nalln); avgdat  = ts_data_selection(ts_matrix2avg([tmp;tmp]','sens',sens),'toi',0); 
subplot(2,3,3); ts_ezplot(avgdat,'layout',params.layout,'topoplot',1,'toprows',1,'topcols',1,'zlim',zlim,'title','neg trials','newfig',0,'chantype','grad1'); fig(3).subplot(3).topodata = avgdat; fig(3).subplot(3).title = 'neg trials';
subplot(2,3,6); ts_ezplot(avgdat,'layout',params.layout,'topoplot',1,'toprows',1,'topcols',1,'zlim',zlim,'title','neg trials','newfig',0,'chantype','grad2'); fig(3).subplot(6).topodata = avgdat; fig(3).subplot(6).title = 'neg trials';

toc

%% 4. Spatial pattern of involvement
tic
ROI(1).labels = {'MEG 0543','MEG 0612','MEG 0312','MEG 0342','MEG 0323','MEG 0332'}; % LF
ROI(2).labels = {'MEG 0933','MEG 1023','MEG 1213','MEG 1223','MEG 1232','MEG 1243'}; % RF (symmetrical with LF)
ROI(3).labels = {'MEG 0242','MEG 1623','MEG 1612','MEG 1642','MEG 1512','MEG 1522'}; % LB
ROI(4).labels = {'MEG 1332','MEG 2413','MEG 2422','MEG 2432','MEG 2612','MEG 2643'}; % RB (symmetrical with LB)
ROI(5).labels = {'MEG 0522','MEG 0533','MEG 0812','MEG 0823','MEG 0943','MEG 1012'}; % CF
ROI(6).labels = {'MEG 1922','MEG 1933','MEG 2112','MEG 2123','MEG 2333','MEG 2342'}; % CB

cd(datapath);
load SO_clusters_noflip.mat;      % clusters1_1, clusters1_2
load SO_clusters_consistency.mat; % clusters4_1, clusters4_2
nchan = length(clusters1_1(1).sensor_info);

clusters = clusters1_2;
delay{1} = {clusters(1).epochs.Delays};
delay{2} = {clusters(2).epochs.Delays};
clear recdelays
nr = length(ROI);
ns = length(ROI(1).labels);
for c = 1:length(clusters)
  for r = 1:nr      % loop over ROIs
    recdelays{c}{r} = [];
    for s = 1:ns    % loop over sensors
      this = ROI(r).labels{s};
      cstind = find(arrayfun(@(x)(ismember(this,x.InvolvedChans)),clusters(c).epochs));
      if isempty(cstind), continue; end
      chind  = arrayfun(@(x)find(ismember(x.InvolvedChans,this)),clusters(c).epochs(cstind));
      recdelays{c}{r} = [recdelays{c}{r} cellfun(@(x,y)(x(y)),delay{c}(cstind),num2cell(chind))];
    end
  end
end
% plot
figure('Name',sprintf('%s: no flipping',SubjID)); fig(4).name=sprintf('%s: no flipping',SubjID);
for c = 1:length(clusters)
  for r = 1:nr
    tmp = recdelays{c}{r};
    subplot(1,2,c),plot(tmp,r*ones(1,length(tmp)),'b.',mean(tmp),r,'r*');
    fig(4).subplot(c).x = tmp; fig(4).subplot(c).y = r*ones(1,length(tmp));
    hold on; clear tmp
  end
  vline(0,'k');
  axis tight
end
subplot(1,2,1),title('cond1'),xlabel('delay (s)'); subplot(1,2,2),title('cond2'),xlabel('delay (s)'); 
fig(4).subplot(1).title = 'cond1'; fig(4).subplot(2).title = 'cond2';

clusters = clusters4_2;
delay{1} = {clusters(1).epochs.Delays};
delay{2} = {clusters(2).epochs.Delays};
clear recdelays
nr = length(ROI);
ns = length(ROI(1).labels);
for c = 1:length(clusters)
  for r = 1:nr      % loop over ROIs
    recdelays{c}{r} = [];
    for s = 1:ns    % loop over sensors
      this = ROI(r).labels{s};
      cstind = find(arrayfun(@(x)(ismember(this,x.InvolvedChans)),clusters(c).epochs));
      if isempty(cstind), continue; end
      chind  = arrayfun(@(x)find(ismember(x.InvolvedChans,this)),clusters(c).epochs(cstind));
      recdelays{c}{r} = [recdelays{c}{r} cellfun(@(x,y)(x(y)),delay{c}(cstind),num2cell(chind))];
    end
  end
end
% plot
figure('Name','consistency');
cn    = length(clusters);
for c = 1:length(clusters)
  for r = 1:nr
    tmp = recdelays{c}{r};
    subplot(1,2,c),plot(tmp,r*ones(1,length(tmp)),'b.',mean(tmp),r,'r*');
    fig(4).subplot(c+cn).x = tmp; fig(4).subplot(c+cn).y = r*ones(1,length(tmp));
    hold on; clear tmp
  end
  vline(0,'k');
  axis tight
end
subplot(1,2,1),title('cond1'),xlabel('delay (s)'); subplot(1,2,2),title('cond2'),xlabel('delay (s)'); 
fig(4).subplot(3).title = 'cond1'; fig(4).subplot(4).title = 'cond2';
toc

%% 5. Co-involvement. Plot 2D image of [# co-involved clusters for channels i & j]
% Nij is the number of clusters in which i & j are co-involved
cd(datapath);
load SO_clusters_noflip.mat;      % clusters1_1, clusters1_2
load SO_clusters_consistency.mat; % clusters4_1, clusters4_2
nchan = length(clusters1_1(1).sensor_info);
sens  = clusters1_1(1).sensor_info;
tic
clusters = clusters1_1;
Nij   = zeros(nchan,nchan,2); % # co-involved clusters b/w chan i & j (pos/neg)
emp   = zeros(1,nchan);
for c = 1:length(clusters)
  for trl = 1:length(clusters(c).epochs)
    [sel1,sel2]   = match_str({sens.label},clusters(c).epochs(trl).InvolvedChans);
    mor           = emp;
    mor(1,sel1)   = 1;
    mor           = repmat(mor,[length(sel1) 1]);
    Nij(sel1,:,c) = Nij(sel1,:,c) + mor;
    clear mor
  end
end
Nij1 = Nij;
toc
clusters = clusters4_1;
Nij   = zeros(nchan,nchan,2); % # co-involved clusters b/w chan i & j (pos/neg)
emp   = zeros(1,nchan);
for c = 1:length(clusters)
  for trl = 1:length(clusters(c).epochs)
    [sel1,sel2] = match_str({sens.label},clusters(c).epochs(trl).InvolvedChans);
    mor         = emp;
    mor(1,sel1) = 1;
    mor         = repmat(mor,[length(sel1) 1]);
    Nij(sel1,:,c)  = Nij(sel1,:,c) + mor;
    clear mor
  end
end
Nij2 = Nij;
toc
% normalized versions
load SO_clustered_detections.mat;             cstdet1 = detections;
load SO_clustered_consistent_detections.mat;  cstdet2 = detections;
% normvec1 = cellfun(@length,{cstdet1.cluster_number});
% normvec2 = cellfun(@length,{cstdet2.cluster_number});
% Nij1norm = Nij1 ./ repmat(normvec1,[nchan 1]);
% Nij2norm = Nij2 ./ repmat(normvec2,[nchan 1]);
Nij1norm = Nij1;
Nij1norm(:,:,1) = Nij1norm(:,:,1) ./ repmat(cellfun(@length,{cstdet1.pospeak_cluster_number}),[nchan 1]);
Nij1norm(:,:,2) = Nij1norm(:,:,2) ./ repmat(cellfun(@length,{cstdet1.negpeak_cluster_number}),[nchan 1]);
Nij2norm = Nij2;
Nij2norm(:,:,1) = Nij2norm(:,:,1) ./ repmat(cellfun(@length,{cstdet2.pospeak}),[nchan 1]);
Nij2norm(:,:,2) = Nij2norm(:,:,2) ./ repmat(cellfun(@length,{cstdet2.negpeak}),[nchan 1]);

% clim1 = [0 .5*round(1.1*max(Nij1(:)))];
% clim2 = [0 .5*round(1.1*max(Nij2(:)))];
% clim3 = [0 .5];
% figure
% subplot(2,2,1),imagesc(Nij1(:,:,1));     caxis(clim1); title('cond1 (N co-involved)');                  fig(5).subplot(1).image = Nij1(:,:,1);     fig(5).subplot(1).title = 'cond1 (N co-involved)';
% subplot(2,2,2),imagesc(Nij1(:,:,2));     caxis(clim1); title('cond2 (N co-involved)');                  fig(5).subplot(2).image = Nij1(:,:,2);     fig(5).subplot(2).title = 'cond2 (N co-involved)';
% subplot(2,2,3),imagesc(Nij1norm(:,:,1)); caxis(clim3); title('cond1 (N co-involved / N-in-chan, row)'); fig(5).subplot(3).image = Nij1norm(:,:,1); fig(5).subplot(3).title = 'cond1 (N co-involved / N-in-chan, row)';
% subplot(2,2,4),imagesc(Nij1norm(:,:,2)); caxis(clim3); title('cond2 (N co-involved / N-in-chan, row)'); fig(5).subplot(4).image = Nij1norm(:,:,2); fig(5).subplot(4).title = 'cond2 (N co-involved / N-in-chan, row)';
% 
% figure
% subplot(2,2,1),imagesc(Nij2(:,:,1));     caxis(clim2); title('consistency cond1 (N co-involved)');                  fig(5).subplot(5).image = Nij2(:,:,1);     fig(5).subplot(5).title = 'consistency cond1 (N co-involved)';
% subplot(2,2,2),imagesc(Nij2(:,:,2));     caxis(clim2); title('consistency cond2 (N co-involved)');                  fig(5).subplot(6).image = Nij2(:,:,2);     fig(5).subplot(6).title = 'consistency cond2 (N co-involved)';
% subplot(2,2,3),imagesc(Nij2norm(:,:,1)); caxis(clim3); title('consistency cond1 (N co-involved / N-in-chan, row)'); fig(5).subplot(7).image = Nij2norm(:,:,1); fig(5).subplot(7).title = 'consistency cond1 (N co-involved / N-in-chan, row)';
% subplot(2,2,4),imagesc(Nij2norm(:,:,2)); caxis(clim3); title('consistency cond2 (N co-involved / N-in-chan, row)'); fig(5).subplot(8).image = Nij2norm(:,:,2); fig(5).subplot(8).title = 'consistency cond2 (N co-involved / N-in-chan, row)';

% grad1 & grad2 together
clim1 = [0 .4*round(1.1*max(Nij1(:)))];
clim2 = [0 .4*round(1.1*max(Nij2(:)))];
clim3 = [0 .3];
figure('name',SubjID); fig(5).name = SubjID;
subplot(2,2,1),imagesc(Nij1(:,:,1));     caxis(clim1); axis square; title('cond1 (N co-involved)');                  fig(5).subplot(1).image = Nij1(:,:,1);     fig(5).subplot(1).title = 'cond1 (N co-involved)';
subplot(2,2,2),imagesc(Nij1(:,:,2));     caxis(clim1); axis square; title('cond2 (N co-involved)');                  fig(5).subplot(2).image = Nij1(:,:,2);     fig(5).subplot(2).title = 'cond2 (N co-involved)';
subplot(2,2,3),imagesc(Nij1norm(:,:,1)); caxis(clim3); axis square; title('cond1 (N co-involved / N-in-chan, row)'); fig(5).subplot(3).image = Nij1norm(:,:,1); fig(5).subplot(3).title = 'cond1 (N co-involved / N-in-chan, row)';
subplot(2,2,4),imagesc(Nij1norm(:,:,2)); caxis(clim3); axis square; title('cond2 (N co-involved / N-in-chan, row)'); fig(5).subplot(4).image = Nij1norm(:,:,2); fig(5).subplot(4).title = 'cond2 (N co-involved / N-in-chan, row)';
figure('name',SubjID);
subplot(2,2,1),imagesc(Nij2(:,:,1));     caxis(clim2); axis square; title('consistency cond1 (N co-involved)');                  fig(5).subplot(5).image = Nij2(:,:,1);     fig(5).subplot(5).title = 'consistency cond1 (N co-involved)';
subplot(2,2,2),imagesc(Nij2(:,:,2));     caxis(clim2); axis square; title('consistency cond2 (N co-involved)');                  fig(5).subplot(6).image = Nij2(:,:,2);     fig(5).subplot(6).title = 'consistency cond2 (N co-involved)';
subplot(2,2,3),imagesc(Nij2norm(:,:,1)); caxis(clim3); axis square; title('consistency cond1 (N co-involved / N-in-chan, row)'); fig(5).subplot(7).image = Nij2norm(:,:,1); fig(5).subplot(7).title = 'consistency cond1 (N co-involved / N-in-chan, row)';
subplot(2,2,4),imagesc(Nij2norm(:,:,2)); caxis(clim3); axis square; title('consistency cond2 (N co-involved / N-in-chan, row)'); fig(5).subplot(8).image = Nij2norm(:,:,2); fig(5).subplot(8).title = 'consistency cond2 (N co-involved / N-in-chan, row)';
% grad1 only
clim1 = [0 .4*round(1.1*max(Nij1(:)))];
clim2 = [0 .4*round(1.1*max(Nij2(:)))];
clim3 = [0 .3];
sel   = 1:2:nchan;
figure('name',SubjID);
subplot(2,2,1),imagesc(Nij1(sel,sel,1));     caxis(clim1); axis square; title('grad1: cond1 (N co-involved)');
subplot(2,2,2),imagesc(Nij1(sel,sel,2));     caxis(clim1); axis square; title('cond2 (N co-involved)');
subplot(2,2,3),imagesc(Nij1norm(sel,sel,1)); caxis(clim3); axis square; title('cond1 (N co-involved / N-in-chan, row)');
subplot(2,2,4),imagesc(Nij1norm(sel,sel,2)); caxis(clim3); axis square; title('cond2 (N co-involved / N-in-chan, row)');
figure('name',SubjID);
subplot(2,2,1),imagesc(Nij2(sel,sel,1));     caxis(clim2); axis square; title('grad1: consistency cond1 (N co-involved)');
subplot(2,2,2),imagesc(Nij2(sel,sel,2));     caxis(clim2); axis square; title('consistency cond2 (N co-involved)');
subplot(2,2,3),imagesc(Nij2norm(sel,sel,1)); caxis(clim3); axis square; title('consistency cond1 (N co-involved / N-in-chan, row)');
subplot(2,2,4),imagesc(Nij2norm(sel,sel,2)); caxis(clim3); axis square; title('consistency cond2 (N co-involved / N-in-chan, row)');
% grad2 only
clim1 = [0 .4*round(1.1*max(Nij1(:)))];
clim2 = [0 .4*round(1.1*max(Nij2(:)))];
clim3 = [0 .3];
sel   = 2:2:nchan;
figure('name',SubjID);
subplot(2,2,1),imagesc(Nij1(sel,sel,1));     caxis(clim1); axis square; title('grad2: cond1 (N co-involved)');
subplot(2,2,2),imagesc(Nij1(sel,sel,2));     caxis(clim1); axis square; title('cond2 (N co-involved)');
subplot(2,2,3),imagesc(Nij1norm(sel,sel,1)); caxis(clim3); axis square; title('cond1 (N co-involved / N-in-chan, row)');
subplot(2,2,4),imagesc(Nij1norm(sel,sel,2)); caxis(clim3); axis square; title('cond2 (N co-involved / N-in-chan, row)');
figure('name',SubjID);
subplot(2,2,1),imagesc(Nij2(sel,sel,1));     caxis(clim2); axis square; title('grad2: consistency cond1 (N co-involved)');
subplot(2,2,2),imagesc(Nij2(sel,sel,2));     caxis(clim2); axis square; title('consistency cond2 (N co-involved)');
subplot(2,2,3),imagesc(Nij2norm(sel,sel,1)); caxis(clim3); axis square; title('consistency cond1 (N co-involved / N-in-chan, row)');
subplot(2,2,4),imagesc(Nij2norm(sel,sel,2)); caxis(clim3); axis square; title('consistency cond2 (N co-involved / N-in-chan, row)');

%% 6. Functions of anterior-posterior distance.
% Need two new distance metrics:
%   i. anterior-posterior distance between pairs of channels
%     => scatter plot delays for all channel pairs vs distance
%   ii. anterior-posterior distance between 1 channel & anterior pole
%     => scatter plot delays for all channels vs distance


%% 7. Streamline analysis of 1 single peak wave & 1 multipeak wave



%% 8. Source localization: 1 single peak wave & 1 multipeak wave (same as 7)



%% - compare percentage of multipeak waves between early & late sleep epochs.



%% 
% % make a list of cluster indices for each channel
% cluster_list = [];
% for ch = 1:nchan
%   this = flipdat.sensor_info(ch).label;
%   c = 1; cluster_list(c).channel(ch).label           = flipdat.sensor_info(ch).label;
%   c = 1; cluster_list(c).channel(ch).cluster_indices = find(arrayfun(@(x)(ismember(this,x.InvolvedChans)),clusters(c).epochs));
%   c = 1; cluster_list(c).channel(ch).channel_indices = arrayfun(@(x)find(ismember(x.InvolvedChans,this)),clusters(c).epochs(cluster_list(c).channel(ch).cluster_indices));
%   c = 2; cluster_list(c).channel(ch).label           = flipdat.sensor_info(ch).label;
%   c = 2; cluster_list(c).channel(ch).cluster_indices = find(arrayfun(@(x)(ismember(this,x.InvolvedChans)),clusters(c).epochs));
%   c = 2; cluster_list(c).channel(ch).channel_indices = arrayfun(@(x)find(ismember(x.InvolvedChans,this)),clusters(c).epochs(cluster_list(c).channel(ch).cluster_indices));  
% end
% 
% % channelcmb = [];
% % band = {[1 3],[4 8],[8 12],[12 25],[25 55],[20 60],[10 80]};
% 
% % create list of all sensor combinations
% cmb = [];
% nchan = length(sens);
% c1  = [repmat(1:nchan,[nchan 1])']'; c1 = c1(:);
% c2  = repmat([1:nchan]',[nchan 1]);
% cmb = [c1 c2];
% cup = cmb(c2> c1,:);  % upper triangular off-diagonal matrix
% cdn = cmb(c2< c1,:);  % lower triangular off-diagonal matrix
% cdg = cmb(c2==c1,:);  % diagonal
% [jnk cupidx] = setdiff(cup,fliplr(cdn),'rows');
% [jnk cdnidx] = setdiff(fliplr(cdn),cup,'rows');
% [jnk ix tss] = intersect(cup,fliplr(cdn),'rows');
% cupidx = [cupidx ix];
% cmb = sortrows([cup(cupidx,:); cdg; cdn(cdnidx,:)]);
% 
