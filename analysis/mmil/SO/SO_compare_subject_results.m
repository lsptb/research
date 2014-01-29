function [res,grand] = SO_compare_subject_results(subjects,ClusterFile,DetectionFile,StructSelection,ArrayIndex,outfile,SetID)
% fixed subject-specific mean delay calculation
if nargin < 7, SetID   = date; end
if nargin < 6, outfile = [];   end
if nargin < 5, error('Some required parameters were not specified'); end

% subjects        = [1 2 4 5 6 8];
% % Example: multipeak, noflip, all trials, earliest detection
% outfile         = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/images/noflip/alltrials/workspace_AllSubjects_Summary_EarliestDetection_noflip_alltrials.mat';
% ClusterFile     = 'SO_clusters_noflip.mat';       % noflip, flipmatrix, consistency
% DetectionFile   = 'SO_clustered_detections.mat';  % clustered_detections, clustered_consistent_detections, clustered_flipmatrix_detections
% StructSelection = 'clusters_Near_all';            % Near_all, Nmax_all, Cear_all, Cmax_all; X_pks(1or2)
% ArrayIndex      = 1;                              % 1 for pos or all, 2 for neg
AllDetFile      = 'SO_detections.mat';            % SO_detections or SO_detections_singlepeaks
SpeedBins       = linspace(0,10,15);                           % m/s, bin centers
RBins           = linspace(-1,1,100);

tic
nsubj       = length(subjects);
res         = [];
grand       = [];
grand.Rplot.Rall_count          = zeros(1,length(RBins));%[];
grand.Rplot.Rsig_count          = zeros(1,length(RBins));%[];
grand.SpeedPlot.Speed_count     = zeros(1,length(SpeedBins));%[];
grand.DelayStdRplot.alltaustd   = [];
grand.DelayStdRplot.allR        = [];
grand.DelayStdRplot.sigtaustd   = [];
grand.DelayStdRplot.sigR        = [];
grand.DelayStdRplot.insigtaustd = [];
grand.DelayStdRplot.insigR      = [];
screensize = get(0,'ScreenSize');
figure('name',sprintf('Group results (%s)',SetID),'color','w','Position',[1 1 1*screensize(3) 1*screensize(4)]);

for s = 1:nsubj
  subj    = subjects(s);
  params  = SO_params(subj);
  SubjID  = sprintf('s%g',subj);
  cd(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/%s/matfiles',SubjID));
  load(AllDetFile,'detections');    alldet = detections;
  load(DetectionFile,'detections'); det    = detections;
  tmp       = load(ClusterFile,StructSelection);
  tmp       = tmp.(StructSelection);
  clusters  = tmp(ArrayIndex);
  clear tmp detections
  % NOTE: alldet, det, & clusters will be used for everything
  
  % define some variables & create list of involved trials for each channel
  Rcheck          = cellfun(@isempty,{clusters.epochs.R});
  clustsiz        = cellfun(@length,{clusters.epochs.InvolvedChans});
  clusters.epochs = clusters.epochs(clustsiz>1 & Rcheck~=1);
  clustsiz        = clustsiz(clustsiz>1 & Rcheck~=1);
  nclust  = length(clusters.epochs);
  sens    = clusters.sensor_info;
  allR    = [clusters.epochs.R];
  allp    = [clusters.epochs.p];
  allN    = [clusters.epochs.N];
  allrefs = {clusters.epochs.RefChan};
  alltaus = {clusters.epochs.Delays};
  alldist = {clusters.epochs.Theta3D};
  alltime = {clusters.epochs.DetectionTimes};
  g1      = strmatch('grad1',{sens.typestring});
  g2      = strmatch('grad2',{sens.typestring});
  labels  = {sens.label};
  nchan   = length(labels);
  tmp     = {clusters.epochs.InvolvedChans}; 
  tmp     = cellfun(@(x)[x{:}],tmp,'uniformoutput',false);
  [trlind{1:nchan}] = deal([]);
  [invind{1:nchan}] = deal([]);
  [delays{1:nchan}] = deal([]);
  for k   = 1:nchan
    trlind{k} = find(~cellfun(@isempty,regexp(tmp,labels{k})));
    if isempty(trlind{k}), continue; end
    invind{k} = cellfun(@(x)strmatch(labels{k},x),{clusters.epochs(trlind{k}).InvolvedChans});
    delays{k} = cellfun(@(x,y)x(y),alltaus(trlind{k}),num2cell(invind{k}));
  end
  clear tmp k
  
  % 1. Detection density (combine grad1 & grad2)
  Ndet              = cellfun(@length,trlind);
  res(s).detection  = (Ndet(g1) + Ndet(g2)) / (2*nclust);
  
  % 2. Origin density
  orgind            = cellfun(@(x,y)strmatch(x,allrefs),labels,'uniformoutput',false);
  Norg              = cellfun(@length,orgind);
  res(s).origin     = (Norg(g1) + Norg(g2)) / nclust;
  
  % 3. Average delay
  avgdelay          = cellfun(@mean,delays);
  res(s).delay      = (Ndet(g1).*avgdelay(g1) + Ndet(g2).*avgdelay(g2)) ./ (Ndet(g1) + Ndet(g2));
  
  % 4. Speed histogram  
  r       = (3*params.BrainVolume/(4*pi))^(1/3);
  c       = (2*pi*r/360)*(1/100);
  v       = cellfun(@(x,y)c*mean(x(x~=0 & y~=0) ./ y(x~=0 & y~=0)),alldist,alltaus);
  [Ni,Vi] = hist(v,SpeedBins);
  res(s).SpeedPlot.Speed_bin    = Vi;
  res(s).SpeedPlot.Speed_count  = Ni;
  grand.SpeedPlot.Speed_bin     = Vi;
  grand.SpeedPlot.Speed_count   = grand.SpeedPlot.Speed_count + Ni;
  
  % 5. R(delay,distance) histogram
  [Ni,Ri] = hist(allR,RBins);
  res(s).Rplot.Rall_bin   = Ri; grand.Rplot.R_bin      = Ri;
  res(s).Rplot.Rall_count = Ni; grand.Rplot.Rall_count = grand.Rplot.Rall_count + Ni;
  [Ni,Ri] = hist(allR(allp<.05),RBins);
  res(s).Rplot.Rsig_bin   = Ri; 
  res(s).Rplot.Rsig_count = Ni; grand.Rplot.Rsig_count = grand.Rplot.Rsig_count + Ni;
  res(s).Rplot.SigFraction= sum(res(s).Rplot.Rsig_count) / sum(res(s).Rplot.Rall_count);

  % 6. R(delay,distance) vs std(delay)
%   taustd = cellfun(@std,alltaus);
%   res(s).DelayStdRplot.alltaustd   = taustd;          grand.DelayStdRplot.alltaustd   = [grand.DelayStdRplot.alltaustd taustd];
%   res(s).DelayStdRplot.allR        = abs(allR);       grand.DelayStdRplot.allR        = [grand.DelayStdRplot.allR abs(allR)];
%   sel                              = allp < .05 & allR > min(.8*max(allR),.5);
%   res(s).DelayStdRplot.sigtaustd   = taustd(sel);     grand.DelayStdRplot.sigtaustd   = [grand.DelayStdRplot.sigtaustd taustd(sel)];
%   res(s).DelayStdRplot.sigR        = abs(allR(sel));  grand.DelayStdRplot.sigR        = [grand.DelayStdRplot.sigR abs(allR(sel))];
%   sel                              = allp > .05 | allR < min(.8*max(allR),.5);
%   res(s).DelayStdRplot.insigtaustd = taustd(sel);     grand.DelayStdRplot.insigtaustd = [grand.DelayStdRplot.insigtaustd taustd(sel)];
%   res(s).DelayStdRplot.insigR      = abs(allR(sel));  grand.DelayStdRplot.insigR      = [grand.DelayStdRplot.insigR abs(allR(sel))];
%   [tmpR,tmpp] = corrcoef([res(s).DelayStdRplot.alltaustd' res(s).DelayStdRplot.allR']);
%   if numel(tmpR)>1
%     res(s).DelayStdRplot.allCorrCoef_R = tmpR(1,2);
%     res(s).DelayStdRplot.allCorrCoef_p = tmpp(1,2);
%   end
%   [tmpR,tmpp] = corrcoef([res(s).DelayStdRplot.sigtaustd' res(s).DelayStdRplot.sigR']);
%   if numel(tmpR)>1
%     res(s).DelayStdRplot.sigCorrCoef_R = tmpR(1,2);
%     res(s).DelayStdRplot.sigCorrCoef_p = tmpp(1,2);
%   end
%   [tmpR,tmpp] = corrcoef([res(s).DelayStdRplot.insigtaustd' res(s).DelayStdRplot.insigR']);
%   if numel(tmpR)>1
%     res(s).DelayStdRplot.insigCorrCoef_R = tmpR(1,2);
%     res(s).DelayStdRplot.insigCorrCoef_p = tmpp(1,2);
%   end
%   clear tmpR tmpp
  
  % 7. Delay vs A-P separation distance
% Angular distance
  clear APdelay APangle APstep
  theta3D = ts_BetweenSensorAngles(sens,.08);
  theta3D(isnan(theta3D)) = 0;  
  midline_g1   = [57 59 69 43 45 75 49 51 55 53 137 171 149 151 155 153 157 159];
  midline_pos1 = [0 1 2 3 4 4 5 5 6 6 7 7 8 8 9 9 10 11];
  midline_g2   = midline_g1 + 1; % highlight = [midline_g1 midline_g2]; figure; dewar;
  midline_pos2 = [0 1 2 3 4 4 5 5 6 6 7 7 8 8 9 9 10 11];
  midline_chan = [midline_g1 midline_g2];
  [midline_chan,ix] = sort(midline_chan);
  midline_pos  = [midline_pos1 midline_pos2];
  midline_pos  = midline_pos(ix);
  % [APdelay{1:length(midline_chan)^2/2}] = deal([]);
  % APdistance = nan(length(midline_chan)^2/2,1);
  cnt = 1;
  for i = 1:length(midline_chan)
    I   = midline_chan(i);
    for j = 1:length(midline_chan)
      if (j > i) || (midline_pos(i)==midline_pos(j)), break; end        % lower-triangular
      J = midline_chan(j);
      [ind,II,JJ] = intersect(trlind{I},trlind{J});
      if isempty(ind), continue; end
      alltimes1   = cellfun(@(x,y)x(y),alltime(ind),num2cell(invind{I}(II)));
      alltimes2   = cellfun(@(x,y)x(y),alltime(ind),num2cell(invind{J}(JJ)));
      if midline_pos(i) < midline_pos(j)
        APdelay{cnt} = alltimes2 - alltimes1;
      else
        APdelay{cnt} = alltimes1 - alltimes2;
      end
      APangle(cnt)   = theta3D(I,J);
      APstep(cnt)    = abs(midline_pos(i) - midline_pos(j));
      cnt = cnt + 1;
    end
  end
  res(s).APplot.APdelay = APdelay; %if s==1, grand.APplot.APdelay = APdelay; else grand.APplot.APdelay = {grand.APplot.APdelay{:} APdelay{:}}; end %grand.APplot.APdelay = cellfun(@(x,y)[x y],grand.APplot.APdelay,APdelay,'uniformoutput',false); end% sec
  res(s).APplot.APangle = APangle; %if s==1, grand.APplot.APangle = APangle; else grand.APplot.APangle = {grand.APplot.APangle{:} APangle{:}}; end %cellfun(@(x,y)[x y],grand.APplot.APangle,APangle,'uniformoutput',false); end% deg
  res(s).APplot.APstep  = APstep;  %if s==1, grand.APplot.APstep = APstep;   else grand.APplot.APstep  = {grand.APplot.APstep{:} APstep{:}};   end %cellfun(@(x,y)[x y],grand.APplot.APstep,APstep,'uniformoutput',false);   end% unit
  res(s).APplot.APavgdelay = cellfun(@mean,APdelay);
  vals  = unique(res(s).APplot.APstep); clear tmp
  for i = 1:length(vals)
    tmp(i) = mean(res(s).APplot.APavgdelay(res(s).APplot.APstep==vals(i)));
  end
  res(s).APplot.APmean_steps    = vals; if s==1, grand.APplot.APmean_steps    = vals; end
  res(s).APplot.APmean_avgdelay = tmp;  if s==1, grand.APplot.APmean_avgdelay = tmp/nsubj; else grand.APplot.APmean_avgdelay = grand.APplot.APmean_avgdelay + tmp/nsubj; end
  clear tmp vals
    
  [tmpR,tmpp] = corrcoef([res(s).APplot.APstep' res(s).APplot.APavgdelay']);
  if numel(tmpR)>1
    res(s).APplot.DelayStep_R = tmpR(1,2);
    res(s).APplot.DelayStep_p = tmpp(1,2);
  end  
  [tmpR,tmpp] = corrcoef([res(s).APplot.APmean_steps' res(s).APplot.APmean_avgdelay']);
  if numel(tmpR)>1
    res(s).APplot.MeanDelayStep_R = tmpR(1,2);
    res(s).APplot.MeanDelayStep_p = tmpp(1,2);
  end  
  [tmpR,tmpp] = corrcoef([res(s).APplot.APangle' res(s).APplot.APavgdelay']);
  if numel(tmpR)>1
    res(s).APplot.DelayAngle_R = tmpR(1,2);
    res(s).APplot.DelayAngle_p = tmpp(1,2);
  end     

  % Plot subject 
  topostyle = 'straight';
  % 1. detections
  X   = res(s).detection'; zlim = [];
  tmp = ts_matrix2avg([X X],'sens',sens(g1));
  tmp.averages.data = X;
  tmp.averages.time = 0;
  subplot(6,7,s); ts_ezplot(tmp,'highlight',[],'electrodes','on','zlim',zlim,'style',topostyle,'topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
  if s == 1, ylabel('detection'); end; title(params.SubjID); axis equal;
    h = get(gca,'Children');
    delete(h(1));
    set(h(2:5),'LineWidth',.2);
    set(h(6),'Marker','.','MarkerSize',.1);  
    set(gca,'clim',[0 .5]);
  % 2. origins
  X   = res(s).origin'; zlim = [];
  tmp = ts_matrix2avg([X X],'sens',sens(g1));
  tmp.averages.data = X;
  tmp.averages.time = 0;
  subplot(6,7,7+s); ts_ezplot(tmp,'highlight',[],'electrodes','on','zlim',zlim,'style',topostyle,'topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
  if s == 1, ylabel('origin'); end; axis equal;
    h = get(gca,'Children');
    delete(h(1));
    set(h(2:5),'LineWidth',.2);
    set(h(6),'Marker','.','MarkerSize',.1);  
    set(gca,'clim',[0 .03]);
  % 3. average delay
  X   = res(s).delay'; zlim = [];
  tmp = ts_matrix2avg([X X],'sens',sens(g1));
  tmp.averages.data = X;
  tmp.averages.time = 0;
  tmpbadchans = find(isnan(X));
  subplot(6,7,14+s); ts_ezplot(tmp,'badchans',tmpbadchans,'highlight',[],'electrodes','on','zlim',zlim,'style',topostyle,'topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
  if s == 1, ylabel('delay'); end; axis equal;
    h = get(gca,'Children');
    delete(h(1));
    set(h(2:5),'LineWidth',.2);
    set(h(6),'Marker','.','MarkerSize',.1);  
    set(gca,'clim',[.08 .1]);
  % 4. speed
%   subplot(6,7,21+s); bar(res(s).SpeedPlot.Speed_bin,res(s).SpeedPlot.Speed_count);
  nval = sum(res(s).SpeedPlot.Speed_count);
  xval = res(s).SpeedPlot.Speed_bin;
  yval = 100*res(s).SpeedPlot.Speed_count/nval;
  subplot(6,7,21+s); bar(xval,yval);
  set(gca,'xlim',[0 10]); if s==1, ylabel('speed'); end
  set(gca,'ylim',[0 60]);
  text(5,.8*max(get(gca,'ylim')),sprintf('%gm/s',sum(res(s).SpeedPlot.Speed_bin.*res(s).SpeedPlot.Speed_count)./sum(res(s).SpeedPlot.Speed_count)));
    set(gca,'FontSize',5);
    set(get(gca,'xlabel'),'fontsize',5);
    set(get(gca,'ylabel'),'fontsize',5);
    h = get(gca,'Children');
    ylim=get(gca,'ylim'); set(h(1),'Position',[7 .8*ylim(2)],'FontSize',6);      
  % 5. R histogram
%   subplot(6,7,28+s),plot(res(s).Rplot.Rall_bin,res(s).Rplot.Rall_count,'b.-',res(s).Rplot.Rsig_bin,res(s).Rplot.Rsig_count,'r.-');
  nval1 = sum(res(s).Rplot.Rall_count);
  xval1 = res(s).Rplot.Rall_bin;
  yval1 = 100*res(s).Rplot.Rall_count/nval1;
%   nval2 = sum(res(s).Rplot.Rsig_count);
  xval2 = res(s).Rplot.Rsig_bin;
  yval2 = 100*res(s).Rplot.Rsig_count/nval1;
  subplot(6,7,28+s),plot(xval1,yval1,'b.-',xval1,yval2,'r.-');
  set(gca,'ylim',[0 5]);
  set(gca,'xlim',[-1 1]); vline(0,'k'); if s==1, ylabel('R(tau,theta)'); end
  text(-1,.9*max(get(gca,'ylim')),sprintf('%g%% sig',res(s).Rplot.SigFraction));
    set(gca,'FontSize',5);
    set(get(gca,'xlabel'),'fontsize',5);
    set(get(gca,'ylabel'),'fontsize',5);
    h = get(gca,'Children');
    ylim=get(gca,'ylim'); set(h(1),'Position',[-.9 .8*ylim(2)],'FontSize',6);
    delete(h(2));
    set(h(3),'MarkerSize',5);
    vline(0,'k');  
%   % 6. R vs std(tau)
%   subplot(6,7,35+s)
%   plot(res(s).DelayStdRplot.sigtaustd,res(s).DelayStdRplot.sigR,'r.',res(s).DelayStdRplot.insigtaustd,res(s).DelayStdRplot.insigR,'b.');
%   if s==1,xlabel('std(delay)'),ylabel('R(delay,dist)'); end
%   h=legend(sprintf('(p<.05): R=%g,p=%g',res(s).DelayStdRplot.sigCorrCoef_R,res(s).DelayStdRplot.sigCorrCoef_p),...
%       sprintf('(p>.05): R=%g,p=%g',res(s).DelayStdRplot.insigCorrCoef_R,res(s).DelayStdRplot.insigCorrCoef_p));
%   set(h,'Location','NorthOutside','FontSize',6);  
  % 7. AP distance vs tau    
  subplot(6,7,35+s),plot(res(s).APplot.APstep,res(s).APplot.APavgdelay,'b.',res(s).APplot.APmean_steps,res(s).APplot.APmean_avgdelay,'ro');
  set(gca,'ylim',[-.03 .03]); set(gca,'xlim',[0 max(res(s).APplot.APstep)+1]); hline(0,'k');
  if s==1, xlabel('A-P separation'); ylabel('delay'); end
%   subplot(6,7,42+s),plot(res(s).APplot.APangle,res(s).APplot.APavgdelay,'b.');
%   set(gca,'ylim',[-.1 .1]); hline(0,'k');
      set(gca,'FontSize',5);
      set(get(gca,'xlabel'),'fontsize',5);
      set(get(gca,'ylabel'),'fontsize',5);
      set(gca,'xlim',[0 12]);
      h = get(gca,'Children');
      set(h(2),'MarkerSize',2);
      set(h(3),'MarkerSize',2);
    % Build grand average
    NinvClust          = cellfun(@length,trlind);
    grand.nclusters(s) = nclust;
%     if s == 1, grand.Ndet = nclust*res(s).detection; else grand.Ndet = grand.Ndet + nclust*res(s).detection; end
    if s == 1, grand.Ndet = (Ndet(g1) + Ndet(g2)); else grand.Ndet = grand.Ndet + (Ndet(g1) + Ndet(g2)); end
    if s == 1, grand.Norg = nclust*res(s).origin;    else grand.Norg = grand.Norg + nclust*res(s).origin;    end
%     if s == 1, grand.Ndel = NinvClust.*res(s).delay; else grand.Ndel = grand.Ndel + NinvClust.*res(s).delay; end
    if s == 1, grand.Ndel = (Ndet(g1) + Ndet(g2)).*res(s).delay; else grand.Ndel = grand.Ndel + (Ndet(g1) + Ndet(g2)).*res(s).delay; end
  toc
end
grand.detection = grand.Ndet / sum(grand.nclusters);
grand.origin    = grand.Norg / sum(grand.nclusters);
grand.delay     = grand.Ndel ./ grand.Ndet;
% tmp             = arrayfun(@(x)x.SpeedPlot.Speed_count,res,'uniformoutput',false);
% grand.SpeedPlot.Speed_count = sum(cat(1,tmp{:}),1);
grand.Rplot.SigFraction= sum(grand.Rplot.Rsig_count) / sum(grand.Rplot.Rall_count);
% [tmpR,tmpp] = corrcoef([grand.DelayStdRplot.alltaustd' grand.DelayStdRplot.allR']);
% if numel(tmpR)>1
%   grand.DelayStdRplot.allCorrCoef_R = tmpR(1,2);
%   grand.DelayStdRplot.allCorrCoef_p = tmpp(1,2);
% end
% [tmpR,tmpp] = corrcoef([grand.DelayStdRplot.sigtaustd' grand.DelayStdRplot.sigR']);
% if numel(tmpR)>1
%   grand.DelayStdRplot.sigCorrCoef_R = tmpR(1,2);
%   grand.DelayStdRplot.sigCorrCoef_p = tmpp(1,2);
% end
% [tmpR,tmpp] = corrcoef([grand.DelayStdRplot.insigtaustd' grand.DelayStdRplot.insigR']);
% if numel(tmpR)>1
%   grand.DelayStdRplot.insigCorrCoef_R = tmpR(1,2);
%   grand.DelayStdRplot.insigCorrCoef_p = tmpp(1,2);
% end
clear tmpR tmpp

%% Plot grand average
% 1. detections
X   = grand.detection'; zlim = [];
tmp = ts_matrix2avg([X X],'sens',sens(g1));
tmp.averages.data = X; tmp.averages.time = 0;
subplot(6,7,7); ts_ezplot(tmp,'highlight',[],'electrodes','on','zlim',zlim,'style',topostyle,'topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
axis equal; title('grand average');
    h = get(gca,'Children');
    delete(h(1));
    set(h(2:5),'LineWidth',.2);
    set(h(6),'Marker','.','MarkerSize',.1);
    set(gca,'clim',[0 1]);
% 2. origins
X   = grand.origin'; zlim = [];
tmp = ts_matrix2avg([X X],'sens',sens(g1));
tmp.averages.data = X; tmp.averages.time = 0;
subplot(6,7,14); ts_ezplot(tmp,'highlight',[],'electrodes','on','zlim',zlim,'style',topostyle,'topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
if s == 1, ylabel('origin'); end; axis equal;
    h = get(gca,'Children');
    delete(h(1));
    set(h(2:5),'LineWidth',.2);
    set(h(6),'Marker','.','MarkerSize',.1);
    set(gca,'clim',[0 .03]);
% 3. average delay
X   = grand.delay'; zlim = [];
tmp = ts_matrix2avg([X X],'sens',sens(g1));
tmp.averages.data = X; tmp.averages.time = 0;
tmpbadchans = find(isnan(X));
subplot(6,7,21); ts_ezplot(tmp,'badchans',tmpbadchans,'highlight',[],'electrodes','on','zlim',zlim,'style',topostyle,'topoplot',1,'toprows',1,'topcols',1,'layout','/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay','newfig',0);
if s == 1, ylabel('delay'); end; axis equal;
    h = get(gca,'Children');
    delete(h(1));
    set(h(2:5),'LineWidth',.2);
    set(h(6),'Marker','.','MarkerSize',.1);
    set(gca,'clim',[.08 .1]);
% 4. speed
  nval = sum(grand.SpeedPlot.Speed_count);
  xval = grand.SpeedPlot.Speed_bin;
  yval = 100*grand.SpeedPlot.Speed_count/nval;
  subplot(6,7,28); bar(xval,yval);
  set(gca,'ylim',[0 60]);
% subplot(6,7,28); bar(grand.SpeedPlot.Speed_bin,grand.SpeedPlot.Speed_count);
set(gca,'xlim',[0 10]);
text(5,.8*max(get(gca,'ylim')),sprintf('%gm/s',sum(grand.SpeedPlot.Speed_bin.*grand.SpeedPlot.Speed_count)./sum(grand.SpeedPlot.Speed_count)));
    set(gca,'FontSize',5);
    set(get(gca,'xlabel'),'fontsize',5);
    set(get(gca,'ylabel'),'fontsize',5);
    h = get(gca,'Children');
    ylim=get(gca,'ylim'); set(h(1),'Position',[7 .8*ylim(2)],'FontSize',6);      
% 5. R histogram
  nval1 = sum(grand.Rplot.Rall_count);
  xval1 = grand.Rplot.R_bin;
  yval1 = 100*grand.Rplot.Rall_count/nval1;
  xval2 = grand.Rplot.R_bin;
  yval2 = 100*grand.Rplot.Rsig_count/nval1;
  subplot(6,7,35),plot(xval1,yval1,'b.-',xval1,yval2,'r.-');
  set(gca,'ylim',[0 5]);
% subplot(6,7,35),plot(grand.Rplot.R_bin,grand.Rplot.Rall_count,'b.-',grand.Rplot.R_bin,grand.Rplot.Rsig_count,'r.-');
set(gca,'xlim',[-1 1]); vline(0,'k');
text(-1,.9*max(get(gca,'ylim')),sprintf('%g%% sig',grand.Rplot.SigFraction));
    set(gca,'FontSize',5);
    set(get(gca,'xlabel'),'fontsize',5);
    set(get(gca,'ylabel'),'fontsize',5);
    h = get(gca,'Children');
    ylim=get(gca,'ylim'); set(h(1),'Position',[-.9 .8*ylim(2)],'FontSize',6);
    delete(h(2));
    set(h(3),'MarkerSize',5);
    vline(0,'k');  
% 6. R vs std(tau)
% subplot(6,7,42)
% plot(grand.DelayStdRplot.sigtaustd,grand.DelayStdRplot.sigR,'r.',grand.DelayStdRplot.insigtaustd,grand.DelayStdRplot.insigR,'b.');
% if s==1,xlabel('std(delay)'),ylabel('R(delay,dist)'); end
% legend(sprintf('(p<.05): R=%g,p=%g',grand.DelayStdRplot.sigCorrCoef_R,grand.DelayStdRplot.sigCorrCoef_p),...
%     sprintf('(p>.05): R=%g,p=%g',grand.DelayStdRplot.insigCorrCoef_R,grand.DelayStdRplot.insigCorrCoef_p));

% 7. AP distance vs tau    
subplot(6,7,42),plot(grand.APplot.APmean_steps,grand.APplot.APmean_avgdelay,'ro');
set(gca,'ylim',[-.003 .003]); set(gca,'xlim',[0 max(grand.APplot.APmean_steps)+1]); hline(0,'k');
if s==1, xlabel('A-P separation'); ylabel('delay'); end
      set(gca,'FontSize',5);
      set(get(gca,'xlabel'),'fontsize',5);
      set(get(gca,'ylabel'),'fontsize',5);
      set(gca,'xlim',[0 12]);
      h = get(gca,'Children');
      set(h(2),'MarkerSize',2);
% Save figure and workspace (return results)
if ~isempty(outfile)
  save(outfile,'res','grand','subjects','ClusterFile','DetectionFile','StructSelection','ArrayIndex');
end

toc
  