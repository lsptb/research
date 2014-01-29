%% run plot_delay_comparisons after this

%%
datafile = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/s1/matfiles/proc_epoch_data_1.mat';
load(datafile);
load SO_clustered_detections.mat

StepSize=.01; IntegrationWindow=.025;
t   = epoch_data.epochs.time;
Fs  =  epoch_data.sfreq;
  L     = floor(IntegrationWindow*Fs/2);
  n     = floor(StepSize*Fs);
  nsmp  = length(t);
  c     = L+1:n:nsmp-L;
  tc    = t(c);
tmp = ts_matrix2data(count,'time',tc);
visualizer(tmp);

  
%% s1 - interval validation
subj         = [1];
CorrCoefFile = 'matfiles/SO_corrcoef.mat';
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run1';
plot_intervals(subj,rootdir,CorrCoefFile);
subplot(4,1,2)
% Draw sleep stage lines based on hypnogram data
w      = 30;
colors = 'yrkgc'; L = [.5 3 3 .5 .5];
stage2 = [1 38 60 92 99 108 158 161 164 186];
for k =1:length(stage2), vline(stage2(k)*w,'Color',colors(1),'LineWidth',L(1)); end
stage3 = [12 17 56 61 64 76 157 160 163 167];
for k =1:length(stage3), vline(stage3(k)*w,'Color',colors(2),'LineWidth',L(2)); end
stage4 = [16 18 63 75];
for k =1:length(stage4), vline(stage4(k)*w,'Color',colors(3),'LineWidth',L(3)); end
stageW = [93 182];
for k =1:length(stageW), vline(stageW(k)*w,'Color',colors(4),'LineWidth',L(4)); end
stageR = 105;
for k =1:length(stageR), vline(stageR(k)*w,'Color',colors(5),'LineWidth',L(5)); end
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run2';
plot_intervals(subj,rootdir,CorrCoefFile);
subplot(4,1,2)
% Draw sleep stage lines based on hypnogram data
w      = 30;
colors = 'yrkgc'; L = [.5 3 3 .5 .5];
stage2 = [1 38 60 92 99 108 158 161 164 186];
for k =1:length(stage2), vline(stage2(k)*w,'Color',colors(1),'LineWidth',L(1)); end
stage3 = [12 17 56 61 64 76 157 160 163 167];
for k =1:length(stage3), vline(stage3(k)*w,'Color',colors(2),'LineWidth',L(2)); end
stage4 = [16 18 63 75];
for k =1:length(stage4), vline(stage4(k)*w,'Color',colors(3),'LineWidth',L(3)); end
stageW = [93 182];
for k =1:length(stageW), vline(stageW(k)*w,'Color',colors(4),'LineWidth',L(4)); end
stageR = 105;
for k =1:length(stageR), vline(stageR(k)*w,'Color',colors(5),'LineWidth',L(5)); end
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/RUN3';
plot_intervals(subj,rootdir,CorrCoefFile);
subplot(4,1,2)
% Draw sleep stage lines based on hypnogram data
w      = 30;
colors = 'yrkgc'; L = [.5 3 3 .5 .5];
stage2 = [1 38 60 92 99 108 158 161 164 186];
for k =1:length(stage2), vline(stage2(k)*w,'Color',colors(1),'LineWidth',L(1)); end
stage3 = [12 17 56 61 64 76 157 160 163 167];
for k =1:length(stage3), vline(stage3(k)*w,'Color',colors(2),'LineWidth',L(2)); end
stage4 = [16 18 63 75];
for k =1:length(stage4), vline(stage4(k)*w,'Color',colors(3),'LineWidth',L(3)); end
stageW = [93 182];
for k =1:length(stageW), vline(stageW(k)*w,'Color',colors(4),'LineWidth',L(4)); end
stageR = 105;
for k =1:length(stageR), vline(stageR(k)*w,'Color',colors(5),'LineWidth',L(5)); end

rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/RUN4';
plot_intervals(subj,rootdir,CorrCoefFile);
subplot(4,1,2)
% Draw sleep stage lines based on hypnogram data
w      = 30;
colors = 'yrkgc'; L = [.5 3 3 .5 .5];
stage2 = [1 38 60 92 99 108 158 161 164 186];
for k =1:length(stage2), vline(stage2(k)*w,'Color',colors(1),'LineWidth',L(1)); end
stage3 = [12 17 56 61 64 76 157 160 163 167];
for k =1:length(stage3), vline(stage3(k)*w,'Color',colors(2),'LineWidth',L(2)); end
stage4 = [16 18 63 75];
for k =1:length(stage4), vline(stage4(k)*w,'Color',colors(3),'LineWidth',L(3)); end
stageW = [93 182];
for k =1:length(stageW), vline(stageW(k)*w,'Color',colors(4),'LineWidth',L(4)); end
stageR = 105;
for k =1:length(stageR), vline(stageR(k)*w,'Color',colors(5),'LineWidth',L(5)); end


%% intervals: s2,s5,s6,s8
subj         = [2 5 6 8];
CorrCoefFile = 'matfiles/SO_corrcoef.mat';
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run1';
plot_intervals(subj,rootdir,CorrCoefFile);
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run2';
plot_intervals(subj,rootdir,CorrCoefFile);
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/RUN3';
plot_intervals(subj,rootdir,CorrCoefFile);
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/RUN4';
plot_intervals(subj,rootdir,CorrCoefFile);

%% corrcoef results: s2, s5, s6, s8
subj         = [2 5 6 8];
CorrCoefFile = 'matfiles/SO_corrcoef.mat';
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run1';
allresults1  = plot_results(subj,rootdir,CorrCoefFile); close all
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/run2';
allresults2  = plot_results(subj,rootdir,CorrCoefFile); close all
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/RUN3';
allresults3  = plot_results(subj,rootdir,CorrCoefFile); close all
rootdir      = '/space/emc2/1/halgdev/projects/sleep/MEG/SO/RUN4';
allresults4  = plot_results(subj,rootdir,CorrCoefFile); close all

allresults   = allresults1;
SubjectColors = 'bgrcmyk';
MethodSymbols = {'.','o','x','*'};
PeakLines     = {'--',':'};

% LineType = [SubjectColors(s) MethodSymbols{} PeakLines{}];
clear grandavg fig1 fig2
for s = 1:length(allresults)+1
  if s <= length(allresults)
    results = allresults(s);
    width   = .5;
  else
    warning off
    grandavg.summary.wrtFirstDet.SigFraction          = grandavg.summary.wrtFirstDet.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtFirstDet.NumTrials,[1 1 size(grandavg.summary.wrtFirstDet.NumSigTrials,3)]);
    grandavg.summary.wrtFirstDet.SigFraction(isnan(grandavg.summary.wrtFirstDet.SigFraction)) = 0;
    grandavg.summary.wrtChanNearHistPeak.SigFraction  = grandavg.summary.wrtChanNearHistPeak.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtChanNearHistPeak.NumTrials,[1 1 1 size(grandavg.summary.wrtChanNearHistPeak.NumSigTrials,4)]);
    grandavg.summary.wrtChanNearHistPeak.SigFraction(isnan(grandavg.summary.wrtChanNearHistPeak.SigFraction)) = 0;
    warning on
    results = grandavg;
    width   = 2;
    SubjectColors(s) = 'k';
  end
  params  = results.params;
  % Plot summary
  lims        = [-1 1 0 1];
  % FIG 1
  if s==1
    fig1      = figure('Name',sprintf('Delay vs distance - fraction of trials with p<%g AND R>Rthresh',params.corrcoef_alpha));
    grandavg  = results;
  elseif s <= length(allresults)
    figure(fig1);
    grandavg.summary.wrtFirstDet.NumTrials            = grandavg.summary.wrtFirstDet.NumTrials            + results.summary.wrtFirstDet.NumTrials;
    grandavg.summary.wrtFirstDet.NumSigTrials         = grandavg.summary.wrtFirstDet.NumSigTrials         + results.summary.wrtFirstDet.NumSigTrials;
    grandavg.summary.wrtChanNearHistPeak.NumTrials    = grandavg.summary.wrtChanNearHistPeak.NumTrials    + results.summary.wrtChanNearHistPeak.NumTrials;
    grandavg.summary.wrtChanNearHistPeak.NumSigTrials = grandavg.summary.wrtChanNearHistPeak.NumSigTrials + results.summary.wrtChanNearHistPeak.NumSigTrials;
  else
    figure(fig1);
  end
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  Npostrl     = squeeze(sum(res.NumTrials(:,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,:),1));  % nlevels
  subplot(3,3,1),plot(Rlevel,Nsigpostrl/Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl/Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); hold on
  ylabel('fraction of trials'); title('wrtFirstChan: peaks (all methods)'); legend('peak1','peak2'); axis(lims); vline(0,'k');
  subplot(3,3,3),plot(Rlevel,(Nsigpostrl+Nsignegtrl)/(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials,2)); vline(0,'k');    % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials,2)); % nmethods x nlevels
  subplot(3,3,2),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
  title('methods (both peaks)'); legend('noflip','flipmatrix','histflip','consistency'); axis(lims); vline(0,'k');hold on
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  % posDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,1),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,1,:),1));  % nlevels
  subplot(3,3,4),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  ylabel('fraction of trials'); title('wrtChanNearHistPeak (posDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,6),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials(:,:,1),2)); vline(0,'k');         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,1,:),2));    % nmethods x nlevels
  subplot(3,3,5),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  title('methods'); axis(lims); vline(0,'k'); hold on% legend('noflip','flipmatrix','histflip','consistency');    
  % negDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,2),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,2,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,2,:),1));  % nlevels
  subplot(3,3,7),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  xlabel('Rthresh'); ylabel('fraction of trials'); title('wrtChanNearHistPeak (negDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,9),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); axis(lims); title('all trials'); xlabel('Rthresh');
  MethN     = squeeze(sum(res.NumTrials(:,:,2),2)); vline(0,'k');hold on         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,2,:),2));    % nmethods x nlevels
  subplot(3,3,8),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  xlabel('Rthresh'); title('methods'); axis(lims); vline(0,'k');hold on % legend('noflip','flipmatrix','histflip','consistency');   
set(gcf,'name','run1');
  % FIG 2
  if s==1
    fig2 = figure('Name',sprintf('Comparison of methods (p<%g AND R>Rthresh)',params.corrcoef_alpha));
  else
    figure(fig2);
  end    
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  subplot(3,4,1),plot(Rlevel,squeeze(res.SigFraction(1,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtFirstDet] noflip'); legend('peak1','peak2'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,2),plot(Rlevel,squeeze(res.SigFraction(2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,3),plot(Rlevel,squeeze(res.SigFraction(3,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,4),plot(Rlevel,squeeze(res.SigFraction(4,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  subplot(3,4,5),plot(Rlevel,squeeze(res.SigFraction(1,1,1,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, posDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,6),plot(Rlevel,squeeze(res.SigFraction(2,1,1,:)), [SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,7),plot(Rlevel,squeeze(res.SigFraction(3,1,1,:)), [SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,8),plot(Rlevel,squeeze(res.SigFraction(4,1,1,:)), [SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,9),plot(Rlevel,squeeze(res.SigFraction(1,1,2,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, negDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,10),plot(Rlevel,squeeze(res.SigFraction(2,1,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,11),plot(Rlevel,squeeze(res.SigFraction(3,1,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,12),plot(Rlevel,squeeze(res.SigFraction(4,1,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');

end
set(gcf,'name','run1');
%%
allresults   = allresults2;
SubjectColors = 'bgrcmyk';
MethodSymbols = {'.','o','x','*'};
PeakLines     = {'--',':'};

% LineType = [SubjectColors(s) MethodSymbols{} PeakLines{}];
clear grandavg fig1 fig2
for s = 1:length(allresults)+1
  if s <= length(allresults)
    results = allresults(s);
    width   = .5;
  else
    warning off
    grandavg.summary.wrtFirstDet.SigFraction          = grandavg.summary.wrtFirstDet.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtFirstDet.NumTrials,[1 1 size(grandavg.summary.wrtFirstDet.NumSigTrials,3)]);
    grandavg.summary.wrtFirstDet.SigFraction(isnan(grandavg.summary.wrtFirstDet.SigFraction)) = 0;
    grandavg.summary.wrtChanNearHistPeak.SigFraction  = grandavg.summary.wrtChanNearHistPeak.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtChanNearHistPeak.NumTrials,[1 1 1 size(grandavg.summary.wrtChanNearHistPeak.NumSigTrials,4)]);
    grandavg.summary.wrtChanNearHistPeak.SigFraction(isnan(grandavg.summary.wrtChanNearHistPeak.SigFraction)) = 0;
    warning on
    results = grandavg;
    width   = 2;
    SubjectColors(s) = 'k';
  end
  params  = results.params;
  % Plot summary
  lims        = [-1 1 0 1];
  % FIG 1
  if s==1
    fig1      = figure('Name',sprintf('Delay vs distance - fraction of trials with p<%g AND R>Rthresh',params.corrcoef_alpha));
    grandavg  = results;
  elseif s <= length(allresults)
    figure(fig1);
    grandavg.summary.wrtFirstDet.NumTrials            = grandavg.summary.wrtFirstDet.NumTrials            + results.summary.wrtFirstDet.NumTrials;
    grandavg.summary.wrtFirstDet.NumSigTrials         = grandavg.summary.wrtFirstDet.NumSigTrials         + results.summary.wrtFirstDet.NumSigTrials;
    grandavg.summary.wrtChanNearHistPeak.NumTrials    = grandavg.summary.wrtChanNearHistPeak.NumTrials    + results.summary.wrtChanNearHistPeak.NumTrials;
    grandavg.summary.wrtChanNearHistPeak.NumSigTrials = grandavg.summary.wrtChanNearHistPeak.NumSigTrials + results.summary.wrtChanNearHistPeak.NumSigTrials;
  else
    figure(fig1);
  end
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  Npostrl     = squeeze(sum(res.NumTrials(:,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,:),1));  % nlevels
  subplot(3,3,1),plot(Rlevel,Nsigpostrl/Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl/Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); hold on
  ylabel('fraction of trials'); title('wrtFirstChan: peaks (all methods)'); legend('peak1','peak2'); axis(lims); vline(0,'k');
  subplot(3,3,3),plot(Rlevel,(Nsigpostrl+Nsignegtrl)/(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials,2)); vline(0,'k');    % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials,2)); % nmethods x nlevels
  subplot(3,3,2),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
  title('methods (both peaks)'); legend('noflip','flipmatrix','histflip','consistency'); axis(lims); vline(0,'k');hold on
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  % posDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,1),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,1,:),1));  % nlevels
  subplot(3,3,4),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  ylabel('fraction of trials'); title('wrtChanNearHistPeak (posDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,6),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials(:,:,1),2)); vline(0,'k');         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,1,:),2));    % nmethods x nlevels
  subplot(3,3,5),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  title('methods'); axis(lims); vline(0,'k'); hold on% legend('noflip','flipmatrix','histflip','consistency');    
  % negDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,2),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,2,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,2,:),1));  % nlevels
  subplot(3,3,7),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  xlabel('Rthresh'); ylabel('fraction of trials'); title('wrtChanNearHistPeak (negDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,9),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); axis(lims); title('all trials'); xlabel('Rthresh');
  MethN     = squeeze(sum(res.NumTrials(:,:,2),2)); vline(0,'k');hold on         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,2,:),2));    % nmethods x nlevels
  subplot(3,3,8),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  xlabel('Rthresh'); title('methods'); axis(lims); vline(0,'k');hold on % legend('noflip','flipmatrix','histflip','consistency');   
set(gcf,'name','run2');
  % FIG 2
  if s==1
    fig2 = figure('Name',sprintf('Comparison of methods (p<%g AND R>Rthresh)',params.corrcoef_alpha));
  else
    figure(fig2);
  end    
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  subplot(3,4,1),plot(Rlevel,squeeze(res.SigFraction(1,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtFirstDet] noflip'); legend('peak1','peak2'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,2),plot(Rlevel,squeeze(res.SigFraction(2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,3),plot(Rlevel,squeeze(res.SigFraction(3,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,4),plot(Rlevel,squeeze(res.SigFraction(4,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  subplot(3,4,5),plot(Rlevel,squeeze(res.SigFraction(1,1,1,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, posDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,6),plot(Rlevel,squeeze(res.SigFraction(2,1,1,:)), [SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,7),plot(Rlevel,squeeze(res.SigFraction(3,1,1,:)), [SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,8),plot(Rlevel,squeeze(res.SigFraction(4,1,1,:)), [SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,9),plot(Rlevel,squeeze(res.SigFraction(1,1,2,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, negDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,10),plot(Rlevel,squeeze(res.SigFraction(2,1,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,11),plot(Rlevel,squeeze(res.SigFraction(3,1,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,12),plot(Rlevel,squeeze(res.SigFraction(4,1,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');

end
set(gcf,'name','run2');

%%

allresults   = allresults3;
SubjectColors = 'bgrcmyk';
MethodSymbols = {'.','o','x','*'};
PeakLines     = {'--',':'};

% LineType = [SubjectColors(s) MethodSymbols{} PeakLines{}];
clear grandavg fig1 fig2
for s = 1:length(allresults)+1
  if s <= length(allresults)
    results = allresults(s);
    width   = .5;
  else
    warning off
    grandavg.summary.wrtFirstDet.SigFraction          = grandavg.summary.wrtFirstDet.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtFirstDet.NumTrials,[1 1 size(grandavg.summary.wrtFirstDet.NumSigTrials,3)]);
    grandavg.summary.wrtFirstDet.SigFraction(isnan(grandavg.summary.wrtFirstDet.SigFraction)) = 0;
    grandavg.summary.wrtChanNearHistPeak.SigFraction  = grandavg.summary.wrtChanNearHistPeak.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtChanNearHistPeak.NumTrials,[1 1 1 size(grandavg.summary.wrtChanNearHistPeak.NumSigTrials,4)]);
    grandavg.summary.wrtChanNearHistPeak.SigFraction(isnan(grandavg.summary.wrtChanNearHistPeak.SigFraction)) = 0;
    warning on
    results = grandavg;
    width   = 2;
    SubjectColors(s) = 'k';
  end
  params  = results.params;
  % Plot summary
  lims        = [-1 1 0 1];
  % FIG 1
  if s==1
    fig1      = figure('Name',sprintf('Delay vs distance - fraction of trials with p<%g AND R>Rthresh',params.corrcoef_alpha));
    grandavg  = results;
  elseif s <= length(allresults)
    figure(fig1);
    grandavg.summary.wrtFirstDet.NumTrials            = grandavg.summary.wrtFirstDet.NumTrials            + results.summary.wrtFirstDet.NumTrials;
    grandavg.summary.wrtFirstDet.NumSigTrials         = grandavg.summary.wrtFirstDet.NumSigTrials         + results.summary.wrtFirstDet.NumSigTrials;
    grandavg.summary.wrtChanNearHistPeak.NumTrials    = grandavg.summary.wrtChanNearHistPeak.NumTrials    + results.summary.wrtChanNearHistPeak.NumTrials;
    grandavg.summary.wrtChanNearHistPeak.NumSigTrials = grandavg.summary.wrtChanNearHistPeak.NumSigTrials + results.summary.wrtChanNearHistPeak.NumSigTrials;
  else
    figure(fig1);
  end
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  Npostrl     = squeeze(sum(res.NumTrials(:,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,:),1));  % nlevels
  subplot(3,3,1),plot(Rlevel,Nsigpostrl/Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl/Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); hold on
  ylabel('fraction of trials'); title('wrtFirstChan: peaks (all methods)'); legend('peak1','peak2'); axis(lims); vline(0,'k');
  subplot(3,3,3),plot(Rlevel,(Nsigpostrl+Nsignegtrl)/(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials,2)); vline(0,'k');    % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials,2)); % nmethods x nlevels
  subplot(3,3,2),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
  title('methods (both peaks)'); legend('noflip','flipmatrix','histflip','consistency'); axis(lims); vline(0,'k');hold on
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  % posDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,1),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,1,:),1));  % nlevels
  subplot(3,3,4),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  ylabel('fraction of trials'); title('wrtChanNearHistPeak (posDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,6),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials(:,:,1),2)); vline(0,'k');         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,1,:),2));    % nmethods x nlevels
  subplot(3,3,5),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  title('methods'); axis(lims); vline(0,'k'); hold on% legend('noflip','flipmatrix','histflip','consistency');    
  % negDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,2),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,2,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,2,:),1));  % nlevels
  subplot(3,3,7),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  xlabel('Rthresh'); ylabel('fraction of trials'); title('wrtChanNearHistPeak (negDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,9),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); axis(lims); title('all trials'); xlabel('Rthresh');
  MethN     = squeeze(sum(res.NumTrials(:,:,2),2)); vline(0,'k');hold on         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,2,:),2));    % nmethods x nlevels
  subplot(3,3,8),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  xlabel('Rthresh'); title('methods'); axis(lims); vline(0,'k');hold on % legend('noflip','flipmatrix','histflip','consistency');   
set(gcf,'name','run3');
  % FIG 2
  if s==1
    fig2 = figure('Name',sprintf('Comparison of methods (p<%g AND R>Rthresh)',params.corrcoef_alpha));
  else
    figure(fig2);
  end    
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  subplot(3,4,1),plot(Rlevel,squeeze(res.SigFraction(1,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtFirstDet] noflip'); legend('peak1','peak2'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,2),plot(Rlevel,squeeze(res.SigFraction(2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,3),plot(Rlevel,squeeze(res.SigFraction(3,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,4),plot(Rlevel,squeeze(res.SigFraction(4,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  subplot(3,4,5),plot(Rlevel,squeeze(res.SigFraction(1,1,1,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, posDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,6),plot(Rlevel,squeeze(res.SigFraction(2,1,1,:)), [SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,7),plot(Rlevel,squeeze(res.SigFraction(3,1,1,:)), [SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,8),plot(Rlevel,squeeze(res.SigFraction(4,1,1,:)), [SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,9),plot(Rlevel,squeeze(res.SigFraction(1,1,2,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, negDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,10),plot(Rlevel,squeeze(res.SigFraction(2,1,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,11),plot(Rlevel,squeeze(res.SigFraction(3,1,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,12),plot(Rlevel,squeeze(res.SigFraction(4,1,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');

end
set(gcf,'name','run3');


%%

allresults   = allresults4;
SubjectColors = 'bgrcmyk';
MethodSymbols = {'.','o','x','*'};
PeakLines     = {'--',':'};

% LineType = [SubjectColors(s) MethodSymbols{} PeakLines{}];
clear grandavg fig1 fig2
for s = 1:length(allresults)+1
  if s <= length(allresults)
    results = allresults(s);
    width   = .5;
  else
    warning off
    grandavg.summary.wrtFirstDet.SigFraction          = grandavg.summary.wrtFirstDet.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtFirstDet.NumTrials,[1 1 size(grandavg.summary.wrtFirstDet.NumSigTrials,3)]);
    grandavg.summary.wrtFirstDet.SigFraction(isnan(grandavg.summary.wrtFirstDet.SigFraction)) = 0;
    grandavg.summary.wrtChanNearHistPeak.SigFraction  = grandavg.summary.wrtChanNearHistPeak.NumSigTrials ...
      ./ repmat(grandavg.summary.wrtChanNearHistPeak.NumTrials,[1 1 1 size(grandavg.summary.wrtChanNearHistPeak.NumSigTrials,4)]);
    grandavg.summary.wrtChanNearHistPeak.SigFraction(isnan(grandavg.summary.wrtChanNearHistPeak.SigFraction)) = 0;
    warning on
    results = grandavg;
    width   = 2;
    SubjectColors(s) = 'k';
  end
  params  = results.params;
  % Plot summary
  lims        = [-1 1 0 1];
  % FIG 1
  if s==1
    fig1      = figure('Name',sprintf('Delay vs distance - fraction of trials with p<%g AND R>Rthresh',params.corrcoef_alpha));
    grandavg  = results;
  elseif s <= length(allresults)
    figure(fig1);
    grandavg.summary.wrtFirstDet.NumTrials            = grandavg.summary.wrtFirstDet.NumTrials            + results.summary.wrtFirstDet.NumTrials;
    grandavg.summary.wrtFirstDet.NumSigTrials         = grandavg.summary.wrtFirstDet.NumSigTrials         + results.summary.wrtFirstDet.NumSigTrials;
    grandavg.summary.wrtChanNearHistPeak.NumTrials    = grandavg.summary.wrtChanNearHistPeak.NumTrials    + results.summary.wrtChanNearHistPeak.NumTrials;
    grandavg.summary.wrtChanNearHistPeak.NumSigTrials = grandavg.summary.wrtChanNearHistPeak.NumSigTrials + results.summary.wrtChanNearHistPeak.NumSigTrials;
  else
    figure(fig1);
  end
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  Npostrl     = squeeze(sum(res.NumTrials(:,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,:),1));  % nlevels
  subplot(3,3,1),plot(Rlevel,Nsigpostrl/Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl/Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); hold on
  ylabel('fraction of trials'); title('wrtFirstChan: peaks (all methods)'); legend('peak1','peak2'); axis(lims); vline(0,'k');
  subplot(3,3,3),plot(Rlevel,(Nsigpostrl+Nsignegtrl)/(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials,2)); vline(0,'k');    % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials,2)); % nmethods x nlevels
  subplot(3,3,2),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
  title('methods (both peaks)'); legend('noflip','flipmatrix','histflip','consistency'); axis(lims); vline(0,'k');hold on
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  % posDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,1),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,1),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,1,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,1,:),1));  % nlevels
  subplot(3,3,4),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  ylabel('fraction of trials'); title('wrtChanNearHistPeak (posDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,6),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); title('all trials'); axis(lims);hold on
  MethN     = squeeze(sum(res.NumTrials(:,:,1),2)); vline(0,'k');         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,1,:),2));    % nmethods x nlevels
  subplot(3,3,5),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  title('methods'); axis(lims); vline(0,'k'); hold on% legend('noflip','flipmatrix','histflip','consistency');    
  % negDelay
  Npostrl     = squeeze(sum(res.NumTrials(:,1,2),1));       % nlevels
  Nnegtrl     = squeeze(sum(res.NumTrials(:,2,2),1));       % nlevels
  Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,2,:),1));  % nlevels
  Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,2,:),1));  % nlevels
  subplot(3,3,7),plot(Rlevel,Nsigpostrl./Npostrl,[SubjectColors(s) '.' PeakLines{1}],Rlevel,Nsignegtrl./Nnegtrl,[SubjectColors(s) '.' PeakLines{2}],'linewidth',width); axis(lims);hold on
  xlabel('Rthresh'); ylabel('fraction of trials'); title('wrtChanNearHistPeak (negDelay)'); vline(0,'k'); % legend('peak1','peak2');
  subplot(3,3,9),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),[SubjectColors(s) '.-'],'linewidth',width); axis(lims); title('all trials'); xlabel('Rthresh');
  MethN     = squeeze(sum(res.NumTrials(:,:,2),2)); vline(0,'k');hold on         % nmethods x 1
  MethSigN  = squeeze(sum(res.NumSigTrials(:,:,2,:),2));    % nmethods x nlevels
  subplot(3,3,8),plot(Rlevel,MethSigN(1,:)/MethN(1),[SubjectColors(s) MethodSymbols{1}],Rlevel,MethSigN(2,:)/MethN(2),[SubjectColors(s) MethodSymbols{2}],Rlevel,MethSigN(3,:)/MethN(3),[SubjectColors(s) MethodSymbols{3}],Rlevel,MethSigN(4,:)/MethN(4),[SubjectColors(s) MethodSymbols{4}],'linewidth',width);
%     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
  xlabel('Rthresh'); title('methods'); axis(lims); vline(0,'k');hold on % legend('noflip','flipmatrix','histflip','consistency');   
set(gcf,'name','run4');
  % FIG 2
  if s==1
    fig2 = figure('Name',sprintf('Comparison of methods (p<%g AND R>Rthresh)',params.corrcoef_alpha));
  else
    figure(fig2);
  end    
  res         = results.summary.wrtFirstDet;
  Rlevel      = res.Rlevels;
  subplot(3,4,1),plot(Rlevel,squeeze(res.SigFraction(1,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtFirstDet] noflip'); legend('peak1','peak2'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,2),plot(Rlevel,squeeze(res.SigFraction(2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,3),plot(Rlevel,squeeze(res.SigFraction(3,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,4),plot(Rlevel,squeeze(res.SigFraction(4,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  res         = results.summary.wrtChanNearHistPeak;
  Rlevel      = res.Rlevels;
  subplot(3,4,5),plot(Rlevel,squeeze(res.SigFraction(1,1,1,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,1,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, posDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,6),plot(Rlevel,squeeze(res.SigFraction(2,1,1,:)), [SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,1,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,7),plot(Rlevel,squeeze(res.SigFraction(3,1,1,:)), [SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,1,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,8),plot(Rlevel,squeeze(res.SigFraction(4,1,1,:)), [SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,1,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,9),plot(Rlevel,squeeze(res.SigFraction(1,1,2,:)), [SubjectColors(s) MethodSymbols{1} PeakLines{1}],Rlevel,squeeze(res.SigFraction(1,2,2,:)),[SubjectColors(s) MethodSymbols{1} PeakLines{2}],'linewidth',width); title('[wrtChanNearHistPeak, negDelay] noflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,10),plot(Rlevel,squeeze(res.SigFraction(2,1,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{1}],Rlevel,squeeze(res.SigFraction(2,2,2,:)),[SubjectColors(s) MethodSymbols{2} PeakLines{2}],'linewidth',width); title('flipmatrix'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,11),plot(Rlevel,squeeze(res.SigFraction(3,1,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{1}],Rlevel,squeeze(res.SigFraction(3,2,2,:)),[SubjectColors(s) MethodSymbols{3} PeakLines{2}],'linewidth',width); title('histflip'); axis(lims);hold on; vline(0,'k');
  subplot(3,4,12),plot(Rlevel,squeeze(res.SigFraction(4,1,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{1}],Rlevel,squeeze(res.SigFraction(4,2,2,:)),[SubjectColors(s) MethodSymbols{4} PeakLines{2}],'linewidth',width); title('consistency'); axis(lims);hold on; vline(0,'k');

end
set(gcf,'name','run4');