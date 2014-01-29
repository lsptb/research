function varargout = plot_intervals(SubjectList,rootdir,CorrCoefFile)

if nargin <= 1 || ~ischar(rootdir)
  rootdir = '/space/emc2/1/halgdev/projects/sleep/MEG/SO';
elseif nargin <=2 || ~ischar(CorrCoefFile)
  CorrCoefFile = 'matfiles/SO_corrcoef.mat';
end
if ~isnumeric(SubjectList), error('Input must be a (subject) number.'); end
cwd     = pwd;
cnt     = 1;
ns      = length(SubjectList);
figure('Name','Detection Count'); 
for s = 1:ns
  SubjID  = sprintf('s%g',SubjectList(s));
  SubjDir = [rootdir '/' SubjID];
  if ~exist(SubjDir,'dir'), continue; end
  fprintf('Processing subject %g of %g (s%g)\n',s,length(SubjectList),SubjectList(s));
  cd(SubjDir);
  file = 'matfiles/SO_clustered_detections.mat';
  if exist(file,'file')
    load(file);
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
    nn    = round((max(tc)-min(tc))/params.AASM_EpochLength);
    [N,X] = hist(tc,nn);
    tmpx  = zeros(1,length(t));
    tmpx(clusterIX) = 1;
    subplot(4,ns,s+0*ns),plot(t,tmpx,'.-'); axis tight; title(sprintf('%s: aggregate SO detection count',SubjID))
    subplot(4,ns,s+1*ns),try hist(tc,nn); end; axis tight
    subplot(4,ns,s+2*ns),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
    xlabel('time (sec)'); ylabel(sprintf('count (%gsec bins)',params.AASM_EpochLength));
    for k=1:length(t0f), vline(t0f(k),'k'); end
    subplot(4,ns,s+3*ns),plot(tc,'.-'); ylabel('cluster time (sec)'); xlabel('cluster number'); axis tight;
    clear nn N X tmpx t0f
    % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end    
%   file = CorrCoefFile;
%   if exist(file,'file')
%     load(file);
%     if ~isfield(results,'summary')
%       % Calculate summary measures of sigificant correlations for this subject
%       alpha     = params.corrcoef_alpha;
%       Rlevels   = 0:.1:1;
%       % First RefType
%       reftype   = 'wrtFirstDet';
%       methods   = {'noflip','flipmatrix','histpeak','consistency'};
%       peaktypes = {'pospeak','negpeak'};
%       trlfrac   = zeros(length(methods),length(peaktypes),length(Rlevels));
%       Nchan     = zeros(size(trlfrac));
%       Nsigtrl   = zeros(size(trlfrac));
%       Ntrl      = zeros(length(methods),length(peaktypes));
%       for m = 1:length(methods)
%         method = methods{m};
%         for pk = 1:length(peaktypes)
%           peaktype  = peaktypes{pk};
%           R         = results.(method).(peaktype).(reftype).R;
%           p         = results.(method).(peaktype).(reftype).p;
%           N         = results.(method).(peaktype).(reftype).N;
%           if isempty(R), continue; end
%           % trials with p < alpha
%           Ntrials     = length(R);
%           sigtrials   = cellfun(@(x)find((p<alpha)&(R>x)),num2cell(Rlevels),'UniformOutput',false);
%           Nsigtrials  = cellfun(@length,sigtrials);
%           if length(Nsigtrials) < length(Rlevels)
%             error('Some R fraction values were lost.');
%           end
%           Ntrl(m,pk)       = Ntrials;
%           Nsigtrl(m,pk,:)  = Nsigtrials;
%           trlfrac(m,pk,:)  = Nsigtrials / Ntrials;
%           Nchan(m,pk,:)    = cellfun(@(x)sum(N(x)),sigtrials);
%         end
%       end
%       results.summary.(reftype).methods       = methods;
%       results.summary.(reftype).peaktypes     = peaktypes;
%       results.summary.(reftype).Rlevels       = Rlevels;
%       results.summary.(reftype).dimorder      = {'methods','peaktypes','Rlevels'};
%       results.summary.(reftype).NumTrials     = Ntrl;
%       results.summary.(reftype).NumSigTrials  = Nsigtrl;
%       results.summary.(reftype).SigFraction   = trlfrac;
%       results.summary.(reftype).NumInvChans   = Nchan;
%       % Second RefType
%       reftype   = 'wrtChanNearHistPeak';
%       methods   = {'noflip','flipmatrix','histpeak','consistency'};
%       peaktypes = {'pospeak','negpeak'};
%       delaystr  = {'posDelay','negDelay'};
%       trlfrac   = zeros(length(methods),length(peaktypes),length(delaystr),length(Rlevels));
%       Nchan     = zeros(size(trlfrac));
%       Nsigtrl   = zeros(size(trlfrac));
%       Ntrl      = zeros(length(methods),length(peaktypes),length(delaystr));
%       for m = 1:length(methods)
%         method = methods{m};
%         for pk = 1:length(peaktypes)
%           peaktype = peaktypes{pk};
%           for d  = 1:length(delaystr)
%             suff = ['_' delaystr{d}];
%             R       = results.(method).(peaktype).(reftype).(['R' suff]);
%             p       = results.(method).(peaktype).(reftype).(['p' suff]);
%             N       = results.(method).(peaktype).(reftype).(['N' suff]);
%             if isempty(R), continue; end
%             % trials with p < alpha
%             Ntrials     = length(R);
%             sigtrials   = cellfun(@(x)find((p<alpha)&(R>x)),num2cell(Rlevels),'UniformOutput',false);
%             Nsigtrials  = cellfun(@length,sigtrials);
%             if length(Nsigtrials) < length(Rlevels)
%               error('Some R fraction values were lost.');
%             end
%             Ntrl(m,pk,d)       = Ntrials;
%             Nsigtrl(m,pk,d,:)  = Nsigtrials;
%             trlfrac(m,pk,d,:)  = Nsigtrials / Ntrials;
%             Nchan(m,pk,d,:)    = cellfun(@(x)sum(N(x)),sigtrials);
%           end
%         end
%       end
%       results.summary.(reftype).methods       = methods;
%       results.summary.(reftype).peaktypes     = peaktypes;
%       results.summary.(reftype).delaytypes    = delaystr;
%       results.summary.(reftype).Rlevels       = Rlevels;
%       results.summary.(reftype).dimorder      = {'methods','peaktypes','delaytypes','Rlevels'};
%       results.summary.(reftype).NumTrials     = Ntrl;
%       results.summary.(reftype).NumSigTrials  = Nsigtrl;
%       results.summary.(reftype).SigFraction   = trlfrac;
%       results.summary.(reftype).NumInvChans   = Nchan;
%       fprintf(fid,'Adding summary to results of the cluster correlation coefficient calculations\n');
%       save(file,'results','params');
%       telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);  
%     end
%     % Plot summary
%     lims        = [0 1 0 .45];
%     % FIG 1
%     figure('Name',sprintf('Subject %s: delay vs distance - fraction of trials with p<%g AND R>Rthresh',params.SubjID,params.corrcoef_alpha));
%     res         = results.summary.wrtFirstDet;
%     Rlevel      = res.Rlevels;
%     Npostrl     = squeeze(sum(res.NumTrials(:,1),1));       % nlevels
%     Nnegtrl     = squeeze(sum(res.NumTrials(:,2),1));       % nlevels
%     Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,:),1));  % nlevels
%     Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,:),1));  % nlevels
%     subplot(3,3,1),plot(Rlevel,Nsigpostrl/Npostrl,'b.',Rlevel,Nsignegtrl/Nnegtrl,'r.');
%     ylabel('fraction of trials'); title('wrtFirstChan: peaks (all methods)'); legend('peak1','peak2'); axis(lims);
%     subplot(3,3,3),plot(Rlevel,(Nsigpostrl+Nsignegtrl)/(Npostrl+Nnegtrl),'k.'); title('all trials'); axis(lims);
%     MethN     = squeeze(sum(res.NumTrials,2));    % nmethods x 1
%     MethSigN  = squeeze(sum(res.NumSigTrials,2)); % nmethods x nlevels
%     subplot(3,3,2),plot(Rlevel,MethSigN(1,:)/MethN(1),'r*',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g^',Rlevel,MethSigN(4,:)/MethN(4),'ko');
%     title('methods (both peaks)'); legend('noflip','flipmatrix','histflip','consistency'); axis(lims);
%     res         = results.summary.wrtChanNearHistPeak;
%     Rlevel      = res.Rlevels;
%     % posDelay
%     Npostrl     = squeeze(sum(res.NumTrials(:,1,1),1));       % nlevels
%     Nnegtrl     = squeeze(sum(res.NumTrials(:,2,1),1));       % nlevels
%     Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,1,:),1));  % nlevels
%     Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,1,:),1));  % nlevels
%     subplot(3,3,4),plot(Rlevel,Nsigpostrl./Npostrl,'b.',Rlevel,Nsignegtrl./Nnegtrl,'r.'); axis(lims);
%     ylabel('fraction of trials'); title('wrtChanNearHistPeak (posDelay)'); % legend('peak1','peak2');
%     subplot(3,3,6),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),'k.'); title('all trials'); axis(lims);
%     MethN     = squeeze(sum(res.NumTrials(:,:,1),2));         % nmethods x 1
%     MethSigN  = squeeze(sum(res.NumSigTrials(:,:,1,:),2));    % nmethods x nlevels
%     subplot(3,3,5),plot(Rlevel,MethSigN(1,:)/MethN(1),'r*',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g^',Rlevel,MethSigN(4,:)/MethN(4),'ko');
% %     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
%     title('methods'); axis(lims); % legend('noflip','flipmatrix','histflip','consistency');    
%     % negDelay
%     Npostrl     = squeeze(sum(res.NumTrials(:,1,2),1));       % nlevels
%     Nnegtrl     = squeeze(sum(res.NumTrials(:,2,2),1));       % nlevels
%     Nsigpostrl  = squeeze(sum(res.NumSigTrials(:,1,2,:),1));  % nlevels
%     Nsignegtrl  = squeeze(sum(res.NumSigTrials(:,2,2,:),1));  % nlevels
%     subplot(3,3,7),plot(Rlevel,Nsigpostrl./Npostrl,'b.',Rlevel,Nsignegtrl./Nnegtrl,'r.'); axis(lims);
%     xlabel('Rthresh'); ylabel('fraction of trials'); title('wrtChanNearHistPeak (negDelay)'); % legend('peak1','peak2');
%     subplot(3,3,9),plot(Rlevel,(Nsigpostrl+Nsignegtrl)./(Npostrl+Nnegtrl),'k.'); axis(lims); title('all trials'); xlabel('Rthresh');
%     MethN     = squeeze(sum(res.NumTrials(:,:,2),2));         % nmethods x 1
%     MethSigN  = squeeze(sum(res.NumSigTrials(:,:,2,:),2));    % nmethods x nlevels
%     subplot(3,3,8),plot(Rlevel,MethSigN(1,:)/MethN(1),'r*',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g^',Rlevel,MethSigN(4,:)/MethN(4),'ko');
% %     plot(Rlevel,MethSigN(1,:)/MethN(1),'r.',Rlevel,MethSigN(2,:)/MethN(2),'b.',Rlevel,MethSigN(3,:)/MethN(3),'g.',Rlevel,MethSigN(4,:)/MethN(4),'k.');
%     xlabel('Rthresh'); title('methods'); axis(lims); % legend('noflip','flipmatrix','histflip','consistency');   
% 
%     % FIG 2
%     figure('Name',sprintf('Subject %s: comparison of methods (p<%g AND R>Rthresh)',params.SubjID,params.corrcoef_alpha));
%     res         = results.summary.wrtFirstDet;
%     Rlevel      = res.Rlevels;
%     subplot(3,4,1),plot(Rlevel,squeeze(res.SigFraction(1,1,:)),'b*',Rlevel,squeeze(res.SigFraction(1,2,:)),'r*'); title('[wrtFirstDet] noflip'); legend('peak1','peak2'); axis(lims);
%     subplot(3,4,2),plot(Rlevel,squeeze(res.SigFraction(2,1,:)),'b.',Rlevel,squeeze(res.SigFraction(2,2,:)),'r.'); title('flipmatrix'); axis(lims);
%     subplot(3,4,3),plot(Rlevel,squeeze(res.SigFraction(3,1,:)),'b^',Rlevel,squeeze(res.SigFraction(3,2,:)),'r^'); title('histflip'); axis(lims);
%     subplot(3,4,4),plot(Rlevel,squeeze(res.SigFraction(4,1,:)),'bo',Rlevel,squeeze(res.SigFraction(4,2,:)),'ro'); title('consistency'); axis(lims);
%     res         = results.summary.wrtChanNearHistPeak;
%     Rlevel      = res.Rlevels;
%     subplot(3,4,5),plot(Rlevel,squeeze(res.SigFraction(1,1,1,:)),'b*',Rlevel,squeeze(res.SigFraction(1,2,1,:)),'r*'); title('[wrtChanNearHistPeak, posDelay] noflip'); axis(lims);
%     subplot(3,4,6),plot(Rlevel,squeeze(res.SigFraction(2,1,1,:)),'b.',Rlevel,squeeze(res.SigFraction(2,2,1,:)),'r.'); title('flipmatrix'); axis(lims);
%     subplot(3,4,7),plot(Rlevel,squeeze(res.SigFraction(3,1,1,:)),'b^',Rlevel,squeeze(res.SigFraction(3,2,1,:)),'r^'); title('histflip'); axis(lims);
%     subplot(3,4,8),plot(Rlevel,squeeze(res.SigFraction(4,1,1,:)),'bo',Rlevel,squeeze(res.SigFraction(4,2,1,:)),'ro'); title('consistency'); axis(lims);
%     subplot(3,4,9),plot(Rlevel,squeeze(res.SigFraction(1,1,2,:)),'b*',Rlevel,squeeze(res.SigFraction(1,2,2,:)),'r*'); title('[wrtChanNearHistPeak, negDelay] noflip'); axis(lims);
%     subplot(3,4,10),plot(Rlevel,squeeze(res.SigFraction(2,1,2,:)),'b.',Rlevel,squeeze(res.SigFraction(2,2,2,:)),'r.'); title('flipmatrix'); axis(lims);
%     subplot(3,4,11),plot(Rlevel,squeeze(res.SigFraction(3,1,2,:)),'b^',Rlevel,squeeze(res.SigFraction(3,2,2,:)),'r^'); title('histflip'); axis(lims);
%     subplot(3,4,12),plot(Rlevel,squeeze(res.SigFraction(4,1,2,:)),'bo',Rlevel,squeeze(res.SigFraction(4,2,2,:)),'ro'); title('consistency'); axis(lims);
%   end
%   if ~isfield(results,'params')
%     results.params = params;
%   end
%   out(cnt)  = results;
%   cnt       = cnt + 1;
end  
cd(cwd);
if nargout > 0, varargout{1} = out; end