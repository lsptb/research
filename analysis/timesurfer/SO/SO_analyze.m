function varargout = SO_analyze(SubjID)
% Purpose: function to produce single-subject results for SO paper I
%
% Structures:
%   params
%   detections
%   clusters
%
% Functions:
%   params      = SO_params(SubjID)
%   detections  = SO_detection(data,key,value,...)
%   detections  = SO_cluster_detections(detections,key,value,...)
%   clusters    = SO_cluster_corrcoef(data,detections,key,value,...)
%   SO_analyze(SubjID)
%   SO_figures(SubjID)
%   Auxiliary:
%     LoadData(params)
%
% MAT-files (in each subject directory):
%   SO_params.mat
%   SO_detections.mat
%   SO_clustered_detections.mat
%   SO_clustered_consistent_detections.mat
%   SO_clusters_noflip.mat
%   SO_clusters_flipmatrix.mat
%   SO_clusters_consistency.mat
%   SO_histflip.mat
% 
% Created by Jason Sherfey on 27-Jul-2010
% Multimodal Imaging Laboratory
% Department of Radiology, UCSD
%
% Last modified on 30-Sep-2010
% 
% Modified on 03-Aug-2010
% - clusters involving < threshold % of channels are now removed after
% automatic N3 detection instead of before. Removed script IntervalSelection.
%
% This functions uses:
% SO_params.m, LoadData.m, SO_cluster_corrcoef.m, SO_cluster_detections.m,
% SO_detection.m, calc_flip_matrix.m, select_peakpairs.m
%
% Latest versions are located in
% /space/emc2/1/halgdev/projects/sleep/MEG/SO/scripts

% addpath /space/emc2/1/halgdev/projects/sleep/MEG/SO/scripts

% check inputs
if nargin == 0
  SubjID = input('Subject ID:');
end
if length(SubjID) > 1 && ~ischar(SubjID)
  tic
  for s    = 1:length(SubjID)
    fprintf('Processing subject %g of %g\n',s,length(SubjID));
    tic
    res(s) = SO_analyze(SubjID(s));
    toc
  end
  if nargout > 0
    varargout{1} = res;
  end
  return
end

% Prep
strID   = '';
cwd     = pwd;
tstart  = tic;
params  = SO_params(SubjID);
if ~exist(params.SubjDir,'dir')
  mkdir(params.SubjDir);
end
cd(params.SubjDir)

% try
  
  if ~exist('matfiles','dir')
    mkdir('matfiles');
  end

  % Logs
  logfile = sprintf('%s.log',date);
  fid     = fopen(logfile,'a');
  [~,b]   = unix('echo $HOST');
  fprintf(fid,'---------------------------------\n');
  fprintf(fid,'Slow oscillation analysis.\n');
  fprintf(fid,'Subject: %s\n',params.SubjID);
  fprintf(fid,'Subjdir: %s\n',params.SubjDir);
  fprintf(fid,'Time:    %s\n',datestr(now));
  fprintf(fid,'Host:    %s',b);
  fprintf(fid,'---------------------------------\n');
  clear b

  %% LOAD DATA
  LoadData
  telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % Summary: this script returns a "data" structure after using manual ICA 
    % to remove EKG artifacts.
    % Purpose: this script loads the raw data from fif files, converts the data
    % into TimeSurfer format (continuous data in an epoch_data structure with
    % one trial), preprocesses the data and runs ICA to remove EKG.  The script
    % will start with the most recently saved files (ex. the post-ICA
    % epoch_data will be loaded directly if the MAT file already exists).

  %% SO DETECTION
  file = sprintf('matfiles/SO_%s_detections.mat',params.chantype);
  if exist(file,'file')
    fprintf(fid,'Loading SO detections: %s\n',file);
    load(file,'detections');
    if ~isempty(params.toilim)
      fprintf(fid,'Selecting detections in %g-%gsec from detections in %g-%gsec\n',params.toilim,detections(1).tstart,detections(1).tstop);
      t   = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
      tmp = {detections.pospeak};
      tmp = cellfun(@(x)x(t(x)>=params.toilim(1)&t(x)<=params.toilim(2)),tmp,'uniformoutput',false);
      [detections.pospeak] = deal(tmp{:});
      tmp = {detections.negpeak};
      tmp = cellfun(@(x)x(t(x)>=params.toilim(1)&t(x)<=params.toilim(2)),tmp,'uniformoutput',false);
      [detections.negpeak] = deal(tmp{:});
      clear tmp
    end
    tmp1 = cellfun(@length,{detections.pospeak});
    tmp2 = cellfun(@length,{detections.negpeak});
    tmp  = mean([tmp1;tmp2],1);
    fprintf(fid,'Set contains %g detections across all channels.\n',round(sum(tmp)));
    clear tmp tmp1 tmp2
  elseif params.SO_detections_flag
    origparams        = params;
    params.monotonic  = 0;
    args              = mmil_parms2args(params);
    fprintf(fid,'Detecting slow oscillations in %s time series: t = %g - %g sec\n',params.chantype,data.epochs.time(1),data.epochs.time(end));
    if params.bpfilter,         fprintf(fid,'\t%g-%gHz zero phase shift frequency domain bandpass filter.\n',params.bpfreq); end
    if params.blc,              fprintf(fid,'\tRemove mean from data.\n'); end
    if params.decimate,         fprintf(fid,'\tDecimate by factor of %g.\n',params.decimate_factor); end
    if params.smooth,           fprintf(fid,'\tSmooth with a linear fit (lowess) in a %s sec sliding window.\n', params.smooth_window); end
    if params.hilbertpeaks,     fprintf(fid,'\tDetect peaks using the instantaneous Hilbert phase.\n'); end
    if params.derivpeaks,       fprintf(fid,'\tRefine peak detection based on zero-crossings of the signal derivative.\n'); end
    if params.onlysinglepeaks,  fprintf(fid,'\tSkip cycles with multiple peaks.\n'); end
    if params.zerocross,        fprintf(fid,'\tSelect cycles with %g-%g sec between zero crossings.\n',params.zero2zero_limits); end
    if params.monotonic,        fprintf(fid,'\tSelect cycles with monotonic phase runs.\n'); end
    if params.gtmedian,         fprintf(fid,'\tSelect cycles with peak amplitudes greater than the median of detected peaks.\n'); end   
    detections        = SO_detection(data,args{:});
    params            = origparams;
    clear origparams
    if params.peakpairs_flag
      fprintf(fid,'\tSelecting pairs of adjacent half-wave detections.\n');
      all_detections = detections;
      % This will eliminate all positive detections without an adjacent negative
      % peak and vice versa.  detections must occur within tau of each other.
      detections  = select_peakpairs(detections,params.peakpairs_tau);
    end
    t             = data.epochs.time;
    events  = [];
      % The events structure can be used by visualizer to mark detections
    for k   = 1:length(detections)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
      events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
    end
    fprintf(fid,'Saving SO detections: %s\n',file);
    if params.peakpairs_flag
      save(file,'events','detections','params','all_detections');
    else
      save(file,'events','detections','params');
    end
    clear events args t k
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % visualizer(data); % load SO_detections.mat
  end
  clear file tmp
%   if exist('detections','var'), all_detections = detections; end
  
  %% CLUSTERING SO DETECTIONS BASED ON AGGREGATE DETECTION COUNT HISTOGRAM
  file = sprintf('matfiles/SO_clustered_%s_%s_detections.mat',params.chantype,params.peaktype);
  if strcmp(params.chantype,'grad') && ~isempty(params.cluster_count)
    file = sprintf('matfiles/SO_clustered_%s_%s_detections_ref-%s.mat',params.chantype,params.peaktype,params.cluster_count);
  end
  if exist(file,'file')
    fprintf(fid,'Loading SO %s %s detection clusters: %s\n',params.chantype,params.peaktype,file);
    load(file,'events','detections','count','CountTime','ClusterTime','CountThresh','clusterID','clusterIX','Ninvolved');
  elseif params.SO_cluster_detections_flag
    fprintf(fid,'Finding clusters of %s %s detections.\n',params.chantype,params.peaktype);
    fprintf(fid,'\tCount the number of detections across channels in a %gsec window every %gsec.\n',params.cluster_IntegrationWindow,params.cluster_StepSize);
    fprintf(fid,'\tFind time intervals exceeding a count threshold of %g.\n',params.cluster_thresh);
    fprintf(fid,'\tCombine intervals separated by less than %gsec.\n',params.cluster_MinSeparation);
    fprintf(fid,'\tIgnore intervals with durations less than %gsec.\n',params.cluster_MinSeparation);
    fprintf(fid,'\tDefine histogram peaks as centers of the resulting time intervals.\n');
    fprintf(fid,'\tCluster detections in %gsec windows centered on histogram peaks.\n',params.cluster_ClusterWindow);    
    if strcmp(params.chantype,'grad') && strcmp(params.cluster_count,'eeg_negpeak')
      cluster_file = sprintf('matfiles/SO_clustered_eeg_negpeak_detections.mat');
      fprintf(fid,'Clustering MEG data using EEG negpeak clusters as reference.\n');
      if exist(cluster_file,'file')
        fprintf(fid,'EEG cluster file: %s\n',cluster_file);
        % load EEG negative half-wave agg count
        load('matfiles/SO_clustered_eeg_negpeak_detections.mat','count','CountTime');
        count = count(CountTime<detections(1).tstop);
        % use EEG peaks to cluster MEG detections
        [detections,count,CountTime,ClusterTime,CountThresh] = SO_cluster_detections(...
          detections,'method','histogram','threshold',params.cluster_thresh,...
          'StepSize',params.cluster_StepSize,'IntegrationWindow',params.cluster_IntegrationWindow,...
          'MinSeparation',params.cluster_MinSeparation,'ClusterWindow',params.cluster_ClusterWindow,...
          'peaktype',params.peaktype,'count',count);        
      else
        fprintf(fid,'Cluster count file does not exist: %s\n',cluster_file);
        return;
      end
    elseif strcmp(params.chantype,'grad') && strcmp(params.cluster_count,'eeg_pospeak')
      cluster_file = sprintf('matfiles/SO_clustered_eeg_pospeak_detections.mat');
      fprintf(fid,'Clustering MEG data using EEG pospeak clusters as reference.\n');
      if exist(cluster_file,'file')
        fprintf(fid,'EEG cluster file: %s\n',cluster_file);
        % load EEG negative half-wave agg count
        load('matfiles/SO_clustered_eeg_pospeak_detections.mat','count','CountTime');
        count = count(CountTime<detections(1).tstop);
        % use EEG peaks to cluster MEG detections
        [detections,count,CountTime,ClusterTime,CountThresh] = SO_cluster_detections(...
          detections,'method','histogram','threshold',params.cluster_thresh,...
          'StepSize',params.cluster_StepSize,'IntegrationWindow',params.cluster_IntegrationWindow,...
          'MinSeparation',params.cluster_MinSeparation,'ClusterWindow',params.cluster_ClusterWindow,...
          'peaktype',params.peaktype,'count',count);        
      else
        fprintf(fid,'Cluster count file does not exist: %s\n',cluster_file);
        return;
      end
    else
      if strcmp(params.peaktype,'pospeak')
        [detections.negpeak] = deal([]);
      elseif strcmp(params.peaktype,'negpeak')
        [detections.pospeak] = deal([]);
      end
      [detections,count,CountTime,ClusterTime,CountThresh] = SO_cluster_detections(detections,'method','histogram',...
        'threshold',params.cluster_thresh,'StepSize',params.cluster_StepSize,...
        'IntegrationWindow',params.cluster_IntegrationWindow,'MinSeparation',params.cluster_MinSeparation,...
        'ClusterWindow',params.cluster_ClusterWindow,'peaktype',params.peaktype);%,'count',count);    
    end
    % 1. this function adds the following fields to the detections structure:
    % cluster_number, cluster_time_index, pospeak_cluster_number,
    % pospeak_cluster_time_index,         negpeak_cluster_number, and
    % negpeak_cluster_time_index.
    % 2. count = windowed aggregate detection count at each step.   
    clustered_detections = detections;
    %% Apply post-cluster criteria
    % remove detections in channels that don't have at least 1 pos and neg
    % detection
%     tmp1 = cellfun(@length,{detections.pospeak});
%     tmp2 = cellfun(@length,{detections.negpeak});
%     tmp  = find(tmp1==0 | tmp2==0);
%     if ~isempty(tmp)
%       fprintf(fid,'Ignoring %g channels with less than 1 pos and 1 neg half-wave peak detection.\n',length(tmp));
%       [detections(tmp).pospeak]                     = deal([]);
%       [detections(tmp).negpeak]                     = deal([]);
%       [detections(tmp).pospeak_cluster_number]      = deal([]);
%       [detections(tmp).negpeak_cluster_number]      = deal([]);
%       [detections(tmp).pospeak_cluster_time_index]  = deal([]);
%       [detections(tmp).negpeak_cluster_time_index]  = deal([]);
%       [detections(tmp).cluster_number]              = deal([]);
%       [detections(tmp).cluster_time_index]          = deal([]);
%     end
    % remove clusters with less than some number of involved channels
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
    Nthresh   = params.MinChansPerCluster;
    ckeep     = find(Ninvolved > Nthresh);
    % Eliminate clusters with (# involved channels) < threshold
    fprintf(fid,'\tRemoving clusters involving fewer than %g channels.\n',params.MinChansPerCluster);
    fprintf(fid,'\tRemoved %g clusters.\n',length(clusterID)-length(ckeep));
    Ninvolved = Ninvolved(ckeep);
    clusterID = clusterID(ckeep);
    clusterIX = clusterIX(ckeep);
    for k = 1:length(detections)
      pkeep   = find(ismember(detections(k).pospeak_cluster_number,clusterID));
      nkeep   = find(ismember(detections(k).negpeak_cluster_number,clusterID));
      bkeep   = find(ismember(detections(k).cluster_number,clusterID));
      detections(k).pospeak                     = detections(k).pospeak(pkeep);
      detections(k).pospeak_cluster_number      = detections(k).pospeak_cluster_number(pkeep);
      detections(k).pospeak_cluster_time_index  = detections(k).pospeak_cluster_time_index(pkeep);
      detections(k).negpeak                     = detections(k).negpeak(nkeep);
      detections(k).negpeak_cluster_number      = detections(k).negpeak_cluster_number(nkeep);
      detections(k).negpeak_cluster_time_index  = detections(k).negpeak_cluster_time_index(nkeep);
      detections(k).cluster_number              = detections(k).cluster_number(bkeep);
      detections(k).cluster_time_index          = detections(k).cluster_time_index(bkeep);
    end
%     tmp1 = cellfun(@length,{detections.pospeak});
%     tmp2 = cellfun(@length,{detections.negpeak});
%     tmp  = find(tmp1==0 | tmp2==0);
%     if ~isempty(tmp)
%       [detections(tmp).pospeak]                     = deal([]);
%       [detections(tmp).negpeak]                     = deal([]);
%       [detections(tmp).pospeak_cluster_number]      = deal([]);
%       [detections(tmp).negpeak_cluster_number]      = deal([]);
%       [detections(tmp).pospeak_cluster_time_index]  = deal([]);
%       [detections(tmp).negpeak_cluster_time_index]  = deal([]);
%       [detections(tmp).cluster_number]              = deal([]);
%       [detections(tmp).cluster_time_index]          = deal([]);
%     end
    proc_clustered_detections = detections;  
    t       = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
    events  = [];
    for k   = 1:length(detections)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
      events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
    end
%     toc
    fprintf(fid,'Detected %g clusters during t = %g-%gsec.\n',length(clusterID),t(1),t(end));
    fprintf(fid,'Saving clustered detections: %s\n',file);
    save(file,'events','detections','count','params','CountTime','ClusterTime','CountThresh','clusterID','clusterIX','Ninvolved');
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  
  file = sprintf('matfiles/SO_clustered_%s_%s_detections_sws.mat',params.chantype,params.peaktype);
  if strcmp(params.chantype,'grad') && ~isempty(params.cluster_count)
    file = sprintf('matfiles/SO_clustered_%s_%s_detections_sws_ref-%s.mat',params.chantype,params.peaktype,params.cluster_count);
  end  
  if exist(file,'file')
    fprintf(fid,'Loading SO %s %s clusters detected during slow-wave sleep: %s\n',params.chantype,params.peaktype,file);
    load(file,'detections','BegPoints','EndPoints','IntervalT0','IntervalTf','t0','tf','count');
    tmp  = arrayfun(@(x)(x.cluster_number),detections,'UniformOutput',false);
    tmp  = unique([tmp{:}]);
    fprintf(fid,'Analyzing %g clusters detected during slow wave sleep:\n',length(tmp));
    clear tmp
    for k = 1:length(IntervalT0)
      fprintf(fid,'t = %g - %g sec (%3.3g min)\n',IntervalT0(k),IntervalTf(k),(IntervalTf(k)-IntervalT0(k))/60);
    end
  elseif params.SO_cluster_detections_flag
    %% Automatic selection of SWS
    fprintf(fid,'Automatically detecting periods of slow-wave sleep based on SO detection clusters in %s.\n',params.chantype);
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
    fprintf(fid,'\tCount SO clusters in %gsec epochs.\n',params.AASM_EpochLength);
    nn    = round((max(tc)-min(tc))/params.AASM_EpochLength);
    [N,X] = hist(tc,nn);
    tmpx  = zeros(1,length(t));
    tmpx(clusterIX) = 1;
    if params.ShowPlots
      figure('Name','Detection Count');
      subplot(4,1,1),plot(t,tmpx,'.-'); axis tight; title('Automatic detection of slow oscillations'); ylabel('SO detected'); xlabel('time (sec)');
      subplot(4,1,2),try hist(tc,nn); end; axis tight; xlim=get(gca,'xlim'); title('Aggregate detection count histogram');
      xlabel(sprintf('time (%gsec bins)',params.AASM_EpochLength)); ylabel('count');
      subplot(4,1,3),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
      title('SO detection count with threshold for SWS based on AASM guidelines');
      xlabel(sprintf('time (%gsec bins)',params.AASM_EpochLength)); ylabel('count');
    end
    % Select time intervals of interest from detection histogram
    if strcmp(params.IntervalSelection_CountThreshold,'mean')
      thresh = mean(N);
    elseif strcmp(params.IntervalSelection_CountThreshold,'meanstd')
      thresh = mean(N) + std(N);
    elseif strcmp(params.IntervalSelection_CountThreshold,'median')
      thresh = median(N);
    elseif isnumeric(params.IntervalSelection_CountThreshold)
      thresh = params.IntervalSelection_CountThreshold;
    else
      thresh = mean(N);
    end
    fprintf(fid,'\tLabel an epoch as slow-wave sleep if it has more than %g clusters.\n',thresh);
      % threshold should be doubled for EEG
    tmp     = N > thresh;                         % 1 if count > thresh; else 0
    if ~any(tmp)
      fprintf(fid,'Error: No intervals satisfy the acceptance criteria.\n');
      error('No intervals satisfy the acceptance criteria.');
    end
    L = find(tmp(2:end)==1 & tmp(1:end-1)==0);    % index to left edge
    R = find(tmp(1:end-1)==1 & tmp(2:end)==0);    % index to right edge
    if tmp(1)==1,   L = [1 L];         end
    if tmp(end)==1, R = [R length(X)]; end
    if L(1)-R(1)==1                               % correct for right edge at beginning
      L = L(1:end-1);
      R = R(2:end);
    end
    if length(R)>1 && (R(1)    <L(1)),   R = R(2:end);   end
    if length(L)>1 && (L(end)  >R(end)), L = L(1:end-1); end
    if length(R)>1 && (R(end-1)>L(end)), R = R(1:end-1); end
    % combine threshold crossings if right edge - left edge < some min time
    fprintf(fid,'\tPatch: relabel a period below threshold for <= %g epoch(s) if it is bounded by slow-wave epochs.\n',floor(params.IntervalSelection_CombineIfLessThan/params.AASM_EpochLength));
    if length(L) > 2
      I       = find(X(L(2:end)) - X(R(1:end-1)) < params.IntervalSelection_CombineIfLessThan);
      L(I+1)  = [];
      R(I)    = [];
    end
    % require that count exceed CountThreshold for > MinCrossingTime
    fprintf(fid,'\tDefine a multi-epoch interval as slow-wave sleep if it contains at least %g slow-wave epochs.\n',ceil(params.IntervalSelection_MinCrossingTime/params.AASM_EpochLength));
    D   = X(R) - X(L);
    I   = find(D > params.IntervalSelection_MinCrossingTime);
    L   = L(I);
    R   = R(I);
    % get start and stop times for each interval
    t0  = X(L);
    tf  = X(R);
    t0f = [t0 tf];
    % interval limits
    nInterval = length(t0);
    nchan     = length(detections);
    newdetections  = detections;
    [newdetections(1:nchan).pospeak]                     = deal([]);
    [newdetections(1:nchan).negpeak]                     = deal([]);
    [newdetections(1:nchan).pospeak_cluster_time_index]  = deal([]);
    [newdetections(1:nchan).pospeak_cluster_number]      = deal([]);
    [newdetections(1:nchan).negpeak_cluster_time_index]  = deal([]);
    [newdetections(1:nchan).negpeak_cluster_number]      = deal([]);
    [newdetections(1:nchan).cluster_time_index]          = deal([]);
    [newdetections(1:nchan).cluster_number]              = deal([]);
    % loop over limits and constrain detections
    for k = 1:nInterval
      % PEAKS
      % convert time limits to indices
      L   = nearest(t,t0(k));
      R   = nearest(t,tf(k));
      fld = 'pospeak';   tmp = {detections.(fld)};  pkeep = cellfun(@(x)(x>=L&x<=R),tmp,'uniformoutput',false);
      tmp = cellfun(@(x)x(x>=L&x<=R),tmp,'uniformoutput',false); 
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
      fld = 'negpeak';   tmp = {detections.(fld)};  nkeep = cellfun(@(x)(x>=L&x<=R),tmp,'uniformoutput',false);
      tmp = cellfun(@(x)x(x>=L&x<=R),tmp,'uniformoutput',false); 
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:}); 
      fld = 'pospeak_cluster_time_index';           keep  = pkeep;% cellfun(@(x)(x>=L&x<=R),{detections.(fld)},'uniformoutput',false);
      tmp = cellfun(@(x,y)x(y),{detections.(fld)},  keep,'uniformoutput',false); 
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
      fld = 'pospeak_cluster_number'; 
      tmp = cellfun(@(x,y)(x(y)),{detections.(fld)},keep,'UniformOutput',false); 
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
      fld = 'negpeak_cluster_time_index';           keep  = nkeep;%keep = cellfun(@(x)(x>=L&x<=R),{detections.(fld)},'uniformoutput',false);
      tmp = cellfun(@(x,y)x(y),{detections.(fld)},  keep,'uniformoutput',false); 
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
      fld = 'negpeak_cluster_number'; 
      tmp = cellfun(@(x,y)(x(y)),{detections.(fld)},keep,'UniformOutput',false);
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
      fld = 'cluster_time_index';      keep = cellfun(@(x)(x>=L&x<=R),{detections.(fld)},'uniformoutput',false);
      tmp = cellfun(@(x,y)x(y),{detections.(fld)},  keep,'uniformoutput',false); 
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
      fld = 'cluster_number'; 
      tmp = cellfun(@(x,y)(x(y)),{detections.(fld)},keep,'UniformOutput',false);
      tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});  
      % DATA
      if k == 1
        seldat = ts_data_selection(data,'toilim',[t0(k) tf(k)]);
        seldat.epochs.IntervalEndPoints(k) = size(seldat.epochs.data,2);
      else
        tmpdat = ts_data_selection(data,'toilim',[t0(k) tf(k)]);
        seldat.epochs.time = [seldat.epochs.time tmpdat.epochs.time];
        seldat.epochs.data = cat(2,seldat.epochs.data,tmpdat.epochs.data);
        seldat.epochs.IntervalEndPoints(k) = size(seldat.epochs.data,2);
        clear tmpdat
      end
    end  
    detections = newdetections; 
    if strcmp(params.chantype,'grad')
      tmp1 = cellfun(@length,{detections.pospeak});
      tmp2 = cellfun(@length,{detections.negpeak});
      tmp  = find(tmp1==0 | tmp2==0);
      if ~isempty(tmp)
        [detections(tmp).pospeak]                     = deal([]);
        [detections(tmp).negpeak]                     = deal([]);
        [detections(tmp).pospeak_cluster_number]      = deal([]);
        [detections(tmp).negpeak_cluster_number]      = deal([]);
        [detections(tmp).pospeak_cluster_time_index]  = deal([]);
        [detections(tmp).negpeak_cluster_time_index]  = deal([]);
        [detections(tmp).cluster_number]              = deal([]);
        [detections(tmp).cluster_time_index]          = deal([]);
      end    
    end
    proc_clustered_detections_SWS = detections;
%     toc
    % update reference times and show clusters which will be processed
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
    NN = zeros(1,length(N));
    for k = 1:length(t0)
      ix      = X>=t0(k) & X<=tf(k);
      NN(ix)  = N(ix);
    end
    if params.ShowPlots
      subplot(4,1,4),plot(X,NN,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r'); set(gca,'xlim',xlim);
      title('AASM guideline-based automatic N3 detection');
      xlabel(sprintf('time (%gsec bins)',params.AASM_EpochLength)); ylabel('count');
    end
    params.t0  = t0; % start time of each interval of interest
    params.tf  = tf; % end time of each interval of interest
    % total_toilim = [params.t0(1) params.tf(end)];
    clear newdetections tmp L R k fld keep nInterval t nchan    
    BegPoints       = [1 seldat.epochs.IntervalEndPoints(1:end-1)+1]; % indices
    EndPoints       = seldat.epochs.IntervalEndPoints;                % indices
    IntervalT0      = seldat.epochs.time(BegPoints);                  % sec
    IntervalTf      = seldat.epochs.time(EndPoints);                  % sec
    fprintf(fid,'Defined %g periods of slow-wave sleep with %g SO clusters in %s.\n',length(IntervalT0),nclusters,params.chantype);
    for k = 1:length(IntervalT0)
      fprintf(fid,'\tt = %g - %g sec\n',IntervalT0(k),IntervalTf(k));
    end
    % The events structure can be used by visualizer to mark detections  
    t       = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
    events  = [];
    for k   = 1:length(detections)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
      events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
    end
    fprintf(fid,'Saving clusters detected in slow-wave sleep: %s\n',file);
    save(file,'events','detections','count','params','BegPoints','EndPoints','IntervalT0','IntervalTf','t0','tf'); % ,'detections_alltime'
    clear events t k seldat
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % visualizer(seldat); % load SO_clustered_detections.mat    
%     if params.ShowPlots
%       % Compare results
%       vars  = {'clustered_detections','proc_clustered_detections','proc_clustered_detections_SWS'};
%       ns    = length(vars);
%       figure('Name','Detection Count'); 
%       for s = 1:ns
%         det = eval(vars{s});
%         t         = det(1).tstart:1/det(1).sfreq:det(1).tstop;
%         typestr   = ''; % {'','pospeak_','negpeak_'} where ''=>both
%         cnumfield = sprintf('%scluster_number',typestr);
%         cindfield = sprintf('%scluster_time_index',typestr);
%         clusterID = arrayfun(@(x)(x.(cnumfield)),det,'UniformOutput',false);
%         tmpID     = [clusterID{:}];
%         clusterID = unique(tmpID); % sorted cluster index, 1:nclusters
%         clusterIX = arrayfun(@(x)(x.(cindfield)),det,'UniformOutput',false);
%         clusterIX = [clusterIX{:}];
%         clusterIX = cellfun(@(x)unique(clusterIX(tmpID==x)),num2cell(unique(tmpID)));
%         nclusters = length(clusterID);
%         cnum      = {det.(cnumfield)};
%         cnum      = [cnum{:}];
%         Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));
%         % create cell array listing channels involved in each cluster
%         cnum          = {det.(cnumfield)};
%         tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
%         tmp           = [tmp{:}]';
%         tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
%         InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);
%         % center times for each cluster
%         tc            = t(clusterIX);
%         t0f           = [IntervalT0 IntervalTf];
%         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % % %  POTENTIAL FIGURE FOR PAPER
%         % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%         nn    = round((max(tc)-min(tc))/params.AASM_EpochLength);
%         [N,X] = hist(tc,nn);
%         tmpx  = zeros(1,length(t));
%         tmpx(clusterIX) = 1;
%         subplot(4,ns,s+0*ns),plot(t,tmpx,'.-'); axis tight; title('aggregate cluster count')
%         subplot(4,ns,s+1*ns),try hist(tc,nn); end; axis tight
%         subplot(4,ns,s+2*ns),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
%         xlabel('time (sec)'); ylabel(sprintf('count (%gsec bins)',params.AASM_EpochLength));
%         for k=1:length(t0f), vline(t0f(k),'k'); end
%         subplot(4,ns,s+3*ns),plot(tc,'.-'); ylabel('cluster time (sec)'); xlabel('cluster number'); axis tight;
%         clear nn N X tmpx
%       end
%     end
%     toc
  end
  clear file

  if ~exist('detections','var')
    % end function because all subsequent processing requires detections.
    return
  end    

  %% Single peak waves
  file = sprintf('matfiles/SO_clustered_%s_%s_detections_singlepeaks.mat',params.chantype,params.peaktype);
  if strcmp(params.chantype,'grad') && ~isempty(params.cluster_count)
    file = sprintf('matfiles/SO_clustered_%s_%s_detections_singlepeaks_ref-%s.mat',params.chantype,params.peaktype,params.cluster_count);
  end    
  if exist(file,'file') && params.monotonic
    fprintf(fid,'Loading single peak slow-wave sleep SO %s %s detections: %s\n',params.chantype,params.peaktype,file);
    load(file,'detections');
    tmp1 = cellfun(@length,{detections.pospeak});
    tmp2 = cellfun(@length,{detections.negpeak});
    tmp  = mean([tmp1;tmp2],1);
    fprintf(fid,'Set contains %g detections across all channels.\n',round(sum(tmp)));
    clear tmp tmp1 tmp2    
  elseif params.monotonic
    args              = mmil_parms2args(params);
    fprintf(fid,'Detecting single peak slow oscillations in slow-wave sleep time series ranging t = %g - %g sec.\n',data.epochs.time(1),data.epochs.time(end));   
%     detections        = SO_detection(data,args{:});
    % Instead of rerunning the detection algorithm, just remove
    % nonmonotonic phase runs from the previous detections
    load(sprintf('matfiles/SO_%s_detections.mat',params.chantype),'detections');
    tmpdet            = detections;
    [tmpdet.pospeak]  = deal([]);
    [tmpdet.negpeak]  = deal([]);
    % loop over channels
    fprintf(fid,'Removing multipeak waves (nonmonotonic phase runs) from prior detections.\n');
    if isempty(params.blcwindow),blcwindow = [-inf inf]; else blcwindow = params.blcwindow; end
    tmptic= tic;
    for k = 1:length(detections)
      % prep data
      % select data
      Fs  = data.sfreq;
      dat = ts_data_selection(data,'chanlabel',data.sensor_info(k).label,'verbose',0);
      t   = dat.epochs.time;
      x   = dat.epochs.data;
      % bandpass filter
      x   = ts_freq_filt(x,Fs,params.bpfreq,[0 0],'bandpass');
      % baseline correction
      if params.blc
        tix = nearest(t,blcwindow(1)):nearest(t,blcwindow(2));
        x   = x - mean(x(tix));
      end
      % decimate
      if params.decimate
        if isa(x,'single')
          x = double(x);
          x = decimate(x,params.decimate_factor);
          x = single(x);
        else
          x = decimate(x,params.decimate_factor);
        end
        Fs = Fs / params.decimate_factor;
        t  = downsample(t,params.decimate_factor);
        if length(x) ~= length(t), error('Something went wrong with the decimation'); end
      end
      % smooth (moving average)
      if params.smooth
        x = smooth(x,ceil(params.smooth_window*Fs),params.smooth_method);
      end      
      % positive half-waves
      pks = detections(k).pospeak;
      % calculate instantaneous SO phase
      phi = angle(hilbert(x));     
      % simplify by working with the absolute value (just pi-pi/12)
      absphi = abs(phi);
      % find pi-pi/12 & -pi+pi/12 crossings
      ind    = crossing(absphi,[],pi-pi/12);
      % find pi-pi/12 crossing before each peak
      L      = cellfun(@(x)(ind(find(ind<x,1,'last'))),num2cell(pks),'UniformOutput',false);
      % find pi-pi/12 crossing after each peak
      R      = cellfun(@(x)(ind(find(ind>x,1,'first'))),num2cell(pks),'UniformOutput',false);
      sel    = ~(cellfun(@isempty,L) | cellfun(@isempty,R));
      pks    = pks(sel);
      L      = [L{sel}];
      R      = [R{sel}];
      % is monotonic ?
      sel    = cellfun(@(x,y)(ismonotonic(phi(x:y))),num2cell(L),num2cell(R));
      pks    = pks(sel);
      tmpdet(k).pospeak = pks;
      % flip x and repeat for negative half-waves
      x      = -x;
      pks = detections(k).negpeak;
      phi = angle(hilbert(x));     
      absphi = abs(phi);
      ind    = crossing(absphi,[],pi-pi/12);
      L      = cellfun(@(x)(ind(find(ind<x,1,'last'))),num2cell(pks),'UniformOutput',false);
      R      = cellfun(@(x)(ind(find(ind>x,1,'first'))),num2cell(pks),'UniformOutput',false);
      sel    = ~(cellfun(@isempty,L) | cellfun(@isempty,R));
      pks    = pks(sel);
      L      = [L{sel}];
      R      = [R{sel}];
      sel    = cellfun(@(x,y)(ismonotonic(phi(x:y))),num2cell(L),num2cell(R));
      pks    = pks(sel);
      tmpdet(k).negpeak = pks;
      telapsed = toc(tmptic); fprintf('Time elapsed after chan %g of %g: %g min\n',k,length(detections),telapsed/60);
    end    
    detections = tmpdet;
    clear pks L R sel absphi phi x tmpdet tmptic
    % only keep detections in intervals of interest
    t                 = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
    tmpdet            = detections;
    [tmpdet.pospeak]  = deal([]);
    [tmpdet.negpeak]  = deal([]);
    [tmppos{1:length(tmpdet)}] = deal([]);
    [tmpneg{1:length(tmpdet)}] = deal([]);
    fprintf(fid,'Selecting single peak slow oscillations during slow-wave sleep:\n');
    for k = 1:length(IntervalT0)
      a   = IntervalT0(k);
      b   = IntervalTf(k);
      fprintf(fid,'\tt = %g - %g sec\n',a,b);
      pks = {detections.pospeak};
      pks = cellfun(@(x)(x(t(x)>=a & t(x)<=b)),pks,'uniformoutput',false);
      tmppos = cellfun(@(x,y)([x y]),tmppos,pks,'uniformoutput',false);
      pks = {detections.negpeak};
      pks = cellfun(@(x)(x(t(x)>=a & t(x)<=b)),pks,'uniformoutput',false);
      tmpneg = cellfun(@(x,y)([x y]),tmpneg,pks,'uniformoutput',false);
      clear pks a b
    end
    [tmpdet.pospeak] = deal(tmppos{:});
    [tmpdet.negpeak] = deal(tmpneg{:});
    detections       = tmpdet;
    clear tmpdet tmppos tmpneg args
    if params.peakpairs_flag
      % This will eliminate all positive detections without an adjacent negative
      % peak and vice versa.  detections must occur within tau of each other.
      detections  = select_peakpairs(detections,params.peakpairs_tau);
    end    
    tmp1 = cellfun(@length,{detections.pospeak});
    tmp2 = cellfun(@length,{detections.negpeak});
    tmp  = find(tmp1==0 | tmp2==0);
    if ~isempty(tmp)
      [detections(tmp).pospeak]                     = deal([]);
      [detections(tmp).negpeak]                     = deal([]);
      [detections(tmp).pospeak_cluster_number]      = deal([]);
      [detections(tmp).negpeak_cluster_number]      = deal([]);
      [detections(tmp).pospeak_cluster_time_index]  = deal([]);
      [detections(tmp).negpeak_cluster_time_index]  = deal([]);
      [detections(tmp).cluster_number]              = deal([]);
      [detections(tmp).cluster_time_index]          = deal([]);
    end    
    events  = [];
      % The events structure can be used by visualizer to mark detections
    for k   = 1:length(detections)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
      events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
    end
    fprintf(fid,'Saving SO single peak detections: %s\n',file);
    save(file,'events','detections','params');
    clear events args t k
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % visualizer(data); % load SO_detections.mat
  end
  clear file
  
  file = sprintf('matfiles/SO_clustered_%s_%s_detections_sws_singlepeaks.mat',params.chantype,params.peaktype);
  if strcmp(params.chantype,'grad') && ~isempty(params.cluster_count)
    file = sprintf('matfiles/SO_clustered_%s_%s_detections_sws_singlepeaks_ref-%s.mat',params.chantype,params.peaktype,params.cluster_count);
  end      
  if exist(file,'file') && params.monotonic
    fprintf(fid,'Loading SO %s %s single peak detection clusters: %s\n',params.chantype,params.peaktype,file);
    load(file,'detections','BegPoints','EndPoints','IntervalT0','IntervalTf','count');
    try
      load(file,'t0','tf');
    catch
      t0 = IntervalT0;
      tf = IntervalTf;
    end
    tmp  = arrayfun(@(x)(x.cluster_number),detections,'UniformOutput',false);
    tmp  = unique([tmp{:}]);
    fprintf(fid,'Analyzing %g single peak clusters detected during slow wave sleep:\n',length(tmp));
    clear tmp
    for k = 1:length(IntervalT0)
      fprintf(fid,'t = %g - %g sec\n',IntervalT0(k),IntervalTf(k));
    end    
  elseif params.monotonic
    fprintf(fid,'Finding single peak detection clusters using prior aggregate count.\n');
    if strcmp(params.chantype,'grad') && ~isempty(params.cluster_count)
      file = sprintf('matfiles/SO_clustered_%s_%s_detections_singlepeaks_ref-%s.mat',params.chantype,params.peaktype,params.cluster_count);
      load(sprintf('matfiles/SO_clustered_%s_%s_detections_singlepeaks_ref-%s.mat',params.chantype,params.peaktype,params.cluster_count),'count');
    else
      load(sprintf('matfiles/SO_clustered_%s_%s_detections_sws.mat',params.chantype,params.peaktype),'count');
    end        
    % remove count steps outside of the bounds on the intervals of interest
    t     = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
    Fs    = data.sfreq;
    L     = floor(params.cluster_IntegrationWindow*Fs/2);
    n     = floor(params.cluster_StepSize*Fs);
    nsmp  = length(t);
    c     = L+1:n:nsmp-L;
    tc    = t(c);
    count = count(tc>=detections(1).tstart & tc<=detections(1).tstop);
    clear L n nsmp c tc
    [detections,count,CountTime,ClusterTime,CountThresh] = SO_cluster_detections(detections,'method',params.cluster_method,...
      'threshold',params.cluster_thresh,'StepSize',params.cluster_StepSize,...
      'IntegrationWindow',params.cluster_IntegrationWindow,'MinSeparation',params.cluster_MinSeparation,...
      'ClusterWindow',params.cluster_ClusterWindow,'count',count);
      % 1. this function adds the following fields to the detections structure:
      % cluster_number, cluster_time_index, pospeak_cluster_number,
      % pospeak_cluster_time_index,         negpeak_cluster_number, and
      % negpeak_cluster_time_index.
      % 2. count = windowed aggregate detection count at each step.   
    % remove detections in channels that don't have at least 1 pos and neg
    % detection
%     detections          = CheckDetections(detections);
    % remove clusters with less than some number of involved channels
    fprintf(fid,'Eliminating clusters with fewer than %g involved channels.\n',params.MinChansPerCluster);
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
    Nthresh   = params.MinChansPerCluster;
    ckeep     = find(Ninvolved > Nthresh);
    % Eliminate clusters with (# involved channels) < threshold
    Ninvolved = Ninvolved(ckeep);
    clusterID = clusterID(ckeep);
    clusterIX = clusterIX(ckeep);
    for k = 1:length(detections)
      pkeep   = find(ismember(detections(k).pospeak_cluster_number,clusterID));
      nkeep   = find(ismember(detections(k).negpeak_cluster_number,clusterID));
      bkeep   = find(ismember(detections(k).cluster_number,clusterID));
      detections(k).pospeak                     = detections(k).pospeak(pkeep);
      detections(k).pospeak_cluster_number      = detections(k).pospeak_cluster_number(pkeep);
      detections(k).pospeak_cluster_time_index  = detections(k).pospeak_cluster_time_index(pkeep);
      detections(k).negpeak                     = detections(k).negpeak(nkeep);
      detections(k).negpeak_cluster_number      = detections(k).negpeak_cluster_number(nkeep);
      detections(k).negpeak_cluster_time_index  = detections(k).negpeak_cluster_time_index(nkeep);
      detections(k).cluster_number              = detections(k).cluster_number(bkeep);
      detections(k).cluster_time_index          = detections(k).cluster_time_index(bkeep);
    end
%     detections          = CheckDetections(detections);
    tmp1 = cellfun(@length,{detections.pospeak});
    tmp2 = cellfun(@length,{detections.negpeak});
    tmp  = find(tmp1==0 | tmp2==0);
    if ~isempty(tmp)
      [detections(tmp).pospeak]                     = deal([]);
      [detections(tmp).negpeak]                     = deal([]);
      [detections(tmp).pospeak_cluster_number]      = deal([]);
      [detections(tmp).negpeak_cluster_number]      = deal([]);
      [detections(tmp).pospeak_cluster_time_index]  = deal([]);
      [detections(tmp).negpeak_cluster_time_index]  = deal([]);
      [detections(tmp).cluster_number]              = deal([]);
      [detections(tmp).cluster_time_index]          = deal([]);
    end
    detections_alltime  = detections;  
    t       = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
    events  = [];
    for k   = 1:length(detections)
      events(k).label = data.sensor_info(k).label;
      events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
      events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
    end
    fprintf(fid,'Saving clustered single peak detections: %s\n',file);
    save(file,'events','detections','count','params','BegPoints','EndPoints','IntervalT0','IntervalTf','detections_alltime','t0','tf');
    clear events t k count seldat ckeep pkeep nkeep bkeep
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  clear file
  
  if params.monotonic
    strID = [strID '_singlepeaks'];
  end

  % check that slow-wave sleep detection has been performed
  if ~exist('IntervalT0','var'), return; end
  
  %% IMPORTANT NOTE:
  % The detection structure obtained at this point will be used for each of
  % the four correction methods below.
  detections_common = detections;
    
  % Select data containing the intervals of interest
  fprintf(fid,'Selecting data subset containing slow-wave sleep: %g-%gsec\n',IntervalT0(1),IntervalTf(end));
  toilim = [IntervalT0(1) IntervalTf(end)];
%   data   = ts_data_selection(data,'toilim',toilim);

  % reprocess cluster detections after selecting time intervals
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
    % for each (x), clusterIX(tmpID==x) should be index to the time of the
    % histpeak=const for all channels involved in cluster (x).
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
  if params.ShowPlots
    % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %  POTENTIAL FIGURE FOR PAPER
    % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    figure('Name','Detection Count'); 
    nn    = round((max(tc)-min(tc))/params.AASM_EpochLength);
    [N,X] = hist(tc,nn);
    tmpx  = zeros(1,length(t));
    tmpx(clusterIX) = 1;
    subplot(4,1,1),plot(t,tmpx,'.-'); axis tight; title('Aggregate slow oscillation detection count')
    subplot(4,1,2),try hist(tc,nn); end; axis tight
    subplot(4,1,3),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
    xlabel('time (sec)'); ylabel(sprintf('count (%gsec bins)',params.AASM_EpochLength));
    for k=1:length(IntervalT0), vline(IntervalT0(k),'r'); end
    for k=1:length(IntervalTf), vline(IntervalTf(k),'k'); end
%     for k=1:length(t0f), vline(t0f(k),'k'); end
    subplot(4,1,4),plot(tc,'.-'); ylabel('cluster time (sec)'); xlabel('cluster number'); axis tight;
    clear nn N X tmpx t0f
    % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  end      
    
  %%
    %% SO FLIPPING: method 1 (use the detection clusters already found w/o flipping data.)
  file = sprintf('matfiles/SO_clusters_%s_%s_noflip%s.mat',params.chantype,params.peaktype,strID);  %['matfiles/SO_clusters_noflip' strID '.mat'];
  if strcmp(params.chantype,'grad') && ~isempty(params.cluster_count)
    file = sprintf('matfiles/SO_clusters_%s_%s_noflip%s_ref-%s.mat',params.chantype,params.peaktype,strID,params.cluster_count);
  end
  if exist(file,'file')
    fprintf(fid,'Loading SO %s %s clusters with corrcoef (no flipping): %s\n',params.chantype,params.peaktype,file);
    load(file,'clusters_Near_pks','clusters_Nmax_pks','clusters_Near_all','clusters_Nmax_all');
  elseif params.SO_cluster_corrcoef1_flag
    fprintf(fid,'Creating %s %s clusters struct and calculating corrcoefs without flipping.\n',params.chantype,params.peaktype);
    % Create clusters structure and calculate correlation coefficients between
    % a cluster-specific reference (params.cluster_RefChanType) and the other
    % channels involved in the cluster.
    fprintf(fid,'Filter data (%g-%gHz) and remove mean.\n',params.bpfreq);
    fprintf(fid,'Analyzing distance vs delay correlations (%g-%gsec).\n',toilim);
    args        = mmil_parms2args(params);
    % prepare data
    procdat     = ts_data_selection(data,'toilim',toilim);
    procdat     = ts_preproc(procdat,'bpfilter','yes','bpfreq',params.bpfreq,'bandpass_detrend_flag',0,'blc','yes','blcwindow',[-inf inf]);
    % Split pos and neg within clusters when calculating R
      % EarliestDetection
    clusters_Near_pks = SO_cluster_corrcoef(procdat,detections,args{:},'cluster_RefChanType','EarliestDetection','preproc',0);
      % MaxR
    clusters_Nmax_pks = SO_cluster_corrcoef(procdat,detections,args{:},'cluster_RefChanType','MaxR','preproc',0);
    % Use all detections within a cluster when calculating R
    % prepare flipdat & args
    tmpdet  = detections;
    tmppos  = {detections.pospeak};
    tmpneg  = {detections.negpeak};
    tmp     = cellfun(@(x,y)[x y],tmppos,tmpneg,'uniformoutput',false);
    [tmpdet.negpeak]  = deal([]);
    [tmpdet.pospeak]  = deal(tmp{:});
      % EarliestDetection
    clusters_Near_all = SO_cluster_corrcoef(procdat,tmpdet,args{:},'cluster_RefChanType','EarliestDetection','preproc',0);
      % MaxR
    clusters_Nmax_all = SO_cluster_corrcoef(procdat,tmpdet,args{:},'cluster_RefChanType','MaxR','preproc',0);
    clear tmp tmppos tmpneg tmpdet args procdat
    fprintf(fid,'Saving clusters struct (no flipping): %s\n',file);
    save(file,'clusters_Near_pks','clusters_Nmax_pks','clusters_Near_all','clusters_Nmax_all','detections','params');
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  clear file

  
  %% fin
  file = ['matfiles/SO_params' strID '.mat'];
  fprintf(fid,'Saving parameter file: %s\n',file);
  save(file,'params');
  telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  fclose(fid);
  cd(cwd);
  if nargout > 0 && ~exist('varargout','var')
    varargout{1} = detections;
  end
% catch
%   fprintf(fid,'%s',rethrow(lasterr));
%   telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
%   fclose(fid);
%   cd(cwd);
% end
%    