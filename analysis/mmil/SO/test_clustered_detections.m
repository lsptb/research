clear all
cd /space/emc2/1/halgdev/projects/sleep/MEG/SO
load s1/matfiles/proc_epoch_data_ICA.mat; data=epoch_data; clear epoch_data;
load s1/matfiles/SO_detections.mat
tic
all_detections = detections;

params.cluster_StepSize           = .010; % sec
params.cluster_MinSeparation      = .05; % sec, MinLength is also set to this
params.cluster_IntegrationWindow  = .05;  % sec
params.cluster_ClusterWindow      = .2;   % sec
params.cluster_thresh             = 5;%'meanstd3';

[detections,count,CountTime,ClusterTime,CountThresh] = SO_cluster_detections(detections,'method','histogram',...
  'threshold',params.cluster_thresh,'StepSize',params.cluster_StepSize,...
  'IntegrationWindow',params.cluster_IntegrationWindow,'MinSeparation',params.cluster_MinSeparation,...
  'ClusterWindow',params.cluster_ClusterWindow);%,'count',count);
toc
clustered_detections = detections;

%%
params.MinChansPerCluster                   = 40;%round(.20*data.num_sensors);%10; % # of involved chans
params.AASM_EpochLength                     = 30; % sec
params.IntervalSelection_CountThreshold     = 10; % # of clusters
params.IntervalSelection_CombineIfLessThan  = 45; % sec (1 epoch)
params.IntervalSelection_MinCrossingTime    = 60; % sec (2 epochs)

% remove detections in channels that don't have at least 1 pos and neg
% detection
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
proc_clustered_detections = detections;  
t       = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
events  = [];
for k   = 1:length(detections)
  events(k).label = data.sensor_info(k).label;
  events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
  events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
end

toc
%% Automatic selection of SWS

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

figure('Name','Detection Count');
nn    = round((max(tc)-min(tc))/params.AASM_EpochLength);
[N,X] = hist(tc,nn);
tmpx  = zeros(1,length(t));
tmpx(clusterIX) = 1;
subplot(4,1,1),plot(t,tmpx,'.-'); axis tight; title('Automatic detection of slow oscillations'); ylabel('SO detected'); xlabel('time (sec)');
subplot(4,1,2),try hist(tc,nn); end; axis tight; xlim=get(gca,'xlim'); title('Aggregate detection count histogram');
xlabel(sprintf('time (%gsec bins)',params.AASM_EpochLength)); ylabel('count');
subplot(4,1,3),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
title('SO detection count with threshold for SWS based on AASM guidelines');
xlabel(sprintf('time (%gsec bins)',params.AASM_EpochLength)); ylabel('count');

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
tmp     = N > thresh;                         % 1 if count > thresh; else 0
if ~any(tmp)
  error('No intervals satisfied the acceptance criteria.');
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
if length(L) > 2
  I       = find(X(L(2:end)) - X(R(1:end-1)) < params.IntervalSelection_CombineIfLessThan);
  L(I+1)  = [];
  R(I)    = [];
end
% require that count exceed CountThreshold for > MinCrossingTime
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
  fld = 'pospeak';   tmp = {detections.(fld)}; tmp = cellfun(@(x)x(x>=L&x<=R),tmp,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
  fld = 'negpeak';   tmp = {detections.(fld)}; tmp = cellfun(@(x)x(x>=L&x<=R),tmp,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:}); 
  fld = 'pospeak_cluster_time_index';      keep = cellfun(@(x)(x>=L&x<=R),{detections.(fld)},'uniformoutput',false);
  tmp = cellfun(@(x,y)x(y),{detections.(fld)},  keep,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
  fld = 'pospeak_cluster_number'; 
  tmp = cellfun(@(x,y)(x(y)),{detections.(fld)},keep,'UniformOutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newdetections.(fld)},'uniformoutput',false); [newdetections(1:nchan).(fld)] = deal(tmp{:});
  fld = 'negpeak_cluster_time_index';      keep = cellfun(@(x)(x>=L&x<=R),{detections.(fld)},'uniformoutput',false);
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
proc_clustered_detections_SWS = detections;
toc
%%

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
subplot(4,1,4),plot(X,NN,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r'); set(gca,'xlim',xlim);
title('AASM guideline-based automatic N3 detection');
xlabel(sprintf('time (%gsec bins)',params.AASM_EpochLength)); ylabel('count');

params.t0  = t0; % start time of each interval of interest
params.tf  = tf; % end time of each interval of interest
% total_toilim = [params.t0(1) params.tf(end)];

clear newdetections tmp L R k fld keep nInterval t nchan    

BegPoints       = [1 seldat.epochs.IntervalEndPoints(1:end-1)+1]; % indices
EndPoints       = seldat.epochs.IntervalEndPoints;                % indices
IntervalT0      = seldat.epochs.time(BegPoints);                  % sec
IntervalTf      = seldat.epochs.time(EndPoints);                  % sec
% The events structure can be used by visualizer to mark detections  
t       = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
events  = [];
for k   = 1:length(detections)
  events(k).label = data.sensor_info(k).label;
  events(k).time  = [t([detections(k).pospeak detections(k).negpeak])];
  events(k).type  = [1*ones(1,length(detections(k).pospeak)) 2*ones(1,length(detections(k).negpeak))];
end

toc


%% Compare results

vars  = {'clustered_detections','proc_clustered_detections','proc_clustered_detections_SWS'};
ns    = length(vars);
figure('Name','Detection Count'); 
for s = 1:ns
  det = eval(vars{s});
  t         = det(1).tstart:1/det(1).sfreq:det(1).tstop;
  typestr   = ''; % {'','pospeak_','negpeak_'} where ''=>both
  cnumfield = sprintf('%scluster_number',typestr);
  cindfield = sprintf('%scluster_time_index',typestr);
  clusterID = arrayfun(@(x)(x.(cnumfield)),det,'UniformOutput',false);
  tmpID     = [clusterID{:}];
  clusterID = unique(tmpID); % sorted cluster index, 1:nclusters
  clusterIX = arrayfun(@(x)(x.(cindfield)),det,'UniformOutput',false);
  clusterIX = [clusterIX{:}];
  clusterIX = cellfun(@(x)unique(clusterIX(tmpID==x)),num2cell(unique(tmpID)));
  nclusters = length(clusterID);
  cnum      = {det.(cnumfield)};
  cnum      = [cnum{:}];
  Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));
  % create cell array listing channels involved in each cluster
  cnum          = {det.(cnumfield)};
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
  subplot(4,ns,s+0*ns),plot(t,tmpx,'.-'); axis tight; title('aggregate cluster count')
  subplot(4,ns,s+1*ns),try hist(tc,nn); end; axis tight
  subplot(4,ns,s+2*ns),plot(X,N,'.-'); axis tight; hline(params.IntervalSelection_CountThreshold,'r');
  xlabel('time (sec)'); ylabel(sprintf('count (%gsec bins)',params.AASM_EpochLength));
  for k=1:length(t0f), vline(t0f(k),'k'); end
  subplot(4,ns,s+3*ns),plot(tc,'.-'); ylabel('cluster time (sec)'); xlabel('cluster number'); axis tight;
  clear nn N X tmpx
end

toc

save('tmpevents4.mat','events','clustered_detections','proc_clustered_detections','proc_clustered_detections_SWS','params');
