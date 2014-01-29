% Purpose: this script applies additional constraints on the detection
% clusters. It requires that clusters have at least some minimum fraction
% (MinChansPerCluster) of the total number of channels being processed. 
% Then it creates a low resolution histogram of the detection center times
% which is used to define time intervals of interest as those satisfying 
% the following criteria:
% 1. Detection count > params.IntervalSelection_CountThreshold
% 2. Combine if separated by less than params.IntervalSelection_CombineIfLessThan
% 3. Duration of threshold crossing > params.IntervalSelection_MinCrossingTime
%
% Created by Jason Sherfey on 27-Jul-2010
% Multimodal Imaging Laboratory
% Department of Radiology, UCSD

% find # of channels involved in each cluster
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
clear typestr

% % Require that some % of the good channels are present in each cluster
% Nthresh   = round(params.MinChansPerCluster*length(detections));
% ckeep     = find(Ninvolved > Nthresh);
% % Eliminate clusters with (# involved channels) < threshold
% Ninvolved = Ninvolved(ckeep);
% clusterID = clusterID(ckeep);
% clusterIX = clusterIX(ckeep);
%   % clusterID now contains only the clusters to process

% Create cell array listing channels involved in each cluster
% cnum          = {detections.(cnumfield)};
% tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
% tmp           = [tmp{:}]';
% tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
% InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);
%
% Check that the number of chans in each cluster is correct:
% all(bsxfun(@eq,Ninvolved,cellfun(@length,InvolvedChans)))
% clear InvolvedChans

clear cnumfield cindfield clusterID tmpID cnum Nthresh ckeep

% center times for each cluster
t           = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;
tc          = t(clusterIX);
nn = round((max(tc)-min(tc))/params.AASM_EpochLength);
% nn          = ceil(length(tc)/10);
[N,X]       = hist(tc,nn);
% [cidx,ctrs] = kmeans(N,2); ix = find(ctrs==max(ctrs));

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
if L(1)-R(1)==1                               % correct for right edge at beginning
  L = L(1:end-1);
  R = R(2:end);
end
if R(1)   < L(1),   R = R(2:end);   end
if L(end) > R(end), L = L(1:end-1); end
if R(end-1) > L(end), R = R(1:end-1); end
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

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %  POTENTIAL FIGURE FOR PAPER
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Detection Count'); 
tmpx  = zeros(1,length(t));
tmpx(clusterIX) = 1;
subplot(3,1,1),plot(t,tmpx,'.-'); axis tight; title('Aggregate slow oscillation detection count')
subplot(3,1,2),try hist(tc,nn); end; axis tight
subplot(3,1,3),plot(X,N,'.-'); axis tight;%,X(cidx==ix),N(cidx==ix),'b.'); axis tight
hline(thresh,'r');xlabel('time (sec)'); ylabel(sprintf('count (%gsec bins)',params.AASM_EpochLength));
for k=1:length(t0f), vline(t0f(k),'k'); end
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear clusterIX tmpx Ninvolved tc nn N X I t0f thresh

%% select data & detections within intervals of interest

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
params.t0  = t0; % start time of each interval of interest
params.tf  = tf; % end time of each interval of interest
% total_toilim = [params.t0(1) params.tf(end)];

clear newdetections tmp L R k fld keep t0 tf nInterval t nchan