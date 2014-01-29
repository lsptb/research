function [detections,params] = SWS_detection(detections,varargin)
% Purpose: find slow-wave sleep given slow wave detections
%
% Input:
%   - detections: output of SO_cluster_detections()
%
% Output:
%   - detections: only detection clusters during SWS
%   - params structure with two new fields:
%       t0 = start times for all intervals of interest
%       tf = end   times for all intervals of interest
% 
% Recommended function sequence:
%   SO_detection()
%   SO_cluster_detections()
%   SWS_detection()
%
% Created by Jason Sherfey on 08-Nov-2010

params = mmil_args2parms( varargin, ...
                         { 'AASM_EpochLength',30,[],...
                           'IntervalSelection_CountThreshold',5,[],...
                           'IntervalSelection_CombineIfLessThan',45,[],...
                           'IntervalSelection_MinCrossingTime',60,[],...
                           'chantype',[],[],...
                           'ShowPlots',0,[],...
                           'fid',1,[],...
                           },false);
            
if isempty(params.chantype)
  if strmatch('MEG',detections(1).label)
    params.chantype = 'grad';
  else
    params.chantype = 'eeg';
  end
end

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
fprintf(params.fid,'\tCount SO clusters in %gsec epochs.\n',params.AASM_EpochLength);
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
fprintf(params.fid,'\tLabel an epoch as slow-wave sleep if it has more than %g clusters.\n',thresh);
  % threshold should be doubled for EEG
tmp     = N > thresh;                         % 1 if count > thresh; else 0
if ~any(tmp)
  fprintf(params.fid,'Error: No intervals satisfy the acceptance criteria.\n');
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
fprintf(params.fid,'\tPatch: relabel a period below threshold for <= %g epoch(s) if it is bounded by slow-wave epochs.\n',floor(params.IntervalSelection_CombineIfLessThan/params.AASM_EpochLength));
if length(L) > 2
  I       = find(X(L(2:end)) - X(R(1:end-1)) < params.IntervalSelection_CombineIfLessThan);
  L(I+1)  = [];
  R(I)    = [];
end
% require that count exceed CountThreshold for > MinCrossingTime
fprintf(params.fid,'\tDefine a multi-epoch interval as slow-wave sleep if it contains at least %g slow-wave epochs.\n',ceil(params.IntervalSelection_MinCrossingTime/params.AASM_EpochLength));
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
%   if k == 1
%     seldat = ts_data_selection(data,'toilim',[t0(k) tf(k)]);
%     seldat.epochs.IntervalEndPoints(k) = size(seldat.epochs.data,2);
%   else
%     tmpdat = ts_data_selection(data,'toilim',[t0(k) tf(k)]);
%     seldat.epochs.time = [seldat.epochs.time tmpdat.epochs.time];
%     seldat.epochs.data = cat(2,seldat.epochs.data,tmpdat.epochs.data);
%     seldat.epochs.IntervalEndPoints(k) = size(seldat.epochs.data,2);
%     clear tmpdat
%   end
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
for k = 1:length(t0)
  fprintf(params.fid,'t = %g - %g sec (%3.3g min)\n',t0(k),tf(k),(tf(k)-t0(k))/60);
end
% total_toilim = [params.t0(1) params.tf(end)];
% BegPoints       = [1 seldat.epochs.IntervalEndPoints(1:end-1)+1]; % indices
% EndPoints       = seldat.epochs.IntervalEndPoints;                % indices
% IntervalT0      = seldat.epochs.time(BegPoints);                  % sec
% IntervalTf      = seldat.epochs.time(EndPoints);                  % sec