function [peaks,count,count_time,HistTime,CountThresh] = SO_cluster_detections(peaks,varargin)
% peaks = select_peakpairs(peaks,tau)
% Purpose: return only pos & neg peaks separated by less than tau [sec]
% Required inputs:
%   peaks - output of SO_detection()
% 
% Optional inputs:
%   method - {'aggregate','overlap','refchan'} (default = 'aggregate')
%
% Additional inputs depend on method:
%   aggregate:
%     w,s,tau
%   overlap:
%     tau
%   refchan:
%     tau
%
% Output:
%   peaks structure containing peaks involved in clusters
%   count: aggregate detection count histogram
%   count_time: times at each count step
%   HistTime: histogram-based cluster center times
%
% Created by JSS on 07-Jul-2010
% Multimodal Imaging Laboratory
% Department of Radiology, UCSD
%
% Modified for Paper I analysis on 27-Jul-2010

tic

parms = mmil_args2parms( varargin, ...
                         { 'method',            'histogram',      {'histogram','aggregate','overlap','refchan'},...
                           'threshold',         'meanstd',        [],...
                           'ClusterWindow',     .2,               [],...
                           'IntegrationWindow', .1,               [],...
                           'MinSeparation',     .1,               [],...
                           'StepSize',          1/peaks(1).sfreq, [],...
                           'peaktype',          'both',           [],...
                           'count',             [],               [],...
                         }, ...
                         false );
                       
% peaktypes: {'both','pospeak','negpeak','all',{'pospeak','negpeak'}}
tau               = parms.ClusterWindow;                       
method            = parms.method;
thresh            = parms.threshold;
StepSize          = parms.StepSize;
MinSeparation     = parms.MinSeparation;
IntegrationWindow = parms.IntegrationWindow;

% common parameters
Fs      = peaks(1).sfreq;
t       = peaks(1).tstart:1/Fs:peaks(1).tstop;
nchan   = length(peaks);
[sel(1:nchan).pospeak] = deal([]);
[sel(1:nchan).negpeak] = deal([]);

% if strcmp(parms.peaktype,'both') || strcmp(parms.peaktype,'all')
if any(strmatch('both',parms.peaktype)) || any(strmatch('all',parms.peaktype))
  parms.peaktype = {'pospeak','negpeak'};
elseif ~iscell(parms.peaktype)
  parms.peaktype = {parms.peaktype};
end

if ~isempty(parms.count)
  count = parms.count;
  L     = floor(IntegrationWindow*Fs/2);
  n     = floor(StepSize*Fs);
  nsmp  = length(t);
  c     = L+1:n:nsmp-L;
  tc    = t(c);
end

switch method
  case {'histogram','aggregate'}
    % calculate sliding window cross-channel aggregate detection count
    fprintf('calculating aggregate detection count...\n')
    if isempty(parms.count)
      [count,tc,c]  = calc_agg_detection_count(peaks,IntegrationWindow,StepSize,parms.peaktype);
      % count = windowed aggregate detection count at each step
      % tc    = interpolated time at each count step
      % c     = map from count-step to data-index
      count_time = tc;
    else
      count_time = tc;
    end
    % determine interval start/stop of each peak > threshold & max therein
    fprintf('finding peaks and defining reference times...\n')
    [Ic,tc,CountThresh] = find_peaks(count,thresh,tc,round(MinSeparation*Fs),MinSeparation);
      % Ic = count-step index for each epoch
      % tc = center time of each epoch
      HistTime = tc;
    % convert count-step index to data index
    refpeaks = c(Ic);
%     figure; subplot(3,1,1),plot(tc);title('tc');axis tight;subplot(3,1,2),plot(t(refpeaks)),title('t(Ipeaks)');axis tight;subplot(3,1,3),plot(1000*(tc-t(refpeaks)));title('diff');axis tight;
    % find all channels involved in each detection-count peak
    peaks = proc_involved_chans(refpeaks,peaks,round(tau*Fs/2));
  case 'overlap'
    % define clusters based on overlap of windows centered on detection
    % times across channels
    
  case 'refchan'
    % define clusters based on window around detection times in a reference
    % channel
    
end

% toc

function [count,tc,c] = calc_agg_detection_count(peaks,w,s,peaktype)
% w     = .1;                         % sec, integration window size
% s     = 1/peaks(1).sfreq;           % sec, step size
Fs    = peaks(1).sfreq;
t     = peaks(1).tstart:1/Fs:peaks(1).tstop;
pks   = [];
for k = 1:length(peaktype)
  tmp = {peaks.(peaktype{k})};
  pks = [pks [tmp{:}]];
end
pks   = sort(pks);        % peak indices (sorted)
tk    = t(pks);           % peak times (sorted)
L     = floor(w*Fs/2);    % padding around window centers
n     = floor(s*Fs);      % step size in indices
nsmp  = length(t);        % # of time points
npks  = length(tk);       % # of peaks
c     = L+1:n:nsmp-L;     % window centers in indices
tc    = t(c);             % window centers in sec
nstep = length(c);
count = zeros(1,nstep);
pkbin      = zeros(1,length(t));
pkbin(pks) = 1;
count      = cellfun(@(x)sum(pkbin(x-L:x+L)),num2cell(c));

function [I,tc,threshold] = find_peaks(count,thresh,tc,MinSeparation,MinLength)
% MinLength:
  % - combine two threshold crossings if they are closer than MinSeparation [index]
  % - discard any threshold crossings which are shorter than MinLength [s]
  
if isnumeric(thresh)
  threshold   = thresh;
elseif isequal(thresh,'mean')
  threshold   = mean(count);
elseif isequal(thresh,'meanstd')
  threshold   = mean(count) + std(count);
elseif isequal(thresh,'meanstd2')
  threshold   = mean(count) + 2*std(count);
elseif isequal(thresh,'meanstd3')
  threshold   = mean(count) + 3*std(count);
elseif isequal(thresh,'meanstd4')
  threshold   = mean(count) + 4*std(count);
elseif isequal(thresh,'meanstd5')
  threshold   = mean(count) + 5*std(count);  
elseif isequal(thresh,'median')
  threshold   = median(count);
elseif isequal(thresh,'medianstd')
  threshold   = median(count) + std(count);
elseif isequal(thresh,'medianstd2')
  threshold   = median(count) + 2*std(count);
elseif isequal(thresh,'medianstd3')
  threshold   = median(count) + 3*std(count);
elseif isequal(thresh,'medianstd4')
  threshold   = median(count) + 4*std(count);
elseif isequal(thresh,'medianstd5')
  threshold   = median(count) + 5*std(count);
end
fprintf('Aggregate detection count threshold: %g\n',threshold);
if threshold > max(count)
  error('Threshold exceeds maximum count');
end
ExceedThresh  = count > threshold; % 1's where agg > threshold, else 0
L = find(ExceedThresh(2:end)==1 & ExceedThresh(1:end-1)==0); % index to left edge
R = find(ExceedThresh(1:end-1)==1 & ExceedThresh(2:end)==0); % index to right edge
if L(1)-R(1)==1 % correct for right edge at beginning
  L = L(1:end-1);
  R = R(2:end);
end
% combine threshold crossings if right edge - left edge < MinLength
sel       = find((L(2:end)-R(1:end-1)) <= MinSeparation);
L(sel+1)  = [];
R(sel)    = [];
% find center time & index for each remaining threshold crossing
T = tc(2:end-1);
D = T(R) - T(L); 
I   = find(D >= MinLength);
t0  = T(L(I));
tf  = T(R(I));
tc  = mean([t0;tf],1);
I   = round(mean([L(I);R(I)],1));
% I   = I + 1; % shift b/c tc(2:end-1)
% figure; 
% try hist(tc,round(length(tc)/10)); end; 
% try axis tight; title('detection histogram'); catch close; end

function inv = proc_involved_chans(tk,peaks,tau)
% involved if either pos or neg peak in tk+/-tau
% tk  = indices to peaks in the detection histogram
% tau = padding around histogram peaks used to define windows
inv = peaks;
pos = {peaks.pospeak};                      % indices to peaks in each chan
neg = {peaks.negpeak};
% loop over channels & find all peaks within window around histogram peak
for k = 1:length(pos)
  ix1 = cellfun(@(x)(find((pos{k}>=x-tau)&(pos{k}<=x+tau))),num2cell(tk),'uniformoutput',false);
    % for each histpeak: indices of pospeaks within window
      % a chan may have no peaks in window, 1, or more
    % for this chan, only consider involved clusters
  s1  = find(~cellfun(@isempty,ix1));     % cluster #'s for involved clusters 
                                          %  = indices to peaks involved in clusters
                                          % (s1 is a list of involved clusters)                                          
  ix1 = ix1(s1);                          % ix1 is now a list of indices to peaks involved in clusters
  tmp = cellfun(@length,ix1);
  % if multiple detections fall in the cluster window, select the detection
  % closest to the histogram reference time
  if any(tmp>1)
    tmp = find(tmp>1);
    ix1(tmp) = cellfun(@(x){ix1{x}(nearest(pos{k}(ix1{x}),tk(x)))},num2cell(tmp));
  end
  % if a detection is currently involved in multiple clusters, remove it
  % from all clusters accept that which has a hist ref time closest to it
  tmp = [ix1{:}];
  tmp = find(diff(tmp)==0);   
    % index to a histpeak followed by a second peak involving the same
    % chanpeak => ix1(tmp)==ix1(tmp+1)
  % is histpeak at tmp(1) or tmp(1)+1 closer to chanpeak ix1(tmp(1))??
  if ~isempty(tmp)
    if tmp(end) >= length(ix1), tmp = tmp(tmp<length(ix1)); end
    % tmp       = cluster histpeak indices (histpeak)
    % ix1(tmp)  = channel peak indices (wavepeak)
    % is the 1st or 2nd cluster histpeak closer to the chan pospeak?
    clusterchoice = cellfun(@(x,y)nearest(tk(s1([x x+1])),pos{k}(y)),num2cell(tmp),ix1(tmp));
    % remove the opposite 
    rmix = [];%zeros(1,length(tmp));
    rmix(clusterchoice==1) = tmp(clusterchoice==1) + 1; % remove 2nd
    rmix(clusterchoice==2) = tmp(clusterchoice==2);     % remove 1st
    s1(rmix)  = [];
    ix1(rmix) = [];
    clear clusterchoice rmix
  end
  % repeat for negative half-wave peaks detected in this channel
  ix2 = cellfun(@(x)(find((neg{k}>=x-tau)&(neg{k}<=x+tau))),num2cell(tk),'uniformoutput',false);
  s2  = find(~cellfun(@isempty,ix2));
  ix2 = ix2(s2);
  tmp = cellfun(@length,ix2);
  % if multiple detections fall in the cluster window, select the detection
  % closest to the histogram reference time
  if any(tmp>1)
    tmp = find(tmp>1);
    ix2(tmp) = cellfun(@(x){ix2{x}(nearest(neg{k}(ix2{x}),tk(x)))},num2cell(tmp));
  end
  % if a detection is currently involved in multiple clusters, remove it
  % from all clusters accept that which has a hist ref time closest to it
  tmp = [ix2{:}];
  tmp = find(diff(tmp)==0);     
  if ~isempty(tmp)
    if tmp(end) >= length(ix2), tmp = tmp(tmp<length(ix2)); end
    clusterchoice = cellfun(@(x,y)nearest(tk(s2([x x+1])),neg{k}(y)),num2cell(tmp),ix2(tmp));
    rmix = [];%zeros(1,length(tmp));
    rmix(clusterchoice==1) = tmp(clusterchoice==1) + 1; % remove 2nd
    rmix(clusterchoice==2) = tmp(clusterchoice==2);     % remove 1st
    s2(rmix)  = [];
    ix2(rmix) = [];
    clear clusterchoice rmix
  end
  inv(k).pospeak = pos{k}(unique([ix1{:}]));
  inv(k).negpeak = neg{k}(unique([ix2{:}]));
  inv(k).pospeak_cluster_number     = s1;
  inv(k).pospeak_cluster_time_index = tk(s1);
  inv(k).negpeak_cluster_number     = s2;
  inv(k).negpeak_cluster_time_index = tk(s2);
  inv(k).cluster_number             = unique([s1 s2]);
  inv(k).cluster_time_index         = tk(unique([s1 s2]));
  clear s1 s2 ix1 ix2
end
% keyboard




