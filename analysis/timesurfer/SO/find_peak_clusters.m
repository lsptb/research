function [peaks,count] = find_peak_clusters(peaks,varargin)
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
%
% Created by JSS on 07-Jul-2010
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

if strcmp(parms.peaktype,'both') || strcmp(parms.peaktype,'all')
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
    end
    % determine interval start/stop of each peak > threshold & max therein
    fprintf('finding peaks and defining reference times...\n')
    [Ic,tc]       = find_peaks(count,thresh,tc,round(MinSeparation*Fs),MinSeparation);
      % Ic = count-step index for each epoch
      % tc = center time of each epoch
    % convert count-step index to data index
    refpeaks = c(Ic);
    figure; subplot(3,1,1),plot(tc);title('tc');axis tight;subplot(3,1,2),plot(t(refpeaks)),title('t(Ipeaks)');axis tight;subplot(3,1,3),plot(1000*(tc-t(refpeaks)));title('diff');axis tight;
    % find all channels involved in each detection-count peak
    peaks = proc_involved_chans(refpeaks,peaks,round(tau*Fs/2));
  case 'overlap'
    % define clusters based on overlap of windows centered on detection
    % times across channels
    
  case 'refchan'
    % define clusters based on window around detection times in a reference
    % channel
    
end

toc

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
pts   = 1:floor(nstep/10):nstep;
tstart= tic;
for k = 1:nstep
  if ismember(k,pts)
    fprintf('progress: report %g of %g (step %g of %g), %g min\n',find(k==pts),length(pts),k,nstep,toc(tstart)/60);
  end
  count(k) = length(find((pks>=c(k)-L)&(pks<=c(k)+L)));
end

function [I tc] = find_peaks(count,thresh,tc,MinSeparation,MinLength)
% MinLength:
  % - combine two threshold crossings if they are closer than MinSeparation [index]
  % - discard any threshold crossings which are shorter than MinLength [s]
  
if 1%strcmp(thresh,'meanstd')
  threshold   = mean(count) + std(thresh);
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
figure; try hist(tc,round(length(tc)/10)); end; axis tight; title('detection histogram')

function inv = proc_involved_chans(tk,peaks,tau)
% involved if either pos or neg peak in tk+/-tau
inv = peaks;
pos = {peaks.pospeak};
neg = {peaks.negpeak};
% loop over channels & find all peaks within window around histogram peak
for k = 1:length(pos)
  ix1 = cellfun(@(x)(find((pos{k}>=x-tau)&(pos{k}<=x+tau))),num2cell(tk),'uniformoutput',false);
  s1  = find(~cellfun(@isempty,ix1)); % % indices of clusters of which this chan is involved
  ix1 = ix1(s1);
  ix2 = cellfun(@(x)(find((neg{k}>=x-tau)&(neg{k}<=x+tau))),num2cell(tk),'uniformoutput',false);
  s2  = find(~cellfun(@isempty,ix2));
  ix2 = ix2(s2);
  sel = unique([s1 s2]);
  inv(k).pospeak = pos{k}(unique([ix1{:}]));
  inv(k).negpeak = neg{k}(unique([ix2{:}]));
  inv(k).pospeak_cluster_number     = s1;
  inv(k).pospeak_cluster_time_index = tk(s1);
  inv(k).negpeak_cluster_number     = s2;
  inv(k).negpeak_cluster_time_index = tk(s2);
  inv(k).cluster_number             = unique([s1 s2]);
  inv(k).cluster_time_index         = tk(unique([s1 s2]));
end




