function peaks = SO_detection(data,varargin)

parms = mmil_args2parms( varargin, ...
                         { 'bpfilter'         ,1,[],...
                           'blc'              ,0,[],...
                           'decimate'         ,0,{0,1},...
                           'smooth'           ,1,{0,1},...
                           'hilbertpeaks'     ,1,{0,1},...
                           'derivpeaks'       ,1,{0,1},...
                           'onlysinglepeaks'  ,0,{0,1},...
                           'zerocross'        ,1,{0,1},...
                           'monotonic'        ,1,{0,1},...
                           'gtmedian'         ,1,{0,1},...
                           'return_zerocross' ,0,{0,1},...
                           'bpfreq'           ,[.3 3],[],...
                           'bptrans'          ,[0 0],[],...
                           'blcwindow'        ,[],[],...                        
                           'decimate_factor'  ,2,[],...
                           'smooth_window'    ,.05,[],...
                           'smooth_method'    ,'moving',{'moving','lowess','loess','sgolay','rlowess','rloess'},...
                           'zero2zero_limits' ,[.25 1],[],...
                           'monotonic_thresh' ,[],[],...
                           'toilim'           ,[],[],...
                           'mindist'          ,[],[],...
                           'debug'            ,0,{0,1},...
                           'debug_toilim'     ,[],[],...
                           'min_abs_peak'     ,[],[],...
                         }, ...
                         false );

[datatype datafield] = ts_object_info(data);
if isequal(parms.bpfilter,'yes'), parms.bpfilter = 1; end
if isequal(parms.bpfilter,'no') , parms.bpfilter = 0; end
if isempty(parms.monotonic_thresh), parms.monotonic_thresh = parms.bpfilter(1); end          
if isempty(parms.toilim),   parms.toilim    = [data.(datafield)(1).time(1) data.(datafield)(1).time(end)]; end
if isempty(parms.mindist),  parms.mindist   = 2*min(parms.zero2zero_limits); end
if isempty(parms.blcwindow),parms.blcwindow = [data.(datafield)(1).time(1) data.(datafield)(1).time(end)]; end
if ischar(parms.blc), if strcmp(parms.blc,'yes'),parms.blc=1; else parms.blc=0; end; end
if ischar(parms.toilim) && strcmp(parms.toilim,'all'), parms.toilim = [data.(datafield)(1).time(1) data.(datafield)(1).time(end)]; end
if parms.debug
  % how many subplots?
  nsubplot = 1 + parms.bpfilter + parms.blc + parms.decimate + parms.smooth + ...
    parms.hilbertpeaks + parms.derivpeaks + parms.zerocross + parms.monotonic + ...
    parms.gtmedian + parms.mindist;
    ncol = ceil(sqrt(nsubplot));
    nrow = (ncol-1) + double(floor(ncol*(ncol-1)/nsubplot)==0);
  if isempty(parms.debug_toilim)
    parms.debug_toilim = [min(data.(datafield)(1).time),min(data.(datafield)(1).time)+20];
  end
  debug_toilim = parms.debug_toilim;
  figure('tag','debugplots');
end
sens  = {data.sensor_info.label};
nchan = data.num_sensors;

% % bpfilter
% if parms.bpfilter
%   dat = ts_freq_filt(data.(datafield).data',Fs,parms.bpfreq,parms.bptrans,'bandpass');
% end
tstart = tic;
% loop over channels
for k = 1:nchan
  fprintf('processing channel %g of %g (%g min)\n',k,nchan,toc(tstart)/60);
  if parms.debug, plotcnt = 0; set(findobj('tag','debugplots'),'Name',sens{k}); end
  
  % select data
  Fs  = data.sfreq;
  dat = ts_data_selection(data,'chanlabel',sens{k},'toilim',parms.toilim,'verbose',0);
  t   = dat.(datafield).time;
  x   = dat.(datafield).data;
  
  if parms.debug
    plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
    tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
    plot(tt,xx,'k-'); axis tight; title('raw data'); hline(0,'k');
  end
  
  % bandpass filter
  if parms.bpfilter
    if isa(x,'single')
      x = double(x);    
      x = ts_freq_filt(x,Fs,parms.bpfreq,parms.bptrans,'bandpass');
      x = single(x);
    else
      x = ts_freq_filt(x,Fs,parms.bpfreq,parms.bptrans,'bandpass');
    end
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      plot(tt,xx,'k-'); axis tight; title(sprintf('bandpass (%g-%gHz)',parms.bpfreq)); hline(0,'k');
    end
  end
  
  % baseline correction
  if parms.blc
    tix = nearest(t,parms.blcwindow(1)):nearest(t,parms.blcwindow(2));
    x   = x - mean(x(tix));
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      plot(tt,xx,'k-'); axis tight; title(sprintf('baseline subtraction (%g-%gsec)',parms.blcwindow)); hline(0,'k');
    end    
  end
  
  % decimate
  if parms.decimate
    if isa(x,'single')
      x = double(x);
      x = decimate(x,parms.decimate_factor);
      x = single(x);
    else
      x = decimate(x,parms.decimate_factor);
    end
    Fs = Fs / parms.decimate_factor;
    t  = downsample(t,parms.decimate_factor);
    if length(x) ~= length(t), error('Something went wrong with the decimation'); end
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      plot(tt,xx); axis tight; title(sprintf('decimate (%gHz/%g)',data.sfreq,parms.decimate_factor)); hline(0,'k');
    end    
  end
  
  % smooth (moving average)
  if parms.smooth
    x = smooth(x,ceil(parms.smooth_window*Fs),parms.smooth_method);
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      plot(tt,xx,'k-'); axis tight; title(sprintf('smooth (%gsec)',parms.smooth_window)); hline(0,'k');
    end    
  end
  
  % Hilbert
  if parms.hilbertpeaks
    % pos/neg peaks (phi zero-crossings)
    tmp = angle(hilbert(x));
    ind = crossing(tmp);
    ind(ind==1 & ind==length(t)) = [];
    % pos if pt left is < 0 and right is > 0 (increase from neg to pos)
    ind = ind(ind>1 & ind<length(tmp));
    pos = ind((tmp(ind-1)<0) & tmp(ind+1)>0);
    % neg if pt left is > 0 and right is < 0 (jump from +pi to -pi)
    neg = ind((tmp(ind-1)>0) & tmp(ind+1)<0);
    clear tmp ind
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      pix=find(pos>tix(1)&pos<tix(end));nix=find(neg>tix(1)&neg<tix(end));
      plot(tt,xx,'k-',t(pos(pix)),x(pos(pix)),'b*',t(neg(nix)),x(neg(nix)),'r*'); axis tight
      title('hilbert peaks'); hline(0,'k');
    end    
  end
  
  % Find peaks as zero-crossings of the signal derivative
  if parms.derivpeaks
    if parms.hilbertpeaks
      % refine peak times from hilbert detection method
      % select peak with max amplitude b/w adjacent signal zero crossings
      pos = sel_derivpeaks(pos,x,parms.onlysinglepeaks,'pos');
      neg = sel_derivpeaks(neg,x,parms.onlysinglepeaks,'neg');
    else
      % find all peaks for the first time
      tmp = diff(x);
      ind = crossing(tmp);
      ind(ind==1 & ind==length(t)) = [];      
      ind = ind(ind>1 & ind<length(tmp));
      % pospeak if crosses zero from pos to neg
      pos = ind((tmp(ind-1)>0) & tmp(ind+1)<0);
      % negpeak if crosses zero from neg to pos
      neg = ind((tmp(ind-1)<0) & tmp(ind+1)>0);
    end
    clear tmp ind
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      pix=find(pos>tix(1)&pos<tix(end));nix=find(neg>tix(1)&neg<tix(end));
      plot(tt,xx,'k-',t(pos(pix)),x(pos(pix)),'b*',t(neg(nix)),x(neg(nix)),'r*'); axis tight
      if parms.hilbertpeaks,title('derivative refinement'),else title('derivative peaks'); end; hline(0,'k');
    end    
  end
    
  % Zero-to-zero criterion
  if parms.zerocross
    pos = sel_zero2zero(pos,x,t,parms.zero2zero_limits);
    neg = sel_zero2zero(neg,x,t,parms.zero2zero_limits);
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      pix=find(pos>tix(1)&pos<tix(end));nix=find(neg>tix(1)&neg<tix(end));
      plot(tt,xx,'k-',t(pos(pix)),x(pos(pix)),'b*',t(neg(nix)),x(neg(nix)),'r*'); axis tight
      title(sprintf('zerocross constraint (%g-%gsec)',parms.zero2zero_limits)); hline(0,'k');
    end    
  end
  
  % Monotonic phase criterion
  if parms.monotonic
    pos = sel_monotonic(pos,x);
    neg = sel_monotonic(neg,-x);
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      pix=find(pos>tix(1)&pos<tix(end));nix=find(neg>tix(1)&neg<tix(end));
      plot(tt,xx,'k-',t(pos(pix)),x(pos(pix)),'b*',t(neg(nix)),x(neg(nix)),'r*'); axis tight
      title('monotonic constraint'); hline(0,'k');
    end    
  end
  
  % Minimum distance
  if parms.mindist ~= 0
    pos = rm_neighbours(pos,x,ceil(parms.mindist*Fs));
    neg = rm_neighbours(neg,x,ceil(parms.mindist*Fs));
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      pix=find(pos>tix(1)&pos<tix(end));nix=find(neg>tix(1)&neg<tix(end));
      plot(tt,xx,'k-',t(pos(pix)),x(pos(pix)),'b*',t(neg(nix)),x(neg(nix)),'r*'); axis tight
      title(sprintf('min dist constraint (%gsec)',parms.mindist)); hline(0,'k');
    end    
  end
  
  % Amplitude criterion
  if parms.gtmedian || (~isempty(parms.min_abs_peak) && strcmp('eeg',dat.sensor_info.typestring))
    if parms.gtmedian
      [pos,neg] = sel_gtmedian(x,pos,neg);
    end
    if ~isempty(parms.min_abs_peak) && strcmp('eeg',dat.sensor_info.typestring)
      pos = pos(x(pos)>parms.min_abs_peak);
      neg = neg(abs(x(neg))>parms.min_abs_peak);
    end
    if parms.debug
      plotcnt = plotcnt+1; subplot(nrow,ncol,plotcnt)
      tix=nearest(t,debug_toilim(1)):nearest(t,debug_toilim(2));xx=x(tix);tt=t(tix);
      pix=find(pos>tix(1)&pos<tix(end));nix=find(neg>tix(1)&neg<tix(end));
      plot(tt,xx,'k-',t(pos(pix)),x(pos(pix)),'b*',t(neg(nix)),x(neg(nix)),'r*'); axis tight
      title('amplitude constraint (>median)'); hline(0,'k');
    end    
  end
    
  if parms.return_zerocross
    fprintf('Not returning adjacent zero-crossings; incomplete implementation\n');
  end

  peaks(k).label    = dat.sensor_info.label;
  peaks(k).pospeak  = pos;
  peaks(k).negpeak  = neg;
  peaks(k).sfreq    = Fs;
  peaks(k).tstart   = t(1);
  peaks(k).tstop    = t(end);
  
  clear tmp pos neg x phi frq
  if parms.debug, pause; end
end

%% SUBFUNCTIONS

function pks = sel_derivpeaks(pks,x,singlepeaks,peaktype)
x0  = crossing(x);            % signal zero crossings
dx0 = crossing(gradient(x));  % signal derivative zero crossing
% group left/right signal zero crossing indices
L   = cellfun(@(x)(x0(find(x0<=x,1,'last'))),num2cell(pks),'UniformOutput',false);
R   = cellfun(@(x)(x0(find(x0>=x,1,'first'))),num2cell(pks),'UniformOutput',false);
% select peaks with adjacent zero crossings
sel = ~(cellfun(@isempty,L) | cellfun(@isempty,R));
pks = pks(sel);
L   = L(sel);
R   = R(sel);
% get indices to all peaks (derivative zero crossings) between the signal zero crossings
groups = cellfun(@(a,b)(dx0(dx0>=a & dx0<=b)),L,R,'UniformOutput',false);
groups = groups(~cellfun(@isempty,groups));
% remove detections with multiple peaks
if singlepeaks
  sel  = cellfun(@length,groups)==1;
  pks  = pks(sel);
else
  % choose peak between zero crossings
  if strcmp(peaktype,'neg')
    % select min peak
    selind = cellfun(@(y)(min(x(y))==x(y)),groups,'UniformOutput',false);
    selind = selind(~cellfun(@isempty,selind));
    tmp    = cellfun(@(y,z)(y(z)),groups,selind,'uniformoutput',false);
    pks    = unique(cellfun(@(x)(x(1)),tmp));    
%     pks    = unique(cellfun(@(y,z)(y(z)),groups,selind));
    % select only peaks < 0
    pks    = pks(x(pks)<0);
  elseif strcmp(peaktype,'pos')
    % select max peak
    selind = cellfun(@(y)(max(x(y))==x(y)),groups,'UniformOutput',false);
    groups = groups(~cellfun(@isempty,selind));
    selind = selind(~cellfun(@isempty,selind));
    tmp    = cellfun(@(y,z)(y(z)),groups,selind,'uniformoutput',false);
    tmp    = cellfun(@(y,z)(y(z)),groups,selind,'uniformoutput',false);
    pks    = unique(cellfun(@(x)(x(1)),tmp));    
%     pks    = unique(cellfun(@(y,z)(y(z)),groups,selind));
    % select only peaks > 0
    pks    = pks(x(pks)>0);
  end
end

function pks = sel_zero2zero(pks,x,t,lims)
a = lims(1);
b = lims(2);
% require time b/w zero-crossings tau: a<=tau<=b
% find zero crossings in the signal
[ind,t0] = crossing(x,t);
% find signal zero-crossings before each peak
L = cellfun(@(x)(find(ind<x,1,'last')),num2cell(pks),'UniformOutput',false);
% find signal zero-crossings after each peak
R = cellfun(@(x)(find(ind>x,1,'first')),num2cell(pks),'UniformOutput',false);
% Select peaks with adjacent zero-crossings
sel = ~(cellfun(@isempty,L) | cellfun(@isempty,R));
pks = pks(sel);
L   = t(ind([L{sel}])); % convert indices to seconds
R   = t(ind([R{sel}])); % convert indices to seconds
% Calculate time between zero-crossings
tau = R - L;
% Apply constraint
sel = (tau>=a) & (tau<=b);
pks = pks(sel);

function pks = sel_monotonic(pks,x)
% select strictly monotonic phase runs
phi = angle(hilbert(x));
% find pi-pi/12 & -pi+pi/12 crossings
% simplify by working with the absolute value (just pi-pi/12)
absphi = abs(phi);
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

function [pos,neg] = sel_gtmedian(x,pos,neg)
% type = [ones(1,length(pos)) 2*ones(1,length(neg))];
% pks  = [pos neg];
% tmp  = abs(x(pks));
% xmu  = median(tmp);
% sel  = find(tmp > xmu);
% pks  = pks(sel);
% type = type(sel);
% pos  = pks(type==1);
% neg  = pks(type==2);

tmp = x(pos);
xmu = median(tmp);
pos = pos(tmp>xmu);

tmp = x(neg);
xmu = median(tmp);
neg = neg(tmp<xmu);

function pks = rm_neighbours(pks,x,tau)
% find pks within tau of each other
groups = cellfun(@(y)(find(abs(pks-y)<tau)),num2cell(pks),'UniformOutput',false);
maxind = cellfun(@(y)(max(x(y))==x(y)),groups,'UniformOutput',false);
groups = groups(~cellfun(@isempty,maxind));
maxind = maxind(~cellfun(@isempty,maxind));
tmp    = cellfun(@(y,z)(y(z)),groups,maxind,'uniformoutput',false);
pks    = pks(unique(cellfun(@(x)(x(1)),tmp)));





