addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/scripts

%% 07-Jul-2010
% waveforms 
  Sleep_SO_MEG_ParamFile_s2
  if parms.clearall,clear all; Sleep_SO_MEG_ParamFile_s2; end
  SensorEventFile         = parms.SensorEventFile;
  SensorEventPhaseFile    = parms.SensorEventPhaseFile;
  FlipFile                = parms.FlipFile;
  CorrSensorEventFile     = parms.CorrSensorEventFile;
  ClusterSensorEventFile  = parms.ClusterSensorEventFile;
  datestring = '07-Jul-2010';%date;
  tstart     = tic;
  tstartorig = tstart;
  toilim     = parms.toilim; %[800 2400]
  writelog   = parms.writelog;
  layout     = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
  outpath    = [parms.rootoutdir '/' parms.SubjectID];
  FigPath    = sprintf('%s/images',outpath);
  logfile    = sprintf('%s/%s.log',outpath,date);
   
  subjects   = {'s1','s2','s3','s4','s5','s6','s7','s8'};
  subj       = strmatch(parms.SubjectID,subjects); % subj = 4
  subject    = subjects{subj};
  script_get_datafiles % fiffiles, badlabels
  
  if writelog
    fid = fopen(logfile,'a');
    fprintf(fid,'---------------------------------\n%s\nSlow oscillation analysis\n---------------------------------\n',datestr(now));
    if parms.clearall
      fprintf(fid,'Memory cleared.\n');
    end
  else
    fid = 1;
  end
  
  %% Load data (grad1 & grad2)
    loadflag = parms.loadflag; %   loadflag = 0;
    toiflag  = parms.toiflag;  %   toiflag  = 0;
    RSSflag  = parms.RSSflag;
  
  if loadflag
    % read fif files, convert data to timesurfer format, and save mat files
    matfiles = {};
    chantype = {'grad1','grad2'};
    findex   = parms.matfile_index;
    fprintf(fid,'Processing %s data\n',[chantype{:}]);
    for f = 1:length(fiffiles)
      fif = fiffiles{f};
      [fpath,fname,fext]  = fileparts(fif);
      outfile             = sprintf('%s/matfiles/%s_grad.mat',outpath,fname);
      matfiles{end+1}     = outfile;
      if exist(outfile,'file') % never overwrite (param independent)
        fprintf(fid,'MAT file already exists. not re-reading FIF: %s\n',fif);
        continue
      else
        fprintf(fid,'Reading FIF file: %s\n',fif);
      end
      data = ts_loadfif(fif,chantype,'epochs');
      fprintf(fid,'Saving MAT file: %s\n',outfile);
      save(outfile,'data');
      clear fif data
    end
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    % read mat files and combine data
    fprintf(fid,'Loading MAT files:\n');
    for k  = 1:length(findex),fprintf(fid,'%s\n',matfiles{findex(k)}); end
    data   = SO_combine_matfiles(matfiles(findex));
  end
  if toiflag
    if ischar(toilim) && strcmpi(toilim,'all')
      toilim = [data.epochs.time(1) data.epochs.time(end)];
      parms.toilim = toilim;
    end
    fprintf(fid,'Selecting times %g-%g sec\n',toilim);
    data   = ts_data_selection(data,'toilim',toilim);
    fprintf(fid,'done.\n');
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  if RSSflag
    fprintf(fid,'Calculating RSS from grad1 & grad2\n');
    sens      = data.sensor_info;
    nchan     = length(sens);
    grad1_ind = strmatch('grad1',{data.sensor_info.typestring});
    grad2_ind = strmatch('grad2',{data.sensor_info.typestring});
    nloc      = length(grad1_ind);
    for k = 1:nloc
      Bx        = data.epochs.data(grad1_ind(k),:);
      By        = data.epochs.data(grad2_ind(k),:);
      data.epochs.data(nchan+k,:) = (Bx.^2 + By.^2).^(1/2);
      sens(nchan+k)               = sens(grad1_ind(k));
      sens(nchan+k).label         = [data.sensor_info(grad1_ind(k)).label data.sensor_info(grad2_ind(k)).label];
      sens(nchan+k).typestring    = 'RSS';
      sens(nchan+k).lognum        = min([data.sensor_info(grad1_ind(k)).lognum data.sensor_info(grad2_ind(k)).lognum]) - 1;
    end
    data.sensor_info = sens;
    data.num_sensors = length(sens);
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end
  toilim        = [data.epochs.time(1) data.epochs.time(end)];
  parms.toilim  = toilim;
  %% Detection
  detection = parms.detectionflag; %   detection = 1;
  tmptype=unique({data.sensor_info.typestring}); 
  if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
  if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
  if parms.peakpairs_flag,  tmpstr = [tmpstr '_pairedpeaks'];       end
  tmpstr = [tmpstr '_' datestring];  
  detectionstring = sprintf('filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s',parms.bpfreq,toilim,[tmptype{:}],'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);
  % peak detection (grad1 & grad2)
  if detection
    if isempty(SensorEventFile) || ~exist(SensorEventFile,'file')
      SensorEventFile = sprintf('%s/%s_SO_init_peaks_%s.mat',outpath,subjects{subj},detectionstring);
    end
    clear tmpstr tmptype
    if exist(SensorEventFile,'file') && ~parms.overwrite
      fprintf(fid,'Loading sensor event file: %s\n',SensorEventFile);
      load(SensorEventFile); % events, peaks
      init_events = events; init_peaks = peaks; clear events peaks
    else
      fprintf(fid,'Resetting stopwatch timer for SO detection\n'); tstart = tic;
      args          = mmil_parms2args(parms);
      init_peaks    = SO_detection(data,args{:});
      if parms.peakpairs_flag
        save(strrep(SensorEventFile,'_pairedpeaks',''),'init_peaks');
        init_peaks  = select_peakpairs(init_peaks,parms.peakpairs_tau);
      end
      t             = data.epochs.time;
      init_events   = [];
      for k   = 1:length(init_peaks)
        init_events(k).label = data.sensor_info(k).label;
        init_events(k).time  = [t([init_peaks(k).pospeak init_peaks(k).negpeak])];
        init_events(k).type  = [1*ones(1,length(init_peaks(k).pospeak)) 2*ones(1,length(init_peaks(k).negpeak))];
      end
      fprintf(fid,'saving %s\n',SensorEventFile);
      peaks  = init_peaks;
      events = init_events;
      save(SensorEventFile,'events','peaks');
      clear peaks events
    end
    telapsed = toc(tstart); fprintf(fid,'time elapsed: %g min\n',telapsed/60);
  end
  
  bpfreq      = parms.bpfreq; if parms.bpfilter, bpfilter = 'yes'; else bpfilter = 'no'; end
  refavgflag  = parms.refavgflag;           %   refavgflag     = 1;
  plotflag    = parms.plotflag;             %   plotflag       = 1;
%   flipflag    = parms.flipflag;             %   flipflag       = 1;
%   noflipflag  = parms.noflipflag;           %   noflipflag     = 1;
  clusterflag = parms.clusterflag;
  flippad     = parms.FlipWindow;
  clusterpad  = parms.ClusterWindow;    
  epochpad    = parms.EpochPad*1000;        % s=>ms
  corrdetflag = parms.corrdetectionflag;    % whether to re-run detection after flipping
  
  if parms.preprocflag
    if ~RSSflag
      init_dat  = ts_data_selection(data,'chantype','grad');
    else
      init_dat  = data;
    end  

    fprintf(fid,'Preprocessing raw data\n');
    init_dat    = ts_preproc(init_dat,'bpfilter',bpfilter,'bpfreq',bpfreq,'blc',parms.blc,'bandpass_detrend_flag',parms.detrend);
    Fs          = init_dat.sfreq;                  % Hz
    t           = init_dat.epochs.time;            % sec
    [sel1 sel2] = match_str({init_peaks.label},{init_dat.sensor_info.label});
    init_peaks  = init_peaks(sel1);
    nchan       = init_dat.num_sensors;
    ntime       = length(init_dat.epochs.time);
    refpeaktype = parms.peaktype;               % ref peak type to use      
    if iscell(refpeaktype), refpeaktype = refpeaktype{1}; end
    npeak       = arrayfun(@(x)length(x.(refpeaktype)),init_peaks);
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  else
    Fs          = data.sfreq;
    t           = data.epochs.time;
    nchan       = data.num_sensors;
    ntime       = length(t);
  end
  % Get flip matrix
  if parms.calc_flip_matrix
    if isempty(FlipFile) || ~exist(FlipFile,'file')
      FlipFile = sprintf('%s/%s_SO_flip_matrix_%s.mat',outpath,subjects{subj},detectionstring);
    end
    if exist(FlipFile,'file') && ~parms.overwrite
      fprintf(fid,'Loading flip matrix: %s\n',FlipFile);
      load(FlipFile); % flip
    else
      flip = calc_flip_matrix(parms,init_dat,init_peaks);
      fprintf(fid,'Saving flip file: %s\n',FlipFile);
      save(FlipFile,'flip');
      telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    end
    % inspect flip matrix
    type = 1; tmpbin=bsxfun(@eq,flip(type).matrix,flip(type).matrix');
    fprintf(fid,'Percent symmetric flip matrix (%s): %g\n',flip(type).peaktype,sum(tmpbin(:))/numel(tmpbin));
    type = 2; tmpbin=bsxfun(@eq,flip(type).matrix,flip(type).matrix');
    fprintf(fid,'Percent symmetric flip matrix (%s): %g\n',flip(type).peaktype,sum(tmpbin(:))/numel(tmpbin));
  end
  
%   % Select a reference sensor
%   REF   = parms.Ref;
%   if ischar(REF)
%     if strcmp(REF,'max')                  % channel w/ the max # of detections
%       REF = find(max(npeak)==npeak);
%     elseif strcmp(REF,'all')              % all channels
%       REF = 1:init_dat.num_sensors;
%     end
%   elseif iscell(REF)                      % cell array of channel labels
%     [REF,jnk] = match_str({init_peaks.label},REF);
%   elseif isnumeric(REF)
%     % already indices
%   else
%     fprintf(fid,'reference not understood\n');
%     return
%   end
telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);

if 0 % disable flip matrix based elimination of inconsistent channels
  figure
  subplot(2,2,1); imagesc(flip(1).matrix); title('pospeak flip'); axis square
  subplot(2,2,2); imagesc(flip(2).matrix); title('negpeak flip'); axis square
  x = flip(1).matrix(:); 
  y = flip(2).matrix;
  % x(x~=y(:)) = -1;
  ix1=x~=y(:); ix2=x==y(:); x(ix1)=-1; x(ix2)=1;
  X = reshape(x,size(y));
  subplot(2,2,3); imagesc(X); title('peaks match (then 1)'); axis square
  ix1=x~=y(:); ix2=x==y(:); x(ix1)=-1; x(ix2)=1;
  ix1=X~=X'; ix2=X==X'; X(ix1 | X~=1)=-1; X(ix2 & X==1)=1;
  % X(X~=X')    = -1;
  subplot(2,2,4); imagesc(X); title('peaks & transpose match (then 1)'); axis square

  k=1; tmpbin=bsxfun(@eq,flip(k).matrix,flip(k).matrix'); sum(tmpbin(:))/numel(tmpbin)
  k=2; tmpbin=bsxfun(@eq,flip(k).matrix,flip(k).matrix'); sum(tmpbin(:))/numel(tmpbin)
  tmpbin=bsxfun(@eq,X,X'); sum(tmpbin(:))/numel(tmpbin)


  %
  x   = flip(1).matrix;
  D   = diag(x);
  n   = length(D);
  per = zeros(n,1);
  for ch = 1:n
    LHS = D==x(:,ch);
    RHS = D(ch)==x(ch,:);
    eqs = LHS==RHS';
    per(ch) = 100*sum(eqs)/length(D);
  end

  th        = 80;
  bchans    = match_str({flip(1).sensor_info.label},{flip(1).sensor_info(per<th).label});
  figure;   highlight = bchans; dewar

  % eliminate inconsistent channels
  %  based on consistency of relative polarity b/w refchan-locked averages
  th  = 90; % threshold
  clear badchans x D n per LHS RHS eqs badgrad1 badgrad2 badchanlabels pktype ch
  for pktype = 1:2 % pktype = 2;
    x   = flip(pktype).matrix;
    D   = diag(x);
    n   = length(D);
    per = zeros(n,1);
    for ch = 1:n
      LHS = D==x(:,ch);
      RHS = D(ch)==x(ch,:);
      eqs = LHS==RHS';
      per(ch) = 100*sum(eqs)/length(D);
    end
    badchans{pktype}  = match_str({flip(pktype).sensor_info.label},{flip(pktype).sensor_info(per<th).label});
    sens = flip(pktype).sensor_info;
    figure; highlight = badchans{pktype}; dewar; view(0,90); title(sprintf('%s: threshold @ %g%',flip(pktype).peaktype,th))
  end
  badchans = intersect(badchans{1},badchans{2});
  badgrad1 = strmatch('grad1',{flip(pktype).sensor_info(badchans).typestring}); nbadgrad1 = length(badgrad1); sen1 = {flip(1).sensor_info(badgrad1).label};
  badgrad2 = strmatch('grad2',{flip(pktype).sensor_info(badchans).typestring}); nbadgrad2 = length(badgrad2); sen2 = {flip(1).sensor_info(badgrad2).label};
  fprintf('Number of bad grad1 sensors: %g (%s)\nNumber of bad grad2 sensors: %g (%s)\n',nbadgrad1,[sen1{:}],nbadgrad2,[sen2{:}]);

  badchanlabels = {flip(1).sensor_info(badchans).label};

  % remove inconsistent channels from data
  [sel1,sel2] = match_str({data.sensor_info.label},badchanlabels);
  data = ts_data_selection(data,'badchans',sel1);

  % remove inconsistent channels from peaks
  if parms.peakpairs_flag
    OutFile = strrep(SensorEventPhaseFile,datestring,['FlipmatConsistencyBased-ChanRejection_' datestring]);
  else
    peaks   = init_peaks;
    OutFile = strrep(SensorEventFile,datestring,['FlipmatConsistencyBased-ChanRejection_' datestring]);
  end
  if exist(OutFile,'file') && ~parms.overwrite
    fprintf(fid,'Loading file: %s\n',OutFile);
    load(OutFile); % events, peaks
  else
    [sel1,sel2] = match_str({peaks.label},badchanlabels);
    peaks(sel1) = [];
    events      = [];
    for k       = 1:length(peaks)
      events(k).label = peaks(k).label;
      events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
      events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
    end
    fprintf(fid,'Saving file: %s\n',OutFile);
    save(OutFile,'events','peaks');
  end

  % make sure peaks & data have the same channels
  [sel1,sel2] = match_str({data.sensor_info.label},{peaks.label});
  data  = ts_data_selection(data,'channels',sel1);
  peaks = peaks(sel2);

  if ~isequal({data.sensor_info.label},{peaks.label}), error('data/peak mismatch!'); end
else
  OutFile = SensorEventFile;
end

%% define SO clusters + involved channels
if ~exist('peaks','var'), peaks = init_peaks; end
method            = 'histogram';
thresh            = 'meanstd';    % threshold for defining histogram peaks
StepSize          = .01;          % step size in sliding aggregate count
IntegrationWindow = .025;         % size of sliding window in aggregate count
MinSeparation     = .025;         % combine peaks closer than this
ClusterWindow     = .3;           % cluster window size
peaktype          = 'both';       % peaktype to use for aggregate count
OutFile = strrep(OutFile,'init_peaks',sprintf('cluster_peaks-%sPeaksBased-%s-AggStep%gs-%sThresh-ClusterWindow%gs-MinSep%gs',peaktype,method,StepSize,thresh,ClusterWindow,MinSeparation));
OutFile = strrep(OutFile,'FlipmatConsistencyBased-','');
if exist(OutFile,'file') && ~parms.overwrite
  fprintf('Loading file: %s\n',OutFile);
  load(OutFile); % events, peaks, count
else
  fprintf(fid,'Resetting stopwatch timer for aggregate detection count & SO cluster definition\n'); tstart = tic;
  [cluster_peaks,count] = find_peak_clusters(peaks,'method',method,'thresh',thresh,...
    'StepSize',StepSize,'IntegrationWindow',IntegrationWindow,'MinSeparation',MinSeparation,...
    'ClusterWindow',ClusterWindow);
  telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  t        = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
  peaks    = cluster_peaks;
  events   = [];
  for k   = 1:length(peaks)
    events(k).label = peaks(k).label;
    events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
    events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
  end
  fprintf(fid,'Saving file: %s\n',OutFile);
  save(OutFile,'events','peaks','count');
end

% find # of channels involved in each cluster
typestr   = ''; % {'','pospeak_','negpeak_'} where ''=>both
cnumfield = sprintf('%scluster_number',typestr);
cindfield = sprintf('%scluster_time_index',typestr);
clusterID = arrayfun(@(x)(x.(cnumfield)),peaks,'UniformOutput',false);
tmpID     = [clusterID{:}];
clusterID = unique(tmpID); % sorted cluster index, 1:nclusters
clusterIX = arrayfun(@(x)(x.(cindfield)),peaks,'UniformOutput',false);
clusterIX = [clusterIX{:}];
clusterIX = cellfun(@(x)unique(clusterIX(tmpID==x)),num2cell(unique(tmpID)));
nclusters = length(clusterID);
cnum      = {peaks.(cnumfield)};
cnum      = [cnum{:}];
Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));

% require that 25% of the good channels are present in each cluster
Nthresh   = round(.25*length(peaks));
ckeep     = find(Ninvolved > Nthresh);
% eliminate clusters with (# involved channels) < threshold
Ninvolved = Ninvolved(ckeep);
clusterID = clusterID(ckeep);
clusterIX = clusterIX(ckeep);
  % clusterID now contains only the clusters to process

% create cell array listing channels involved in each cluster
cnum          = {peaks.(cnumfield)};
tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
tmp           = [tmp{:}]';
tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);

% check that the number of chans in each cluster is correct:
% all(bsxfun(@eq,Ninvolved,cellfun(@length,InvolvedChans)))

% center times for each cluster
tc          = t(clusterIX);
nn          = round(length(tc)/10);
[N,X]       = hist(tc,nn);
% [cidx,ctrs] = kmeans(N,2); ix = find(ctrs==max(ctrs));

% Select time intervals of interest from detection histogram
MinTime = 60;                                 % sec, min time above thresh
Tpatch  = 45;                                 % sec, combine peaks if sep by less than Tpatch
thresh  = mean(N);                            % detection count threshold
tmp     = N > thresh;                         % 1 if count > thresh; else 0
L = find(tmp(2:end)==1 & tmp(1:end-1)==0);    % index to left edge
R = find(tmp(1:end-1)==1 & tmp(2:end)==0);    % index to right edge
if L(1)-R(1)==1                               % correct for right edge at beginning
  L = L(1:end-1);
  R = R(2:end);
end
% combine threshold crossings if right edge - left edge < MinLength=45sec
if length(L) > 2
  I       = find(X(L(2:end)) - X(R(1:end-1)) < Tpatch);
  L(I+1)  = [];
  R(I)    = [];
end
% require that the count exceed thresh for > 1 minute
D   = X(R) - X(L);
I   = find(D > MinTime);
L   = L(I);
R   = R(I);
% get start and stop times for each interval
t0  = X(L);
tf  = X(R);
t0f = [t0 tf];

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %  POTENTIAL FIGURE FOR PAPER
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Detection Count'); 
tmpx  = zeros(1,length(t));
tmpx(clusterIX) = 1;
subplot(3,1,1),plot(t,tmpx,'.-'); axis tight; title('Aggregate slow oscillation detection count')
subplot(3,1,2),try hist(tc,nn); end; axis tight
subplot(3,1,3),plot(X,N,'-'); axis tight;%,X(cidx==ix),N(cidx==ix),'b.'); axis tight
hline(thresh,'r');xlabel('time (sec)'); ylabel('count (40sec bins)');
for k=1:length(t0f), vline(t0f(k),'k'); end
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms.t0 = t0;
parms.tf = tf;
% total_toilim = [parms.t0(1) parms.tf(end)];

%% select data & peaks within intervals of interest

% interval limits
t0        = parms.t0;
tf        = parms.tf;
nInterval = length(t0);
t         = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
nchan     = length(peaks);
newpeaks  = peaks;
[newpeaks(1:nchan).pospeak]                     = deal([]);
[newpeaks(1:nchan).negpeak]                     = deal([]);
[newpeaks(1:nchan).pospeak_cluster_time_index]  = deal([]);
[newpeaks(1:nchan).pospeak_cluster_number]      = deal([]);
[newpeaks(1:nchan).negpeak_cluster_time_index]  = deal([]);
[newpeaks(1:nchan).negpeak_cluster_number]      = deal([]);
[newpeaks(1:nchan).cluster_time_index]          = deal([]);
[newpeaks(1:nchan).cluster_number]              = deal([]);
% loop over limits and constrain peaks
for k = 1:nInterval
  % PEAKS
  % convert time limits to indices
  L   = nearest(t,t0(k));
  R   = nearest(t,tf(k));
  fld = 'pospeak';   tmp = {peaks.(fld)}; tmp = cellfun(@(x)x(x>=L&x<=R),tmp,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});
  fld = 'negpeak';   tmp = {peaks.(fld)}; tmp = cellfun(@(x)x(x>=L&x<=R),tmp,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:}); 
  fld = 'pospeak_cluster_time_index';      keep = cellfun(@(x)(x>=L&x<=R),{peaks.(fld)},'uniformoutput',false);
  tmp = cellfun(@(x,y)x(y),{peaks.(fld)},  keep,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});
  fld = 'pospeak_cluster_number'; 
  tmp = cellfun(@(x,y)(x(y)),{peaks.(fld)},keep,'UniformOutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});
  fld = 'negpeak_cluster_time_index';      keep = cellfun(@(x)(x>=L&x<=R),{peaks.(fld)},'uniformoutput',false);
  tmp = cellfun(@(x,y)x(y),{peaks.(fld)},  keep,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});
  fld = 'negpeak_cluster_number'; 
  tmp = cellfun(@(x,y)(x(y)),{peaks.(fld)},keep,'UniformOutput',false);
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});
  fld = 'cluster_time_index';      keep = cellfun(@(x)(x>=L&x<=R),{peaks.(fld)},'uniformoutput',false);
  tmp = cellfun(@(x,y)x(y),{peaks.(fld)},  keep,'uniformoutput',false); 
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});
  fld = 'cluster_number'; 
  tmp = cellfun(@(x,y)(x(y)),{peaks.(fld)},keep,'UniformOutput',false);
  tmp = cellfun(@(x,y)[x y],tmp,{newpeaks.(fld)},'uniformoutput',false); [newpeaks(1:nchan).(fld)] = deal(tmp{:});  
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
orig_peaks      = peaks;
peaks           = newpeaks; clear newpeaks
BegPoints       = [1 seldat.epochs.IntervalEndPoints(1:end-1)+1];
EndPoints       = seldat.epochs.IntervalEndPoints;
IntervalT0      = seldat.epochs.time(BegPoints);
IntervalTf      = seldat.epochs.time(EndPoints);
IntervalLengths = IntervalTf - IntervalT0;

% Add code to SAVE PEAKS!!
% ......


% check that these lengths are correct
% isequal(round(diff([t0' tf'],[],2)),round(IntervalLengths)')

% % visualize selected data and peaks
% events  = [];
% for k   = 1:length(peaks)
%   events(k).label = peaks(k).label;
%   events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
%   events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
% end
% save('tmppeaks.mat','events','peaks');
% visualizer(seldat); % set xlim = [0 6000]

% reprocess cluster peaks after selecting time intervals
typestr   = ''; % {'','pospeak_','negpeak_'} where ''=>both
cnumfield = sprintf('%scluster_number',typestr);
cindfield = sprintf('%scluster_time_index',typestr);
clusterID = arrayfun(@(x)(x.(cnumfield)),peaks,'UniformOutput',false);
tmpID     = [clusterID{:}];
clusterID = unique(tmpID); % sorted cluster index, 1:nclusters
clusterIX = arrayfun(@(x)(x.(cindfield)),peaks,'UniformOutput',false);
clusterIX = [clusterIX{:}];
clusterIX = cellfun(@(x)unique(clusterIX(tmpID==x)),num2cell(unique(tmpID)));
nclusters = length(clusterID);
cnum      = {peaks.(cnumfield)};
cnum      = [cnum{:}];
Ninvolved = cellfun(@(x)sum(cnum==x),num2cell(clusterID));
% require that 25% of the good channels are present in each cluster
Nthresh   = round(.25*length(peaks));
ckeep     = find(Ninvolved > Nthresh);
% eliminate clusters with (# involved channels) < threshold
Ninvolved = Ninvolved(ckeep);
clusterID = clusterID(ckeep);
clusterIX = clusterIX(ckeep);
  % clusterID now contains only the clusters to process
% create cell array listing channels involved in each cluster
cnum          = {peaks.(cnumfield)};
tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
tmp           = [tmp{:}]';
tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);
InvolvedChans_orig = InvolvedChans;
% center times for each cluster
tc          = t(clusterIX);
nn          = round(length(tc)/10);
[N,X]       = hist(tc,nn);


% Define refchan as the channel with a peak closest to each tc
if strcmp(peaktype,'both') || length(peaktype)==2
  pks = cellfun(@(x,y)[x,y],{peaks.pospeak},{peaks.negpeak},'uniformoutput',false);
elseif ischar(peaktype)
  pks = {peaks.(peaktype)};
elseif iscell(peaktype)
  pks = {peaks.(peaktype{1})};
end

% For each cluster, only consider pks in InvolvedChans{thiscluster}
% combine all chan peaks for each cluster
% t      = seldat.epochs.time; % can't use peak time from tstart to tstop b/c of discontinuities in time
t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
    % time vector for indices stored in peaks
invpks  = cellfun(@(x)pks(x),InvolvedChans,'UniformOutput',false);
allpks = cellfun(@(x)unique([x{:}]),invpks,'UniformOutput',false);
    % |allpks| = nclusters
    %  arrays of peak indices combined across all chans involved in each cluster
% find peak closest to each tc for each array in allpks
refpks  = cellfun(@(x,y)nearest(t(x),y),allpks,num2cell(tc));
    %  index into each [allpks array of combined peak indices] for peak
    %  closest to each cluster ref time tc. |refpks| = |tc| = nclusters.
refpks  = cellfun(@(x,y)x(y),allpks,num2cell(refpks));
    %  indices for peaks closest to ref times
    %   note: these indices are present for some chan(s), but which chans
    %   was lost when the unique members of the combined sets was taken).
tref    = t(refpks); % channel detection times closest to histogram-based ref
    % the channel detection times closest to the histogram ref times
    
% % check time-interval constrained detection count histogram
% figure; try hist(t(refpks),round(length(refpks)/10)); end

% For each reference time, find the chan(s) with that detection time:
tpks            = cellfun(@(x)(t(x)),pks,'uniformoutput',false);                                % map indices to times
    % pks  contains arrays of peak indices for each channel
    % tpks contains arrays of the corresponding peak times for each channel
[tmptref{1:nchan}] = deal(tref);
    % copy of the list of channel detection times nearest to ref for each chan
ChanRefIndex    = cellfun(@(x,y)ismember(x,y)',tmptref,tpks,'uniformoutput',false);   % for which clusters each chan is the ref
    % logical arrays that indicate whether each nearest detection time
    % belongs to a particular channel; if the array for chan j has a one
    % for cluster k, then j is a refchan of k.
% convert ChanRefIndex into a list of reference channels per cluster
tmp             = [ChanRefIndex{:}]';
    % create matrix: nchan x ncluster
tmp             = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
    % group columns into cell array; one per cluster
RefChans        = cellfun(@(x)find(x),tmp,'uniformoutput',false);
    % find 1's in each column => indices to chans w peaks nearest ref times
    % note: multiple chans may have simultaneous peaks nearest ref times
SingleRefChans  = cellfun(@(x)x(1),RefChans);
    % select only one arbitrary peak nearest the ref times (smallest index)
    % = 1st chan
  
% Note: For each cluster ref time tc(k), SingleRefChans(k) is the index to
% a channel with the closest peak detection.

%% Visualization, verification, and trial rejection
% Overlay involved channels around reference times and highlight the
% reference sensor.
procdat = ts_preproc(seldat,'bpfilter','yes','bpfreq',[.1 4],'bandpass_detrend_flag',0,'blc','yes','blcwindow',[-inf inf]);
T       = procdat.epochs.time;
sens    = procdat.sensor_info;
Fs      = procdat.sfreq;
showtrl = 1:length(tc);%([1 11 44 53 64 77 81 93 101]);% showtrl = [37 39 40 44:47 55 65:67 76 81 91:97 99 101]; % []
plot_flag = 0;
normalize = 0;
winsize   = 2 ;                   % size of window to display centered at ref time
pad       = round(winsize*Fs/2);  % number of indices to display around tk
if ~isempty(showtrl)
  meanflipwave = zeros(2*pad+1,length(showtrl));
  sumflipwave  = zeros(1,length(showtrl)); 
  AIC          = zeros(length(showtrl),2);  
end

% what are the involved channels (code repeated from above)
cnum          = {peaks.(cnumfield)};
tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
tmp           = [tmp{:}]';
tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);

% tic
if plot_flag
  for c = 1:length(tc)
    if ~isempty(showtrl)
      if ~ismember(c,showtrl), continue; elseif c>max(showtrl), break; end    
    end
    % t       = procdat.epochs.time;  % time vector
    t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
    cpad    = ClusterWindow/2;      % 1/2 window used for clustering
    chans   = InvolvedChans{c};     % indices to involved channels
    refchan = SingleRefChans(c);    % index to reference channel in all data
    ref     = find(chans==refchan); % index to 1st refchan in refchan array
    tk      = tc(c);                % histogram-based reference time
    % peak times/amplitudes for each involved channel
    s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
    s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
    ts = zeros(1,length(chans));
    ix = zeros(1,length(chans));
    ix(ix1) = s1(ix1);
    ix(ix2) = s2(ix2);
    if isempty(find(ix>0,1)), continue; end
      % peak index for each chan nearest to tk in the peaks time vector
    ix = cellfun(@(x)nearest(T,t(x)),num2cell(ix)); % convert to procdat time indices
      % peak index for each chan nearest to tk in the procdat time vector
    ts = cellfun(@(x)T(x),num2cell(ix));
      % peak time for each chan from the procdat time vector
    xs = cellfun(@(x,y)procdat.epochs.data(x,y),num2cell(chans),num2cell(ix)');
    % NOTE: chan k has peak at (ts(k),xs(k))
    ind = nearest(T,tk);
    sel = ind + [-pad:pad];
    if sel(1) < 0 || sel(end) > length(T), continue; end
    tt  = T(sel);
    xx  = procdat.epochs.data(chans,sel); % data for involved chans in window
    if normalize
      xmax = max(abs(xx),[],2);
      xs   = xs ./ xmax;
      xmax = repmat(xmax,1,length(sel));
      xx   = xx ./ xmax;
    end
    % all refchans
    allreflabels = {peaks(RefChans{c}).label};
    ylim         = [-1E-8 1E-8];
    xlim         = [tk-winsize/2 tk+winsize/2];
    % PLOTS
    % overlay
    reglinewidth  = .5;
    reflinewidth  = 4;
    winlinewidth  = 3;
    peaklinewidth = 0;
    tmpcat = cat(1,xx(xs>=0,:),-xx(xs<0,:));
    tmpavg = mean(tmpcat,1);
    if plot_flag    
      figure(1); set(gcf,'Name','Before trial rejection'); subplot(3,2,1),cla
      plot(tt,xx,'-','LineWidth',reglinewidth); hold on; axis tight; set(gca,'ylim',ylim); set(gca,'xlim',xlim); xlabel('time (s)'); 
      if peaklinewidth~=0, for k=1:length(ts),vline(ts(k),'Color','y','linewidth',peaklinewidth); end; end  %   plot(ts,xs,'.','markersize',12,'Color','k','linewidth',2); 
      plot(tt,xx(ref,:),'k.-','LineWidth',reflinewidth);
      vline(tk-cpad,'Color','k','linewidth',winlinewidth); vline(tk+cpad,'Color','k','linewidth',winlinewidth);
      vline(tk,'Color','b','linewidth',winlinewidth); hline(0,'k');
      plot(ts(ref),xs(ref),'o','markersize',16,'linewidth',reflinewidth,'Color','k');
      title(sprintf('cluster %g (%gs), ref=[%s]',c,tk,[allreflabels{:}])); 
      % flip
      subplot(3,2,2);plot(tt,xx(xs>0,:),'b-',tt,-xx(xs<0,:),'r-'); set(gca,'linewidth',reglinewidth); title('flipped (red)'); hold on
      plot(tt,tmpavg,'k','LineWidth',reflinewidth); set(gca,'ylim',ylim); set(gca,'xlim',xlim);
      plot(ts(xs<0),-xs(xs<0),'.',ts(xs>=0),xs(xs>=0),'.','markersize',16,'Color','k','linewidth',reglinewidth); 
      hline(0,'k'); vline(tk,'k'); vline(tk+cpad,'k'); vline(tk-cpad,'k'); xlabel('time (s)'); hold off; 
      % dewar
      subplot(3,2,[3 5]),cla,highlight=chans(xs<0); hlcolor='r'; dewar; hold on
      highlight=chans(xs>=0); hlcolor='b'; dewar; view(0,90); 
      highlight=chans(ref); dewar; set(htext,'FontSize',12,'Color','k'); view(0,90); hold off; title('r-flip, b-no');
      % histograms
      subplot(3,2,4), try hist(ts,20); end; set(gca,'ylim',[0 15]); xlabel('detection time'); ylabel('count'); set(gca,'xlim',[tk-cpad tk+cpad]); vline(tk,'k');
      subplot(3,2,6), try hist(xs,20); end; set(gca,'ylim',[0 15]); xlabel('amplitude'); ylabel('count'); set(gca,'xlim',ylim)
      pause%(.25)
      hold off
    end
    if ~isempty(showtrl)
      meanflipwave(:,showtrl==c) = tmpavg;
      sel = tt>=(tk-cpad) & tt<=(tk+cpad);
      sumflipwave(showtrl==c) = sum(sum(tmpcat(:,sel)));
  %     options = statset('MaxIter',500);
  %     try
  %       f1 = gmdistribution.fit(ts',1,'Options',options); AIC(showtrl==c,1) = f1.AIC;
  %       f2 = gmdistribution.fit(ts',2,'Options',options); AIC(showtrl==c,2) = f2.AIC;
  %     end
  %     clear f1 f2    
    end  
  end
end
% toc
% save Workspace_s2_plot_clusters_AggHistogram_ClusterWindow0.6sec_filt0.1-4Hz_toi0-8855.14_grad1grad2.mat

%% use cluster to establish how to flip
RiThreshold    = .3;
NegEffectLevel = 0;                           % zero should be standard

outfile = sprintf('%s/s2_cluster-peaks_post-rejection_chan-then-trial_Ri%g_NegEffLevel%g_%s.mat',outpath,RiThreshold,NegEffectLevel,date); % peaks, events, flipvec, flipref, pij, rij, Ri, Eik, clusterstatus, InvolvedChans
if exist(outfile,'file') && ~parms.overwrite
  fprintf(fid,'Loading file: %s\n',outfile);
  load(outfile); % events, peaks, flipvec, flipref, pij, rij, Ri, Eik, S, clusterstatus, InvolvedChans
else
  
  % create state matrix S = [Sik] where
  %   Sik = 1  if chan i peak in cluster k is > 0
  %   Sik = -1 if chan i peak in cluster k is < 0
  %   Sik = 0  if chan i is not involved in cluster k

  % Nij is the number of clusters in which i & j are co-involved

  tic
  cnum          = {peaks.(cnumfield)};
  tmp           = cellfun(@(x)ismember(clusterID,x)',cnum,'uniformoutput',false);
  tmp           = [tmp{:}]';
  tmp           = mat2cell(tmp,size(tmp,1),ones(1,size(tmp,2)));
  InvolvedChans = cellfun(@(x)find(x),tmp,'uniformoutput',false);
  % peaks = peaks_orig

  ncluster  = length(tc);
  S         = zeros(nchan,ncluster);
  T         = procdat.epochs.time;
  sens      = procdat.sensor_info;
  Fs        = procdat.sfreq;
  normalize = 0;
  winsize   = 2 ;                                                                 % size of window to display centered at ref time
  pad       = round(winsize*Fs/2);                                                % number of indices to display around tk
  % Loop over clusters to define the polarity state matrix [s] = nchan x nclusters
  for c = 1:ncluster
    if ~isempty(showtrl)
      if ~ismember(c,showtrl), continue; elseif c>max(showtrl), break; end    
    end
    t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
    cpad    = ClusterWindow/2;                                                    % 1/2 window used for clustering
    chans   = InvolvedChans{c};                                                   % indices to involved channels
    refchan = SingleRefChans(c);                                                  % index to reference channel in all data
    ref     = find(chans==refchan);                                               % index to 1st refchan in refchan array
    tk      = tc(c);                                                              % histogram-based reference time
    % peak times/amplitudes for each involved channel
    s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
    s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
    ts = zeros(1,length(chans));
    ix = zeros(1,length(chans));
    ix(ix1) = s1(ix1);
    ix(ix2) = s2(ix2);                                                            % peak index for each chan nearest to tk in the peaks time vector
    if isempty(find(ix>0,1)), continue; end
    ix = cellfun(@(x)nearest(T,t(x)),num2cell(ix));                               % convert to procdat time indices
                                                                                  % peak index for each chan nearest to tk in the procdat time vector
    ts = cellfun(@(x)T(x),num2cell(ix));                                          % peak time for each chan from the procdat time vector
    xs = cellfun(@(x,y)procdat.epochs.data(x,y),num2cell(chans),num2cell(ix)');   % NOTE: chan k has peak at (ts(k),xs(k))
    ind = nearest(T,tk);
    sel = ind + [-pad:pad];
    if sel(1) < 0 || sel(end) > length(T), continue; end
    tt  = T(sel);
    xx  = procdat.epochs.data(chans,sel);                                         % data for involved chans in window
    if normalize
      xmax = max(abs(xx),[],2);
      xs   = xs ./ xmax;
      xmax = repmat(xmax,1,length(sel));
      xx   = xx ./ xmax;
    end
    allreflabels = {peaks(RefChans{c}).label};                                    % all refchans
    % set state for channels involved in this cluster & increment co-count
    S(chans(xs>=0),c) =  1;
    S(chans(xs< 0),c) = -1;
  end
  toc

  % -------------------------------------------------------
  % INSERT CLUSTER REJECTION
  % reject if (1) detection time histogram is bimodal, or (2) overall noise 
  % within cluster is too greatnote: this should be done BEFORE next calcs.
  % -------------------------------------------------------

  % Loop over channels to calculate a measure rij of the consistency of the
  % relative polarity between two channels i & j and the mean relative
  % polarities pij between them.

  % [s] = nchan x nclusters
  nclusters = length(tc);
  nchan     = procdat.num_sensors;
  rij = zeros(nchan,nchan);
  Nij = zeros(nchan,nchan);
  for i = 1:nchan
    for j = 1:nchan
      sel = (S(i,:)~=0) & (S(j,:)~=0);
      tmp = S(i,sel).*S(j,sel);
      if isempty(tmp), continue; end
      Nij(i,j) = length(tmp);
      rij(i,j) = abs(sum(tmp)) ./ Nij(i,j);
      pij(i,j) = mode(tmp,2);
    end
  end
  toc
  Ri = zeros(nchan,1);
  for i = 1:nchan
    sel   = rij(i,:) ~= 0;
    tmp   = sum(Nij(i,sel));
    if tmp==0, continue; end
    Ri(i) = sum(rij(i,sel) .* Nij(i,sel)) ./ tmp;
  end
  toc
  Ri_1 = Ri;
  % -------------------------------------------------------
  % OPTION: threshold Ri to eliminate inconsistent channels at this point
  % instead of using the flip matrix before clustering.
  % -------------------------------------------------------

  badchans            = find(Ri < RiThreshold);
  % remove badchans from InvolvedChans and RefChans
  InvolvedChans       = cellfun(@(x)x(~ismember(x,badchans)),InvolvedChans,'UniformOutput',false);
  RefChans            = cellfun(@(x)x(~ismember(x,badchans)),RefChans,'UniformOutput',false);
  % make sure that every cluster still has a RefChan
  sel                 = cellfun(@isempty,RefChans);
  RefChans(sel)       = cellfun(@(x)x(1),InvolvedChans(sel),'UniformOutput',false);
  % replace badchans in SingleRefChans with an element of an array in RefChans
  sel                 = ismember(SingleRefChans,badchans);
  SingleRefChans(sel) = cellfun(@(x)x(1),RefChans(sel));
  % remove i from all clusters
  S(badchans,:)       = 0;
  % Loop over clusters and channels to calculate effect of individual 
  % clusters on individual channels
  clusterstatus = zeros(nchan,nclusters);
    % 0   = not involved in cluster
    % 1   = good cluster
    % -1  = rejected cluster
  Eik   = zeros(nchan,nclusters);
  for k = 1:nclusters
    chans = InvolvedChans{k};
    nk    = length(chans);
    keep  = ones(1,nk);
    for i = 1:nk
      ch  = chans(i); % index into clusterstatus; index out of nchan channels
      p   = pij(ch,chans)';
      r   = rij(ch,chans)';
      sik = S(ch,k)*ones(nk,1);
      sjk = S(chans,k);
      Eik(ch,k) = sum(r.*(p.*(sik.*sjk))) / nk;
      if Eik(ch,k) < NegEffectLevel
        clusterstatus(ch,k) = -1;
        keep(i)             = 0;
      else
        clusterstatus(ch,k) = 1;
      end
    end
    InvolvedChans{k} = chans(keep==1);
  end
  toc
  % recalculate S, rij, pij, Ri
  ncluster  = length(tc);
  S         = zeros(nchan,ncluster);
  T         = procdat.epochs.time;
  sens      = procdat.sensor_info;
  Fs        = procdat.sfreq;
  normalize = 0;
  winsize   = 2 ;                                                                 % size of window to display centered at ref time
  pad       = round(winsize*Fs/2);                                                % number of indices to display around tk
  % Loop over clusters to define the polarity state matrix [s] = nchan x nclusters
  for c = 1:ncluster
    if ~isempty(showtrl)
      if ~ismember(c,showtrl), continue; elseif c>max(showtrl), break; end    
    end
    t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
    cpad    = ClusterWindow/2;                                                    % 1/2 window used for clustering
    chans   = InvolvedChans{c};                                                   % indices to involved channels
    refchan = SingleRefChans(c);                                                  % index to reference channel in all data
    ref     = find(chans==refchan);                                               % index to 1st refchan in refchan array
    tk      = tc(c);                                                              % histogram-based reference time
    % peak times/amplitudes for each involved channel
    s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
    s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
    ts = zeros(1,length(chans));
    ix = zeros(1,length(chans));
    ix(ix1) = s1(ix1);
    ix(ix2) = s2(ix2);                                                            % peak index for each chan nearest to tk in the peaks time vector
    if isempty(find(ix>0,1)), continue; end
    ix = cellfun(@(x)nearest(T,t(x)),num2cell(ix));                               % convert to procdat time indices
                                                                                  % peak index for each chan nearest to tk in the procdat time vector
    ts = cellfun(@(x)T(x),num2cell(ix));                                          % peak time for each chan from the procdat time vector
    xs = cellfun(@(x,y)procdat.epochs.data(x,y),num2cell(chans),num2cell(ix)');   % NOTE: chan k has peak at (ts(k),xs(k))
    ind = nearest(T,tk);
    sel = ind + [-pad:pad];
    if sel(1) < 0 || sel(end) > length(T), continue; end
    tt  = T(sel);
    xx  = procdat.epochs.data(chans,sel);                                         % data for involved chans in window
    if normalize
      xmax = max(abs(xx),[],2);
      xs   = xs ./ xmax;
      xmax = repmat(xmax,1,length(sel));
      xx   = xx ./ xmax;
    end
    allreflabels = {peaks(RefChans{c}).label};                                    % all refchans
    % set state for channels involved in this cluster & increment co-count
    S(chans(xs>=0),c) =  1;
    S(chans(xs< 0),c) = -1;
  end
  toc
  % Loop over channels to calculate a measure rij of the consistency of the
  % relative polarity between two channels i & j and the mean relative
  % polarities pij between them.

  % [s] = nchan x nclusters
  nclusters = length(tc);
  nchan     = procdat.num_sensors;
  rij = zeros(nchan,nchan);
  Nij = zeros(nchan,nchan);
  for i = 1:nchan
    for j = 1:nchan
      sel = (S(i,:)~=0) & (S(j,:)~=0);
      tmp = S(i,sel).*S(j,sel);
      if isempty(tmp), continue; end
      Nij(i,j) = length(tmp);
      rij(i,j) = abs(sum(tmp)) ./ Nij(i,j);
      pij(i,j) = mode(tmp,2);
    end
  end
  toc
  Ri = zeros(nchan,1);
  for i = 1:nchan
    sel   = rij(i,:) ~= 0;
    tmp   = sum(Nij(i,sel));
    if tmp==0, continue; end
    Ri(i) = sum(rij(i,sel) .* Nij(i,sel)) ./ tmp;
  end
  Ri_2 = Ri;
  toc

  %% visualize
  T       = procdat.epochs.time;
  sens    = procdat.sensor_info;
  Fs      = procdat.sfreq;
  showtrl = 1:length(tc);%[1 11 44 53 64 77 81 93 101]);% showtrl = [37 39 40 44:47 55 65:67 76 81 91:97 99 101]; % []
  plot_flag = 0;
  pij_flag  = 1;
  normalize = 0;
  winsize   = 2 ;                   % size of window to display centered at ref time
  pad       = round(winsize*Fs/2);  % number of indices to display around tk
  if ~isempty(showtrl)
    meanflipwave = zeros(2*pad+1,length(showtrl));
    sumflipwave  = zeros(1,length(showtrl)); 
  end
  if plot_flag
    for c = 1:length(tc)
      if ~isempty(showtrl)
        if ~ismember(c,showtrl), continue; elseif c>max(showtrl), break; end    
      end
      % t       = procdat.epochs.time;  % time vector
      t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
      cpad    = ClusterWindow/2;      % 1/2 window used for clustering
      chans   = InvolvedChans{c};     % indices to involved channels
      if isempty(chans), continue; end
      refchan = SingleRefChans(c);    % index to reference channel in all data
      ref     = find(chans==refchan); % index to 1st refchan in refchan array
      if isempty(ref)
        refchan = chans(1);
        ref     = 1;
      end
      tk      = tc(c);                % histogram-based reference time
      % peak times/amplitudes for each involved channel
      s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
      s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
      ts = zeros(1,length(chans));
      ix = zeros(1,length(chans));
      ix(ix1) = s1(ix1);
      ix(ix2) = s2(ix2);
        % peak index for each chan nearest to tk in the peaks time vector
      if isempty(find(ix>0,1)), continue; end
      ix = cellfun(@(x)nearest(T,t(x)),num2cell(ix)); % convert to procdat time indices
        % peak index for each chan nearest to tk in the procdat time vector
      ts = cellfun(@(x)T(x),num2cell(ix));
        % peak time for each chan from the procdat time vector
      xs = cellfun(@(x,y)procdat.epochs.data(x,y),num2cell(chans),num2cell(ix)');
      % NOTE: chan k has peak at (ts(k),xs(k))
      ind = nearest(T,tk);
      sel = ind + [-pad:pad];
      if sel(1) < 0 || sel(end) > length(T), continue; end
      tt  = T(sel);
      xx  = procdat.epochs.data(chans,sel); % data for involved chans in window
      if normalize
        xmax = max(abs(xx),[],2);
        xs   = xs ./ xmax;
        xmax = repmat(xmax,1,length(sel));
        xx   = xx ./ xmax;
      end
      % all refchans
      allreflabels = {peaks(RefChans{c}).label};
      ylim         = [-1E-8 1E-8];
      xlim         = [tk-winsize/2 tk+winsize/2];
      % PLOTS
      % overlay
      if pij_flag
        % use pij for flipping
        disp('using mean relative polarity vector for flipping')
        if xs(ref) < 0, flipcorr = -1; else flipcorr = 1; end
        flipvec = pij(chans(ref),chans)*flipcorr;
        noflip  = (flipvec ==  1);
        flip    = (flipvec == -1);
      else
        % use direct relative polarity for flipping
        disp('using directly calculated relative polarity for flipping')
        noflip = xs>=0;
        flip   = xs<0;
      end
      reglinewidth  = .5;
      reflinewidth  = 4;
      winlinewidth  = 3;
      peaklinewidth = 0;  
      tmpcat = cat(1,xx(noflip,:),-xx(flip,:));
      tmpavg = mean(tmpcat,1);
      if plot_flag    
        figure(1); set(gcf,'Name','After trial rejection');  subplot(3,2,1),cla
        plot(tt,xx,'-','LineWidth',reglinewidth); hold on; axis tight; set(gca,'ylim',ylim); set(gca,'xlim',xlim); xlabel('time (s)'); 
        if peaklinewidth~=0, for k=1:length(ts),vline(ts(k),'Color','y','linewidth',peaklinewidth); end; end  %   plot(ts,xs,'.','markersize',12,'Color','k','linewidth',2); 
        plot(tt,xx(ref,:),'k.-','LineWidth',reflinewidth);
        vline(tk-cpad,'Color','k','linewidth',winlinewidth); vline(tk+cpad,'Color','k','linewidth',winlinewidth);
        vline(tk,'Color','b','linewidth',winlinewidth); hline(0,'k');
        plot(ts(ref),xs(ref),'o','markersize',16,'linewidth',reflinewidth,'Color','k');
        title(sprintf('cluster %g (%gs), ref=[%s]',c,tk,[allreflabels{:}])); 
        % flip
        subplot(3,2,2);plot(tt,xx(noflip,:),'b-',tt,-xx(flip,:),'r-'); set(gca,'linewidth',reglinewidth); title('flipped (red)'); hold on
        plot(tt,tmpavg,'k','LineWidth',reflinewidth); set(gca,'ylim',ylim); set(gca,'xlim',xlim);
        plot(ts(flip),-xs(flip),'.',ts(noflip),xs(noflip),'.','markersize',16,'Color','k','linewidth',reglinewidth); 
        hline(0,'k'); vline(tk,'k'); vline(tk+cpad,'k'); vline(tk-cpad,'k'); xlabel('time (s)'); hold off; 
        % dewar
        subplot(3,2,[3 5]),cla,highlight=chans(flip); hlcolor='r'; dewar; hold on
        highlight=chans(noflip); hlcolor='b'; dewar; view(0,90); 
        highlight=chans(ref); dewar; set(htext,'FontSize',12,'Color','k'); view(0,90); hold off; title('r-flip, b-no');
        % histograms
        subplot(3,2,4), try hist(ts,20); end; set(gca,'ylim',[0 15]); xlabel('detection time'); ylabel('count'); set(gca,'xlim',[tk-cpad tk+cpad]); vline(tk,'k');
        subplot(3,2,6), try hist(xs,20); end; set(gca,'ylim',[0 15]); xlabel('amplitude'); ylabel('count'); set(gca,'xlim',ylim)
        pause%(.25)
        hold off
      end
      if ~isempty(showtrl)
        meanflipwave(:,showtrl==c) = tmpavg;
        sel = tt>=(tk-cpad) & tt<=(tk+cpad);
        sumflipwave(showtrl==c) = sum(sum(tmpcat(:,sel)));
        clear f1 f2    
      end  
    end
  end
  toc
  % fraction of channels removed from each cluster
  % [sum(clusterstatus==-1,1) ./ cellfun(@length,InvolvedChans_orig)]'
  % fraction of clusters removed from each channel
  % [sum(clusterstatus==-1,2) ./ (sum(clusterstatus==-1,2) + sum(clusterstatus==1,2))]
  % [sum(clusterstatus==-1,1)' (cellfun(@length,InvolvedChans_orig)-cellfun(@length,InvolvedChans))']

  % Loop over clusters for final rejection
  toc
  perthresh = [.4 .6]; % reject cluster if (% up) is between 40% & 60%
  ncluster  = length(tc);
  polarity  = zeros(1,ncluster);
  T         = flipdat.epochs.time;
  sens      = flipdat.sensor_info;
  Fs        = flipdat.sfreq;
  rej_flag  = 0;
  normalize = 0;
  winsize   = 2 ;                                                                 % size of window to display centered at ref time
  pad       = round(winsize*Fs/2);                                                % number of indices to display around tk
  for k = 1:ncluster
    t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
    cpad    = ClusterWindow/2;                                                    % 1/2 window used for clustering
    chans   = InvolvedChans{k};                                                   % indices to involved channels
    if isempty(chans), continue; end
    refchan = SingleRefChans(k);                                                  % index to reference channel in all data
    ref     = find(chans==refchan);                                               % index to 1st refchan in refchan array
    if isempty(ref), ref = chans(1); end
    tk      = tc(k);                                                              % histogram-based reference time
    % peak times/amplitudes for each involved channel
    s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
    s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
    ts = zeros(1,length(chans));
    ix = zeros(1,length(chans));
    ix(ix1) = s1(ix1);
    ix(ix2) = s2(ix2);                                                            % peak index for each chan nearest to tk in the peaks time vector
    if isempty(find(ix>0,1)), continue; end
    ix = cellfun(@(x)nearest(T,t(x)),num2cell(ix));                               % convert to procdat time indices
                                                                                  % peak index for each chan nearest to tk in the procdat time vector
    ts = cellfun(@(x)T(x),num2cell(ix));                                          % peak time for each chan from the procdat time vector
    xs = cellfun(@(x,y)flipdat.epochs.data(x,y),num2cell(chans),num2cell(ix)');   % NOTE: chan k has peak at (ts(k),xs(k))
    ind = nearest(T,tk);
    sel = ind + [-pad:pad];
    if sel(1) < 0 || sel(end) > length(T), continue; end
    tt  = T(sel);
    xx  = flipdat.epochs.data(chans,sel);                                         % data for involved chans in window
    if normalize
      xmax = max(abs(xx),[],2);
      xs   = xs ./ xmax;
      xmax = repmat(xmax,1,length(sel));
      xx   = xx ./ xmax;
    end
    % calc the percent of involved channels with positive deflections
    per = length(find(xs>=0)) / length(xs);
    % compare the threshold
    if (per > perthresh(1)) && (per < perthresh(2))
      clusterstatus(:,k)  = -1;
      InvolvedChans{k}    = [];
    else
      % remove channels with deflections opposite of the mode
      if per > .5
        if rej_flag
          keep = xs >= 0;
        else
          keep = ones(1,length(xs))==1; 
        end
        polarity(k) =  1;
      else
        if rej_flag
          keep = xs <  0;
        else
          keep = ones(1,length(xs))==1;
        end
        polarity(k) = -1;
      end
      S(chans(~keep),k)             = 0;
      clusterstatus(chans(~keep),k) = -1;
      InvolvedChans{k}              = chans(keep); 
    end
  end
  if rej_flag
    % update Ri
    nclusters = length(tc);
    nchan     = procdat.num_sensors;
    rij = zeros(nchan,nchan);
    Nij = zeros(nchan,nchan);
    for i = 1:nchan
      for j = 1:nchan
        sel = (S(i,:)~=0) & (S(j,:)~=0);
        tmp = S(i,sel).*S(j,sel);
        if isempty(tmp), continue; end
        Nij(i,j) = length(tmp);
        rij(i,j) = abs(sum(tmp)) ./ Nij(i,j);
        pij(i,j) = mode(tmp,2);
      end
    end
    toc
    Ri = zeros(nchan,1);
    for i = 1:nchan
      sel   = rij(i,:) ~= 0;
      tmp   = sum(Nij(i,sel));
      if tmp==0, continue; end
      Ri(i) = sum(rij(i,sel) .* Nij(i,sel)) ./ tmp;
    end  
  end
  toc
  % Create new peaks & events structures
  % note: event codes - pos/neg now corresponds to UP/DOWN
  newpeaks  = rmfield(peaks,{'pospeak_cluster_time_index','pospeak_cluster_number','negpeak_cluster_time_index','negpeak_cluster_number'});
  [newpeaks(1:nchan).pospeak]             = deal([]);
  [newpeaks(1:nchan).negpeak]             = deal([]);
  [newpeaks(1:nchan).cluster_time_index]  = deal([]);
  [newpeaks(1:nchan).cluster_number]      = deal([]);

  t       = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;
  for i   = 1:nchan
    pk    = find(clusterstatus(i,:) == 1);
    ix    = cellfun(@(x)nearest(t,tc(x)),num2cell(pk));
    pos   = polarity(pk) ==  1;
    neg   = polarity(pk) == -1;
    pkind = [peaks(i).pospeak peaks(i).negpeak];
    newpeaks(i).cluster_number     = pk;
    newpeaks(i).cluster_time_index = ix;
    newpeaks(i).pospeak = cellfun(@(x)(pkind(nearest(t(pkind),t(x)))),num2cell(ix(pos)));
    newpeaks(i).negpeak = cellfun(@(x)(pkind(nearest(t(pkind),t(x)))),num2cell(ix(neg)));
  end
  orig_peaks      = peaks;
  peaks           = newpeaks; clear newpeaks
  toc
  % create events
  events  = [];
  for i   = 1:length(peaks)
    events(i).label = peaks(i).label;
    events(i).time  = [t([peaks(i).pospeak peaks(i).negpeak])];
    events(i).type  = [1*ones(1,length(peaks(i).pospeak)) 2*ones(1,length(peaks(i).negpeak))];
  end

  % save PEAKS - these peaks will be the focus of subsequent analyses
  fprintf(fid,'Saving file: %s\n',outfile);
  save(outfile,'events','peaks','flipvec','flipref','pij','rij','Ri','Eik','S','clusterstatus','InvolvedChans');
end

% flip data matrix
flipref = find(max(Ri)==Ri);
flipvec = pij(flipref,:)';
flipmat = repmat(flipvec,1,length(procdat.epochs.time));
flipdat = procdat;
flipdat.epochs.data = flipdat.epochs.data .* flipmat;
clear flipmat

% visualize flipped data matrix w/ SO clusters for selected intervals
visualizer(ts_data_selection(flipdat,'badchans',badchans)); % load interval-constrained peak clusters w/ new codes & rejection

% W=whos;
% W={W.name};
% W=setdiff(W,{'seldat','procdat','init_dat','data'});
% save('workspace_s2_awesome_success_finally.mat',W{:});

% figure
% plot(cellfun(@length,{peaks.cluster_number}),Ri,'.')
% title('The effect of N on a channel''s consistency')
% xlabel('Number of involved clusters'); ylabel('Consistency (Ri)');
% set(gca,'ylim',[0 1]); set(gca,'xlim',[0 350]);

% figure; ii = 1:nchan; 
% plot(ii,cellfun(@length,{peaks.pospeak}),'b.',ii,cellfun(@length,{peaks.negpeak}),'r.'); 
% legend('pospeak','negpeak'); xlabel('channel'); ylabel('count'); 
% title('positive & negative clusters');

% toi = [t0' tf']; 
% disp([num2str(60 * length(tc) / sum(toi(:,2)-toi(:,1))) 's = mean time bw clusters over toilims'])


%% restoring state to when things were optimistic 
% clear all
% load workspace_s2_AWESOME_SUCCESS_finally_1.mat
% load s2_cluster-peaks_post-rejection_AWESOME-SUCCESS_1.mat
%% -----------------------------------------------------------


%% analyze clusters
clear procdat seldat data invpks allpks
t = peaks(1).tstart:1/peaks(1).sfreq:peaks(1).tstop;                            % peaks time vector
T = flipdat.epochs.time;                                                        % flipdat time vector

cfg.layout =[];
cfg.layout = layout;
cfg.layout = prepare_layout(cfg); % create 2d layout from 3d positions
distance2D = dist(cfg.layout.pos(:,1:2)');
              % distance(i,j) = sqrt( (x(i)-x(j)).^2 + (y(i) - y(j)).^2 );
% Angular distance
zshift     = .08; % shifts the Z center to almost the middle of the helmet
theta3D    = ts_BetweenSensorAngles(flipdat.sensor_info,zshift);
              % theta(i,j) = real(acos(dot(Pi,Pj)/(norm(Pi)*norm(Pj)))*180/pi)
              %              where Pk = (x,y,z)-coord of k-th sensor
theta3D(isnan(theta3D)) = 0;
  % [theta] = nchan x nchan (same for distance2D)

% remove rejected chans from distance matrices
[sel1,sel2] = match_str({data.sensor_info.label},{flipdat.sensor_info.label});
distance2D  = distance2D(sel1,sel1);
theta3D     = theta3D(sel1,sel1);

% make sure distance matrices have the right number of channels
if size(theta3D,1) ~= nchan, error('size of distance matrix is incorrect.'); end

% Update info on clusters post-rejection
% create array of cluster numbers
cnum      = {peaks.cluster_number};
cnum      = unique([cnum{:}]);

ncluster      = length(cnum);         % number of clusters with non-rejected channels
tc            = tc(cnum);             % update cluster reference times
InvolvedChans = InvolvedChans(cnum);  % update list of involved channels

% split clusters into UP & DOWN (don't know which is UP, but one is)
clusters = [];
clusters(1).code = 1; % this will be pospeak
clusters(2).code = 2; % this will be negpeak

err   = zeros(1,ncluster);
for k = 1:ncluster
  chans = InvolveChans{k};
    % will not be empty b/c empty were just removed
  tk    = tc(k);
  % find indices for each involved channel
  s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
  s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
  % should be either pospeak or negpeak for all chans
  if all(s1)
    code = 1;
  elseif all(s2)
    code = 2;
  else
    err(k) = 1;
    continue;
  end

% Loop over clusters
% find channel w detection nearest to the reference time for each cluster
%   => call that the given cluster's reference channel
delay = inf(nchan,ncluster);
  % inf => not in cluster
  % < 0 => detected before the reference
  % > 0 => detected after the reference
clusters = [];
[clusters(1:ncluster).Reference]      = deal('');
[clusters(1:ncluster).ChanLabels]     = deal({});
[clusters(1:ncluster).Distance2D]     = deal([]);
[clusters(1:ncluster).Theta3D]        = deal([]);
[clusters(1:ncluster).DetectionTime]  = deal([]);
[clusters(1:ncluster).DetectionAmp]   = deal([]);
[clusters(1:ncluster).Delay]          = deal([]);
[clusters(1:ncluster).CorrCoef2D]     = deal([]);
[clusters(1:ncluster).CorrCoef3D]     = deal([]);

for k = 1:ncluster
  chans   = InvolvedChans{k};                                                   % indices to involved channels
  if isempty(chans), continue; end
  tk      = tc(k);                                                              % histogram-based reference time
  % peak times/amplitudes for each involved channel
  s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
  s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
  ts = zeros(1,length(chans));
  ix = zeros(1,length(chans));
  ix(ix1) = s1(ix1);
  ix(ix2) = s2(ix2);                                                            % peak index for each chan nearest to tk in the peaks time vector
  if isempty(find(ix>0,1)), continue; end
  refinv  = nearest(t(ix),tk);                                                  % channel index of ref in InvolvedChans{k}
  ind     = chans(refinv);                                                      % channel index of ref out of nchan
%   % note: nearest will always return a max of 1 index => length(ind) never > 1
%   if length(ind) > 1
%     tmp = flipdat.epochs.data(ind,ix(refinv));
%     ind = ind(find(tmp == max(tmp)));
%   end
  ref(k)      = ind;                                                            % channel index of ref out of nchan
  reflabel{k} = peaks(ref(k)).label;                                            % channel label of ref
  
  ix = cellfun(@(x)nearest(T,t(x)),num2cell(ix));                               % convert to flipdat time indices
                                                                                % peak index for each chan nearest to tk in the flipdat time vector
  ts = cellfun(@(x)T(x),num2cell(ix));                                          % peak time for each chan from the flipdat time vector
  xs = cellfun(@(x,y)flipdat.epochs.data(x,y),num2cell(chans),num2cell(ix)');   % NOTE: chan k has peak at (ts(k),xs(k))
  % Store basic info
  clusters(k).Reference     = reflabel{k};
  clusters(k).ChanLabels    = {peaks(chans).label};
  clusters(k).DetectionTime = ts;
  clusters(k).DetectionAmp  = xs;
  clusters(k).DetectionIndex= ix;
  % Calculate R(delay,distance)
  clusters(k).Delay         = ts - ts(refinv);
  clusters(k).Distance2D    = distance2D(ref(k),chans);
  tmp                       = corrcoef([clusters(k).Distance2D' clusters(k).Delay']);
  clusters(k).CorrCoef2D    = tmp(1,2); clear tmp
  clusters(k).Theta3D       = theta3D(ref(k),chans);
  tmp                       = corrcoef([clusters(k).Theta3D' clusters(k).Delay']);
  clusters(k).CorrCoef3D    = tmp(1,2); clear tmp
  delay(chans,k)            = ts - ts(refinv);
end

res.reflabel  = reflabel;
res.channels  = {clusters.ChanLabels};
res.theta3D   = theta3D;
res.delay     = delay; clear delay

% quick overview of the correlation coefficients
% ii = 1:length([clusters.CorrCoef2D]); 
% figure; plot(ii,[clusters.CorrCoef2D],'b.',ii,[clusters.CorrCoef3D],'r.')

% delay(delay(:,1)~=0,k) => delays of channels involved in cluster k

% single-trial PLV (gamma), coherence, cross-correlogram

% demean before calculating xcorr


%% plots

figure
ii = 1:length(Ri_1);
subplot(1,3,1),plot(ii,Ri_1,'k.'),title('no rejection','fontsize',14); axis([0 204 0 1]); xlabel('channel','fontsize',14); ylabel('consistency of relative polarities','fontsize',14);
subplot(1,3,2),plot(ii,Ri_2,'b.'),title('relative polarity-based channel & trial rejection','fontsize',14); axis([0 204 0 1]); xlabel('channel','fontsize',14);
% subplot(1,3,3),plot(ii,Ri_3,'r.'),title('rejected inconsistent trials in indiv chans if inconsistent','fontsize',14); axis([0 204 0 1]); xlabel('channel','fontsize',14);


figure; try hist(Eik(:),150); end

