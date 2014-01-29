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
  datestring = date;
  tstart     = tic;
  tstartorig = tstart;
  toilim     = parms.toilim;
  writelog   = parms.writelog;
  layout     = '/space/mdeh1/10/halgdev/projects/jsherfey/sleep/vv_planar.lay';
  outpath    = [parms.rootoutdir '/' parms.SubjectID];
  FigPath    = sprintf('%s/images',outpath);
  logfile    = sprintf('%s/%s.log',outpath,date);
   
  subjects   = {'s1','s2','s3','s4','s5','s6','s7','s8'};
  subj       = strmatch(parms.SubjectID,subjects);
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
  
  %% ICA - EKG artifact removal
  if parms.ICAflag
    icafile    = sprintf('%s/matfiles/proc_epoch_data_ICA.mat',outpath);
    if exist(icafile,'file')
      fprintf(fid,'Loading post-ICA file: %s\n',icafile);
      load(icafile);
      data = epoch_data; clear epoch_data
    else     
      fprintf(fid,'Preprocessing data before ICA\n');
      epoch_data = ts_preproc(data,   'dsfact',     3,      ...
                                      'bpfilter',   'yes',  'bpfreq',   [.1 100], 'bandpass_detrend_flag',0,...
                                      'notch_flag', 0,      ...
                                      'blc',        'yes',  'blcwindow',[-inf inf]  );  
      telapsed   = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
      matfile    = sprintf('%s/matfiles/proc_epoch_data_1.mat',outpath); 
      fprintf(fid,'Saving matfile before running ICA: %s\n',matfile);
      save(matfile,'epoch_data');
      data       = ts_manualICA(epoch_data,'maxsteps',20,'ntrial',.01);    
      fprintf(fid,'New ICA epoch_data file: %s\n',icafile);
      clear epoch_data
    end    
  end
  
  %% SO Detection
  detection = parms.detectionflag;
  tmptype=unique({data.sensor_info.typestring}); 
  if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
  if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
  if parms.peakpairs_flag,  tmpstr = [tmpstr '_pairedpeaks'];       end
  tmpstr = [tmpstr '_' datestring];  
  detectionstring = sprintf('filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s',parms.bpfreq,toilim,[tmptype{:}],'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);
  if parms.ICAflag, detectionstring = strrep(detectionstring,'filt','ICA_filt'); end
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
  fprintf(fid,'Loading file: %s\n',OutFile);
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

outfile = sprintf('%s/s2_cluster-peaks_post-rejection_chan-then-trial_Ri%g_NegEffectLevel%g_%s.mat',outpath,RiThreshold,NegEffectLevel,date); % peaks, events, flipvec, flipref, pij, rij, Ri, Eik, clusterstatus, InvolvedChans
if exist(outfile,'file') && ~parms.overwrite
  fprintf(fid,'Loading file: %s\n',outfile);
  load(outfile); % events, peaks, flipvec, flipref, pij, rij, Ri, Eik, S, clusterstatus, InvolvedChans
else
  fprintf(fid,'Resetting stopwatch timer for consistency analysis\n'); tstart = tic;  
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
% flip data matrix

  flipref = find(max(Ri)==Ri);
  flipvec = pij(flipref,:)';
  flipmat = repmat(flipvec,1,length(procdat.epochs.time));
  flipdat = procdat;
  flipdat.epochs.data = flipdat.epochs.data .* flipmat;
  flipdat1 = flipdat;
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
    badchans = find(Ri < RiThreshold);
    S(badchans,:)             = 0;
    clusterstatus(badchans,:) = 0;
    InvolvedChans       = cellfun(@(x)x(~ismember(x,badchans)),InvolvedChans,'UniformOutput',false);
    RefChans            = cellfun(@(x)x(~ismember(x,badchans)),RefChans,'UniformOutput',false);
    % make sure that every cluster still has a RefChan
    sel                 = cellfun(@isempty,RefChans);
    RefChans(sel)       = cellfun(@(x)x(1),InvolvedChans(sel),'UniformOutput',false);
    % replace badchans in SingleRefChans with an element of an array in RefChans
    sel                 = ismember(SingleRefChans,badchans);
    SingleRefChans(sel) = cellfun(@(x)x(1),RefChans(sel));      
  end
  Ri_3 = Ri;
  
  figure
  ii = 1:length(Ri_1);
  subplot(1,3,1),plot(ii,Ri_1,'k.'),title('no rejection','fontsize',14); axis([0 204 0 1]); xlabel('channel','fontsize',14); ylabel('consistency of relative polarities','fontsize',14);
  subplot(1,3,2),plot(ii,Ri_2,'b.'),title('relative polarity-based channel & trial rejection','fontsize',14); axis([0 204 0 1]); xlabel('channel','fontsize',14);
  subplot(1,3,3),plot(ii,Ri_3,'r.'),title('rejected inconsistent trials in indiv chans if inconsistent','fontsize',14); axis([0 204 0 1]); xlabel('channel','fontsize',14);
  
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
  save(outfile,'events','peaks','flipvec','flipref','pij','rij','Ri','Eik','S','clusterstatus','InvolvedChans','badchans');
  telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
end

% tmp = {InvolvedChans{:}};
% tmp = cat(1,tmp{:});
% goodchans = unique(tmp);
% % or
% badchans  = find(Ri<RiThreshold);

if ~exist('badchans','var')
  badchans = find(Ri<RiThreshold);
end

flipref = find(max(Ri)==Ri);
flipvec = pij(flipref,:)';
flipmat = repmat(flipvec,1,length(procdat.epochs.time));
flipdat = procdat;
flipdat.epochs.data = flipdat.epochs.data .* flipmat;
clear flipmat
flipdat2 = flipdat;

% flip data matrix
flipref = find(max(Ri)==Ri);
flipvec = pij(flipref,:)';
flipmat = repmat(flipvec,1,length(seldat.epochs.time));
flipdat = seldat;
flipdat.epochs.data = flipdat.epochs.data .* flipmat;
clear flipmat
flipdat3 = flipdat;

%% some plots

figure; try hist(Eik(:),150); end; set(gca,'xlim',[.01 1]);

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


%% analyze clusters
if 0
  clear procdat invpks allpks seldat data
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
  [sel1,sel2] = match_str({peaks.label},{flipdat.sensor_info.label});
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


  ClusterFile = strrep(outfile,'cluster-peaks','cluster-structure');
  if exist(ClusterFile,'file') && ~parms.overwrite
    fprintf(fid,'Loading cluster file: %s\n',ClusterFile);
    load(ClusterFile); % clusters
  else
    % split clusters into UP & DOWN (don't know which is UP, but one is)
    clusters = [];
    clusters(1).code = 1; % this will be pospeak
    clusters(2).code = 2; % this will be negpeak
    [clusters(1:2).cluster_number] = deal([]);

    % Loop over clusters
    err   = zeros(1,ncluster);
    for k = 1:ncluster
      chans = InvolvedChans{k};
        % will not be empty b/c empty were just removed
      tk    = tc(k);
      % find indices for each involved channel
      s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),peaks(chans)); ix1 = abs(t(s1)-tk)<=cpad;
      s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),peaks(chans)); ix2 = abs(t(s2)-tk)<=cpad;
      % should be either pospeak or negpeak for all chans
      if      all(ix1) && ~all(ix2)
        code = 1;
        ix   = s1(ix1);
      elseif ~all(ix1) &&  all(ix2)
        code = 2;
        ix   = s2(ix2);
      else
        err(k) = 1;
        continue;
      end
      AbsIndex = ix;
      SelIndex = cellfun(@(x)nearest(T,t(x)),num2cell(ix));                         % convert to flipdat time indices
                                                                                    % peak index for each chan nearest to tk in the flipdat time vector
      relref    = nearest(t(AbsIndex),tk);
      ref       = chans(relref);
      reflabel  = peaks(ref).label;
      reftime   = t(ref);

      evnt  = find([clusters.code]==code);                                          % UP or DOWN
      cnt   = length(clusters(evnt).cluster_number) + 1;                            % # of this type

      if cnt == 1
        clusters(evnt).sensor_info  = flipdat.sensor_info;
        clusters(evnt).matfiles     = flipdat.epochs.matfiles;
        clusters(evnt).duration     = flipdat.epochs.duration;
        clusters(evnt).abs_toilim   = [peaks(1).tstart peaks(1).tstop];
        clusters(evnt).sel_toilim   = [flipdat.epochs.time(1) flipdat.epochs.time(end)];
        clusters(evnt).sfreq        = flipdat.sfreq;
        clusters(evnt).IntervalEndPoints = flipdat.epochs.IntervalEndPoints;
      end

      ts  = cellfun(@(x)T(x),num2cell(SelIndex));                                        % peak time for each chan from the flipdat time vector
      xs  = cellfun(@(x,y)flipdat.epochs.data(x,y),num2cell(chans),num2cell(SelIndex)'); % NOTE: chan k has peak at (ts(k),xs(k))
      D   = theta3D(ref,chans);   % angular distance from reference
      tau = ts - ts(relref);      % delay wrt reference
      tmp = corrcoef([D' tau']);  % correlation coefficient b/w distance & delay wrt ref

      clusters(evnt).cluster_number(cnt)        = cnum(k);
      clusters(evnt).epochs(cnt).HistTime       = tk;
      clusters(evnt).epochs(cnt).RefChan        = reflabel;
      clusters(evnt).epochs(cnt).RefTime        = reftime;
      clusters(evnt).epochs(cnt).InvolvedChans  = {peaks(chans).label};
      clusters(evnt).epochs(cnt).DetectionTimes = ts;
      clusters(evnt).epochs(cnt).DetectionAmps  = xs;
      clusters(evnt).epochs(cnt).Theta3D        = D;
      clusters(evnt).epochs(cnt).Delays         = tau;
      clusters(evnt).epochs(cnt).CorrCoef       = tmp(1,2);
      clusters(evnt).epochs(cnt).AbsIndex       = AbsIndex;
      clusters(evnt).epochs(cnt).SelIndex       = SelIndex;
    end
    fprintf(fid,'Saving cluster file: %s\n',ClusterFile);
    save(ClusterFile,'clusters');
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
  end

  % % Plot correlation coefficients
  % R1 = [clusters(1).epochs.CorrCoef];
  % R2 = [clusters(2).epochs.CorrCoef];
  % figure; plot(1:length(R1),R1,'b*',1:length(R2),R2,'ro'); 
  % axis tight; hline(0,'k');
  % 
  % % Plot delay vs distance with correlation coefficients
  % c = 2; th = .5;
  % R = [clusters(c).epochs.CorrCoef];    
  % ind = 1:length(R); % ind = find(abs(R) > th);
  % nr  = 8; % nr  = floor(sqrt(length(ind))) + 1;
  % nc  = 8; % nc  = floor(nr);
  % cnt = 1;
  % for k = ind
  %   D   = clusters(c).epochs(k).Theta3D;
  %   tau = clusters(c).epochs(k).Delays;
  %   subplot(nr,nc,cnt),plot(D,tau,'.'); title(sprintf('R = %g',R(k)));
  %   lsline; hline(0,'k'); %xlabel('distance'); ylabel('delay');
  %   cnt = cnt + 1; axis([0 max(D) -.25 .25]);
  %   if cnt > 64
  %     pause
  %     cnt = 1;
  %     clf
  %   end
  % end

  % ICA Artifact Removal (EKG)
  epoch_data = ts_preproc(flipdat,'dsfact',     3,      ...
                                  'bpfilter',   'yes',  'bpfreq',   [.1 100], 'bandpass_detrend_flag',0,...
                                  'notch_flag', 0,      ...
                                  'blc',        'yes',  'blcwindow',[-inf inf]  );  
  
  matfile     = sprintf('%s/matfiles/proc_epoch_data_1.mat',outpath); save(matfile,'epoch_data');
  icafile     = sprintf('%s/matfiles/proc_epoch_data_ICA.mat',outpath);
  flipdatica  = ts_manualICA(epoch_data,'maxsteps',20,'ntrial',.01);

  flipdat = flipdatica; clear flipdatica

  % Analyze clusters

  % remove badchans
  TF_RiThreshold  = .33;
  badchans        = find(Ri < TF_RiThreshold);
  peaks(badchans) = [];
  
  [sel1 sel2] = match_str({flipdat.sensor_info.label},{peaks.label});
  flipdat     = ts_data_selection(flipdat,'channels',sel1);
  
  goodchans   = {flipdat.sensor_info.label};
  
%   EpochTimesA = [clusters(1).epochs.HistTime];
%   EpochTimesB = [clusters(2).epochs.HistTime];
  Fs      = flipdat.sfreq;
  bpfreq  = [10 100];
  foi     = [10:5:100];
  width   = 7;  sf = foi / width;
%   sf      = [2 2 3 5 5 5 5 5 5]; width = foi ./ sf;
  tpad    = .7;
  tcut    = .5;
  st      = width ./ foi;
  
  t       = flipdat.epochs.time;
  pad     = round(tpad*flipdat.sfreq); pad = -pad:pad;
  sel     = round(Fs*(tpad-tcut));     sel = sel:length(pad)-sel;
  nfreq   = length(foi);
  ntime   = length(sel);
  
  WaveletFile = strrep(ClusterFile,'cluster-structure',sprintf('ICA_TF-Ri%g_cluster_structure_epochs-%g-%gs',TF_RiThreshold,tcut,tcut));
  
  if exist(WaveletFile,'file') && ~parms.overwrite
    fprintf(fid,'Loading wavelet file: %s\n',WaveletFile);
    load(WaveletFile);
  else
    warning off
    fprintf(fid,'Resetting stopwatch timer for wavelet analysis\n'); tstart = tic;
    % preallocation
    for c = 1:length(clusters)
      for k = 1:length(clusters(c).epochs)
        clusters(c).epochs(k).TFR = complex(zeros(length(clusters(c).epochs(k).InvolvedChans),ntime,nfreq));
      end
    end
    for c = 1:length(clusters)
      ntrl  = length(clusters(c).epochs);    
      for k = 1:ntrl
        fprintf('Type %g of %g: cluster %g of %g (%g min)\n',c,length(clusters),k,ntrl,toc(tstart)/60);
        tk  = clusters(c).epochs(k).HistTime;
        ind = nearest(t,tk) + pad;
        [sel1,sel2] = match_str({flipdat.sensor_info.label},clusters(c).epochs(k).InvolvedChans);
        xx  = flipdat.epochs.data(sel1,ind); % chan x time
        % bandpass
        xx  = ts_freq_filt(xx',Fs,bpfreq,[0 0],'bandpass')';
        xx  = ts_matrix2epoch(xx,'continuous',1,'time',flipdat.epochs.time(ind));
        tfr = ts_freqanalysis_fieldtrip(xx,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0,'verbose',0);
        clusters(c).epochs(k).TFR = tfr.timefreq.cmplx(:,sel,:);
      end
      clusters(c).parms.RiThreshold = TF_RiThreshold;
      clusters(c).parms.BandpassFc  = bpfreq;
      clusters(c).parms.MorletFreqs = foi;
      clusters(c).parms.MorletSF    = sf;
      clusters(c).parms.MorletST    = st;
      clusters(c).parms.MorletWidth = width;
      clusters(c).parms.DataPad     = tpad;
      clusters(c).parms.DataSel     = tcut;
    end
    fprintf(fid,'Saving wavelet file: %s\n',WaveletFile);
    save(WaveletFile,'clusters');
    telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
    warning on
  end
tic  
  nchan               = length(peaks);
  tfr.sensor_info     = flipdat.sensor_info;
  tfr.num_sensors     = flipdat.num_sensors;
  tfr.timefreq        = rmfield(tfr.timefreq,{'power','cmplx'});
  tfr.timefreq.time   = tfr.timefreq.time(sel);
  tfr.timefreq.time   = tfr.timefreq.time - median(tfr.timefreq.time);
  tfr.timefreq.power  = zeros(nchan,ntime,nfreq);
  tfr(2) = tfr;
  % average over channels
  for c = 1:length(clusters)
    % create cell array fo TFR matrices for this cluster type
    alldat = arrayfun(@(x)x.TFR,clusters(c).epochs,'uniformoutput',false);
    tfr(c).timefreq.event_code = c;
    for i  = 1:nchan
      % find clusters involving this channel
      cix  = find(arrayfun(@(x)(ismember(peaks(i).label,x.InvolvedChans)),clusters(c).epochs));
      % get indices to the matrix in those clusters for this channel
      ind  = arrayfun(@(x)find(ismember(x.InvolvedChans,peaks(i).label)),clusters(c).epochs(cix));
      % get data from those matrices for this channel
      tmp = cellfun(@(x,y)x(y,:,:),alldat(cix),num2cell(ind),'uniformoutput',false);
      % concatenate data along fourth dimension
      tmp = cat(4,tmp{:});
      % calculate power from complex spectra
      tmp = abs(double(tmp)).^2;
      % average POWER along fourth dimension
      tfr(c).timefreq.power(i,:,:) = mean(tmp,4);
      trials{c}{i}                 = cix;
      clear tmp
    end
  end
  TFR = ts_combine_data(tfr);
  clear alldat
toc 

% % TF rejection
% ntime  = length(sel);
% alldat = {};
% allmax = [];
% allmu  = [];
% warning off
% for c  = 1:length(clusters)
%   ntrl  = length(clusters(c).epochs);    
%   for k = 1:ntrl
%     % convert complex spectra to power
%     tmp = abs(double(clusters(c).epochs(k).TFR)).^2;
%     % calculate power z-score
%     tmp = (tmp - repmat(mean(tmp,2),[1 ntime 1])) ./ repmat(std(tmp,0,2),[1 ntime 1]);
%     % store this result
%     alldat{end+1} = tmp;
%     allmax(end+1) = max(tmp(:));
%     allmu(end+1)  = mean(tmp(:));
%   end
% end
% warning on   
% alldat = arrayfun(@(x)[x{1}(:)],alldat,'uniformoutput',false);
% alldat = cat(1,alldat{:});
% 
% % Mean & Max
% figure;
% subplot(1,2,1),plot(allmax,'.');title('max');
% subplot(1,2,2),plot(allmu,'.');title('mean');
% 
% ZSCORE Histogram
% figure; try hist(alldat,1000); end
% title('power z-score histogram, 10-100Hz, all clusters')
% xlabel('z-score'); ylabel('count'); axis tight


% ts_ezplot(TFR,'events',[1 2],'blc',1,'baselinetype','zscore','zlim',[-10 10],'verbose',0,'layout',layout);
% TFR = ts_combine_conditions(TFR,'combinations',{'1-2'},'neweventcodes',3);
% ts_ezplot(TFR,'events',3,'blc',1,'baselinetype','zscore','zlim',[-5 5],'verbose',0,'layout',layout);
% tfz = ts_zscore(TFR,'baselinetype','zscore','blcwindow','all','verbose',0);
% tfz = SO_freqband_average(tfz,bpfreq);
% ts_ezplot(tfz,'events',[1 2],'zlim',[-3 3],'layout',layout);
% ts_ezplot(tfz,'events',3,'zlim',[-3 3],'layout',layout);

% Multiplot average z-scores and difference
tfz = ts_zscore(ts_data_selection(TFR,'events',[1 2]),'baselinetype','zscore','blcwindow','all','verbose',0);
tfz = ts_combine_conditions(tfz,'combinations',{'1-2'},'neweventcodes',3);
xlim = [-.3 .3]; saveflag = 0; prefix = 's2_sleep';
ts_ezplot(tfz,'zlim',[-8 8],'toilim',xlim,'layout',layout,'events',[1 2],'showlabels','yes','chantype','grad','cond_labels',{'pospeak','negpeak'},'title','Spectral Power (z-score)','prefix',prefix,'save',saveflag);
ts_ezplot(tfz,'zlim',[-8 8],'toilim',xlim,'layout',layout,'events',3,'showlabels','yes','chantype','grad','cond_labels',{'pos-neg'},'title','Spectral Power (z-score)','prefix',prefix,'save',saveflag);

% Calc gamma power band average
tfz = ts_zscore(ts_data_selection(TFR,'events',[1 2]),'baselinetype','zscore','blcwindow','all','verbose',0);
tfz = ts_combine_conditions(tfz,'combinations',{'1-2'},'neweventcodes',3);
tfb = SO_freqband_average(tfz,[10 100]);
ts_ezplot(tfb,'events',[1 2],'zlim',[-3 3],'layout',layout,'chantype','grad');
ts_ezplot(tfb,'events',3,'zlim',[-3 3],'layout',layout,'chantype','grad');
% smooth the difference
for k=1:tfb.num_sensors,tfb.averages(3).data(k,:)=smooth(tfb.averages(3).data(k,:),20,'lowess'); end
ts_ezplot(tfb,'events',3,'zlim',[-3 3],'layout',layout,'chantype','grad','title','difference smoothed');

% Calc diff in t<0 and t>0
tfp = tfz; tt = tfz.averages.time;
c=1; tfp.averages(c).data = mean(tfz.averages(c).data(:,tt<0),2) - mean(tfz.averages(c).data(:,tt>0),2); tfp.averages(c).time = 0;
c=2; tfp.averages(c).data = mean(tfz.averages(c).data(:,tt<0),2) - mean(tfz.averages(c).data(:,tt>0),2); tfp.averages(c).time = 0;
c=3; tfp.averages(c).data = mean(tfz.averages(c).data(:,tt<0),2) - mean(tfz.averages(c).data(:,tt>0),2); tfp.averages(c).time = 0;
ts_ezplot(tfp,'events',1,'topoplot',1,'toprows',1,'topcols',1,'title','pos','style','straight','chantype','grad2');
ts_ezplot(tfp,'events',2,'topoplot',1,'toprows',1,'topcols',1,'title','neg','style','straight','chantype','grad2');
ts_ezplot(tfp,'events',3,'topoplot',1,'toprows',1,'topcols',1,'title','pos-neg','style','straight','chantype','grad1');

% figure; 
% subplot(2,1,1),plot(tfz.averages(1).time,mean(tfz.averages(1).data,1)); vline(0,'k'); hline(0,'k');
% subplot(2,1,2),plot(tfz.averages(2).time,mean(tfz.averages(2).data,1)); vline(0,'k'); hline(0,'k');
% 
% display N
% tmp = [cellfun(@length,trials{1})' cellfun(@length,trials{2})'];
% lab = {peaks.label};
% for k = 1:length(lab),fprintf('%s:  %4g  %4g\n',lab{k},tmp(k,:)); end

% % Multiplot individual clusters
% c = 1; k = 1;
% for k = 1:min(length(clusters(1).epochs),length(clusters(2).epochs))
%   for c = 1:2
%     tmpTFR = ts_data_selection(TFR,'chanlabel',clusters(c).epochs(k).InvolvedChans);
%     tmpTFR.timefreq = tmpTFR.timefreq(1);
%     tmpTFR.timefreq.power = abs(clusters(c).epochs(k).TFR).^2;
%     ts_ezplot(tmpTFR,'showlabels','yes','title',num2str(k),'layout',layout,'chantype','grad','blc',1,'baselinetype','zscore','zlim',[-10 10],'verbose',0);
%   end
%   pause
%   close(1); close(2);
% end

if 0 
  chansel             = 1:15; % channels to analayze
  tfr                 = tfr(1);
  tfr.sensor_info     = flipdat.sensor_info(chansel);
  tfr.num_sensors     = length(chansel);
  tfr.timefreq.cmplx  = complex(zeros(length(chansel),length(flipdat.epochs.time),nfreq));
  tstart  = tic;
  for ch  = 1:length(chansel)
    k     = chansel(ch);
    fprintf('processing channel %g of %g (%g min)\n',k,nchan,toc(tstart)/60);
    tmp = ts_data_selection(flipdat,'channels',k);
    tmp = ts_preproc(tmp,'bpfilter','yes','bpfreq',bpfreq,'bandpass_detrend_flag',0);
    tmp = ts_freqanalysis_fieldtrip(tmp,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0,'verbose',0);
    tfr.timefreq.cmplx(k,:,:) = tmp.timefreq.cmplx;
    clear tmp
  end
  tfr.sfreq = flipdat.sfreq;
  tfr.timefreq.power  = abs(tfr.timefreq.cmplx).^2;
  tfr.timefreq.time   = flipdat.epochs.time;
  tfz = ts_zscore(ts_data_selection(tfr,'events',[1 2]),'baselinetype','zscore','blcwindow','all','verbose',0);
  visualizer(flipdat,tfz);%visualizer(flipdat,tfr);

  morebadchans = {'MEG 0113','MEG 0143'};
end
  % Loop over clusters
    % create raw wave data:   [Ninvolved x EpochSize]
    % create filt wave data:  [Ninvolved x EpochSize]
    % create TF data:         [Ninvolved x EpochSize x Nfreq]
      % create from raw data with data padding
      % remove padding after calculation

  % What about TF normalization? ...
  % What about single-trial PLV? ...


  % Create average raw data, filt data, and TF data ... How? ...

end

% demean before calculating xcorr


