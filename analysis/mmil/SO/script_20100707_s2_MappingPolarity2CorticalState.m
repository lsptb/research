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
  
  % Select a reference sensor
  REF   = parms.Ref;
  if ischar(REF)
    if strcmp(REF,'max')                  % channel w/ the max # of detections
      REF = find(max(npeak)==npeak);
    elseif strcmp(REF,'all')              % all channels
      REF = 1:init_dat.num_sensors;
    end
  elseif iscell(REF)                      % cell array of channel labels
    [REF,jnk] = match_str({init_peaks.label},REF);
  elseif isnumeric(REF)
    % already indices
  else
    fprintf(fid,'reference not understood\n');
    return
  end
telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);


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


%%
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



%%






%% Sitting upright at MGH & laying down at UCSD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 05-Jul-2010, JSS
% TODO: Test various methods for mapping MEG polarity to cortical state.
%
% Methods of defining SO epochs:
  % 1. Fixed window size centered on detection times
  % 2. Variable size defined by level-crossings before & after detection
  % times (potential levels: zero, mean or median of max & min over window
  % around individual peaks or over all times).
  % 3. Fixed size interpolation derived from variable level-crossing based
  % start and stop times.  
% Methods of classifying peaks:
  % 1. Normalized, event-related TFR-based mean gamma differences.
  % 2. Sliding window gamma-filtered level-crossing analysis.
  % 3. SO phase-locked averaging of smoothed, Hilbert-based gamma envelope.

%% Test #1 => Epoch method #1 + Classification method #3.
% Load data
tstart    = tic;
T         = data.epochs.time;
X         = data.epochs.data;
testchans = {data.sensor_info.label}; %{'MEG 0143','MEG 0142','MEG 0213','MEG 0212'};
PhiNodes  = linspace(-pi/2,pi/2,40);%-pi/2:pi/48:pi/2;
tlim      = toilim;%[1700 2700];%[500 2700];%[2000 2500];
sofreq    = parms.bpfreq;
% gmfreq    = [20 100];
% GammaRes  = []; % gamma will be smoothed over a window = length(epoch)/GammRes
% GammaEnv  = 0;
winsize   = 2.5;                                    % window size in seconds
% % note: window should be large enough to include peak and surrounding troughs
% select paired peaks
peaks     = init_peaks;
peaks     = select_peakpairs(peaks,1);
% select subset of peaks & data (test channels & time limits)
% script_select_data
  % T & X         => all times & select chans
  % data & peaks  => select times & select chans
  % others: skipchans (no peaks) & pkt (time vector corresponding to peak indices)
% toc
t         = data.epochs.time;
ntime     = length(t);
npeak     = cellfun(@length,{peaks.pospeak});
nchan     = data.num_sensors;
nNode     = length(PhiNodes);
Fs        = data.sfreq;
pad       = ceil(winsize*Fs/2);
peaktypes = {'pospeak','negpeak'};
debugflag = 0;
% xx=ts_freq_filt(X(ix),Fs,[.1 3],[0 0],'bandpass'); phi=angle(hilbert(xx));

tmptype=unique({data.sensor_info.typestring}); 
if parms.derivpeaks,      tmpstr = '_derivpks'; else tmpstr = ''; end
if parms.onlysinglepeaks, tmpstr = [tmpstr '_nodoublepks'];       end
tmpstr = [tmpstr '_pairedpeaks'];
tmpstr = [tmpstr '_with-SO-phase-locked-analysis'];
tmpstr = [tmpstr '_' datestring];  
detectionstring = sprintf('filt%g-%gHz_toi%g-%g_%s_ref%s_smooth%gsec_zero2zero%g-%gsec_mindist%gsec%s',parms.bpfreq,toilim,[tmptype{:}],'all',parms.smooth_window,parms.zero2zero_limits,parms.mindist,tmpstr);
% peak detection (grad1 & grad2)
if isempty(SensorEventPhaseFile) || ~exist(SensorEventPhaseFile,'file')
  SensorEventPhaseFile = sprintf('%s/%s_SO_init_peaks_%s.mat',outpath,subjects{subj},detectionstring);
end
clear tmpstr tmptype
if exist(SensorEventPhaseFile,'file') && ~parms.overwrite
  fprintf(fid,'Loading sensor event file with phase-locked results: %s\n',SensorEventPhaseFile);
  load(SensorEventPhaseFile); % events, peaks
else
  fprintf(fid,'Resetting stopwatch timer for SO phase-locked analysis\n'); tstart = tic;
  % for each channel
  for ch = 1:nchan
    fprintf('Channel %g of %g (%g min)\n',ch,nchan,toc(tstart)/60);
    label    = testchans{ch};
    for type = 1:length(peaktypes)
      pktype = peaktypes{type};
      pks    = peaks(ch).(pktype);
      % bandpass filter for slow oscillation
      x = ts_freq_filt(X(ch,:)',Fs,sofreq,[0 0],'bandpass')';
      x = x(nearest(T,t(1)):nearest(T,t(end)));
      % bandpass filter for gamma activity    
  %     g = ts_freq_filt(X(ch,:)',Fs,gmfreq,[0 0],'bandpass')';
  %     g = ts_freq_filt(g',Fs,60,3,'notch')';
  %     g = g(nearest(T,t(1)):nearest(T,t(end)));

      % if negpeak => flip x
      if strcmp(pktype,'negpeak'), flipval = -1; else flipval = 1; end
      x = flipval*x;

      % shift pks since length(t) < length(T) & pks comes from T
      pks      = pks - nearest(T,t(1)) + 1;
      % eliminate pks with insufficient data padding
      pks      = pks((pks-pad>=1)&(pks+pad<=ntime));
      % preallocation
      npk      = length(pks);
      xnodes   = zeros(nNode,npk);
      ynodes   = zeros(nNode,npk);
  %     gnodes   = zeros(nNode,npk);
      skip     = zeros(1,npk);
  %     GammaMinusYlag = zeros(1,npk);
      % loop over detection times {tk}
      for k = 1:npk
        tki = pks(k);
        tk  = t(tki);
        ii  = tki + (-pad:pad);
        xx  = x(ii);
  %       gg  = g(ii);
  %       if GammaEnv,           gg = abs(hilbert(gg)); end
  %       if ~isempty(GammaRes), gg = smooth(gg,ceil(length(ii)/GammaRes),'lowess'); end
        tt  = t(ii);
        L   = (max(xx)+min(xx))/2;
        xL  = xx - L;
        yL  = xL / xL(floor(length(ii)/2)+1); %max(xL);
        phi = angle(hilbert(xL));

        % find the last -pi/2 crossing before phi=0
        [ind,t0] = crossing(phi,tt,-pi/2); ind1 = max(ind(t0<tk));
        % find the first pi/2 crossing after phi=0
        [ind,t0] = crossing(phi,tt,pi/2);  ind2 = min(ind(t0>tk));
        % skip this peak if either is empty or there are more than 1 zero-crossings between them
        if isempty(ind1) || isempty(ind2) || length(crossing(phi(ind1:ind2)))>1
          skip(k) = 1; continue;
        end
        % find all crossings in window for each phinode
        ix  = cellfun(@(x)(crossing(phi,[],x)),num2cell(PhiNodes),'UniformOutput',false);
        % select only crossings with positive slope
        ix  = cellfun(@(x)(x(phi(x+1)-phi(x)>0)),ix,'UniformOutput',false);
        % select crossings closest to tk
        ix  = cellfun(@(x)x(nearest(tt(x),tk)),ix);

        xL = flipval*xL;
        yL = flipval*yL;

        if debugflag
          % plot to verify monotonic phase; at least verify the crossings surround tk
          figure; 
          subplot(8,1,1),plot(tt,xx,'b'); axis tight; title('so'); hline(L,'k'); vline(tk,'k');
          subplot(8,1,2),plot(tt,xL,'b',tt(ix),xL(ix),'b*'); axis tight; title('xL: so shifted');hline(0,'k');vline(tk,'k');
          subplot(8,1,3),plot(tt,yL,'g',tt(ix),yL(ix),'g*'); axis tight; title('yL: so shifted scaled');hline(0,'k');vline(tk,'k');
          subplot(8,1,4),plot(tt,phi,'k',tt(ix),phi(ix),'k*'); axis tight; title('so phase');hline(0,'k');vline(tk,'k');
  %         subplot(8,1,5),plot(tt,gg,'r',tt(ix),gg(ix),'r*'); axis tight; title('smooth gamma filt');hline(0,'k');vline(tk,'k');
          subplot(8,1,6),plot(PhiNodes,xL(ix),'b*-'); axis tight; title('xL on phi');hline(0,'k');vline(median(PhiNodes),'k');
          subplot(8,1,7),plot(PhiNodes,yL(ix),'g*-'); axis tight; title('yL on phi');hline(0,'k');vline(median(PhiNodes),'k');
  %         subplot(8,1,8),plot(PhiNodes,gg(ix),'r*-'); axis tight; title('gamma on phi');hline(0,'k');vline(median(PhiNodes),'k');
  %         figure('Name','v = predefined phi values');
  %         [tmp,lags] = xcorr(yL(ix),gg(ix),'coeff'); dphi=PhiNodes(2)-PhiNodes(1);
  %         subplot(3,1,1),plot(lags*dphi,tmp); axis tight; title('xcorr(yL(v),gamma(v))');
  %         [tmp,lags] = xcorr(phi(ix),yL(ix),'coeff');
  %         subplot(3,1,2),plot(lags*dphi,tmp); axis tight; title('xcorr(phi(v),yL(v))');
  %         tmp1=lags(max(tmp)==tmp); vline(tmp1*dphi,'k');
  %         [tmp,lags] = xcorr(phi(ix),gg(ix),'coeff');
  %         subplot(3,1,3),plot(lags*dphi,tmp); axis tight; title('xcorr(phi(v),gamma(v))');
  %         tmp2=lags(max(tmp)==tmp); vline(tmp2*dphi,'k');
  %         xlabel(sprintf('gamma\\_lag - yL\\_lag = %g radians\n',(tmp2-tmp1)*dphi));
          pause; close
        end
        % get x-values at PhiNodes
        xnodes(:,k) = xL(ix);
        ynodes(:,k) = yL(ix);
  %       gnodes(:,k) = gg(ix);

        % cross-correlations
        dphi       = PhiNodes(2)-PhiNodes(1);
  %       [tmp,lags] = xcorr(phi(ix),yL(ix),'coeff'); tmp1=min(lags(max(tmp)==tmp));
  %       [tmp,lags] = xcorr(phi(ix),gg(ix),'coeff'); tmp2=min(lags(max(tmp)==tmp));
  %       GammaMinusYlag(k) = (tmp2-tmp1)*dphi;
    %     [tmp,lags] = xcorr(yL(ix),gg(ix),'coeff');
        clear tmp tmp1 tmp2
      end
      peaks(ch).([pktype '_skip'])           = skip;
  %     peaks(ch).([pktype '_GammaMinusYlag']) = GammaMinusYlag;
      peaks(ch).PhiNodes             = PhiNodes;
      peaks(ch).([pktype '_xnodes']) = xnodes;
      peaks(ch).([pktype '_ynodes']) = ynodes;
  %     peaks(ch).([pktype '_gnodes']) = gnodes;
      % calculate tentative average (no outlier rejection)
      peaks(ch).([pktype '_xnode_avg']) = sum(xnodes(:,~skip),2)/sum(~skip);
      peaks(ch).([pktype '_ynode_avg']) = sum(ynodes(:,~skip),2)/sum(~skip);
  %     peaks(ch).([pktype '_gnode_avg']) = sum(gnodes(:,~skip),2)/sum(~skip);
      clear xnodes ynodes gnodes xx gg tt xL yL phi x g pks ix GammaMinusYlag skip
    end % end loop over peak types
  %   toc
  end
  events = init_events;
  fprintf(fid,'saving %s\n',SensorEventPhaseFile);
  save(SensorEventPhaseFile,'events','peaks');
  telapsed = toc(tstart); fprintf(fid,'Time elapsed: %g min\n',telapsed/60);
end
nkeep   = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
thresh  = median(nkeep) - .5*std(nkeep);
figure
highlight = find(nkeep>thresh); subplot(1,2,1); dewar; title('many detections');
highlight = find(nkeep<thresh); subplot(1,2,2); dewar; title('few detections');

% set outliers to zero
% get nNode x N's from # of non-zero elements along 2nd dimension
% average over remaining peaks
chans    = 1:nchan; 
% posgamma_integral = arrayfun(@(x)(sum(x.pospeak_gnode_avg-min([x.pospeak_gnode_avg;x.pospeak_gnode_avg]))),peaks);
% posgamma_integral(posgamma_integral==0) = [];
% neggamma_integral = arrayfun(@(x)(sum(x.negpeak_gnode_avg-min([x.pospeak_gnode_avg;x.pospeak_gnode_avg]))),peaks);
% neggamma_integral(neggamma_integral==0) = [];
% %   figure; plot(1:nchan,posgamma_integral,'b.',1:nchan,neggamma_integral,'ro'); title('gamma integrals');
% %   legend('pos','neg'); [posgamma_integral' neggamma_integral']
% diffgam  = posgamma_integral-neggamma_integral; 
% thresh   = median(abs(diffgam))+2*std(abs(diffgam));
% badchans = abs(diffgam)>(thresh);
%   figure('Name','Gamma Integral: pos - neg'); subplot(2,1,1),plot(chans,diffgam,'.-',chans(badchans),diffgam(badchans),'ro');
%   axis tight; hline(0,'k'); hline(-thresh,'r'); hline(thresh,'r');
%   subplot(2,1,2),plot(chans(~badchans),diffgam(~badchans),'.-'); title('without outlier channels');
%   axis tight; if ~isnan(thresh),set(gca,'ylim',[-1 1]*(1.2*thresh)); hline(0,'k'); hline(-thresh,'r'); hline(thresh,'r'); end
%   {peaks(badchans).label}

nposkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
nnegkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
[nposkeep' nnegkeep'] % equivalent
nkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
figure; subplot(1,2,1),plot(nkeep,'.'); level1 = median(nkeep); hline(level1,'r');
tmp = .5*std(nkeep); hline(level1+tmp,'g'); hline(level1-tmp,'g');
title('nkeep'); xlabel('chan'); axis tight
subplot(1,2,2),try hist(nkeep,round(length(nkeep)/15)); end;

% two subplots: xnode & ynode; each w/ a pospeak/negpeak overlay
% each curve has nNode points
thresh = median(nkeep) - .5*std(nkeep);
nn = 10; startchan = 1; highlight = find(nkeep>thresh);
if 0
  nsel     = length(highlight);
  sodiff   = zeros(nchan,2);
  MaxPos   = zeros(nchan,1); % is the max neg or pos
  MaxLeft  = zeros(nchan,1); % is the max left or right
  evcode   = zeros(nchan,1);
  areacode = zeros(nchan,1);
  evcodes  = [0 11 21 12 22];
%   [gamsumpos{1:5}] = deal([]);
%   [gamsumneg{1:5}] = deal([]);
  for p=1:max(3,ceil(nsel/nn))
    figure
    ch0 = startchan + (p-1)*nn;
    if isempty(highlight)
      thesechans = [ch0:ch0 + nn - 1]; 
    else
      thesechans = highlight(ch0:min(ch0 + nn - 1,length(highlight)));
    end
    cnt = 0;
    for ch = thesechans
      cnt = cnt + 1;
      off = 0;%min([peaks(ch).pospeak_xnode_avg;peaks(ch).pospeak_xnode_avg]);
      y11 = peaks(ch).negpeak_xnode_avg - off; y12 = peaks(ch).pospeak_xnode_avg - off;
      off = 0;%min([peaks(ch).pospeak_gnode_avg;peaks(ch).pospeak_gnode_avg]);
%       y21 = peaks(ch).negpeak_gnode_avg - off; y22 = peaks(ch).pospeak_gnode_avg - off;
      x   = peaks(ch).PhiNodes*(180/pi);
      % take absolute value of y12 so that pos & neg are both up
      y11abs = abs(y11);
      % subtract y(t=0) from y11 & y12
      y11abs = y11abs - y11abs(round((1+length(y11abs))/2));
      y12 = y12 - y12(round((1+length(y12))/2));
      tmp = y12-y11abs;
      sodiff(ch,1) = sum(tmp(x<0));
      sodiff(ch,2) = sum(tmp(x>0));
      tmpcode      = sodiff(ch,1) - sodiff(ch,2);
      if tmpcode > 1, areacode(ch,1) = -1; else areacode(ch,1) = 1; end
        % -1 => left area is larger; 1 => right area is larger
      % get max vals
%       ix           = find(abs(tmp)==max(abs(tmp)));
      ix           = find(tmp==max(tmp));
      MaxLeft(ch)  = PhiNodes(ix)<0;
      MaxPos(ch)   = tmp(ix)>0;
      evcode(ch,:) = str2num(sprintf('%g%g',[MaxLeft(ch)+1 MaxPos(ch)+1]));
      if MaxLeft(ch), str='left';       else str='right';       end
      if MaxPos(ch),  str=[str '-pos']; else str=[str '-neg'];  end
      str = [str ': ' num2str(evcode(ch))];
      subplot(4,nn,cnt),plot(x,y12,'b',x,y11,'r'); axis tight; 
      title(sprintf('%s:n=%g/%g',peaks(ch).label,sum(~peaks(ch).pospeak_skip),sum(~peaks(ch).negpeak_skip))); 
      vline(0,'k'); if cnt==1, legend('pos','neg'); end
      subplot(4,nn,nn+cnt),plot(x,y12-y11abs,'k'); axis tight; hline(0,'k'); vline(0,'k'); title(str); axis off % title('SO: pos-|neg|)');
%       subplot(4,nn,2*nn+cnt),plot(x,y22,'b',x,y21,'r'); axis tight; title('mean gamma'); vline(0,'k');
%       subplot(4,nn,3*nn+cnt),plot(x,smooth(y22-y21,5,'lowess'),'k'); axis tight; title('gamma: pos-neg'); hline(0,'k'); vline(0,'k');
%       xlabel('instantaneous phase of SO (Hz)'); % maybe smooth it
%       gamsumpos{evcodes==evcode(ch)}(end+1,:) = y22;
%       gamsumneg{evcodes==evcode(ch)}(end+1,:) = y21;
      clear y11 y12 y21 y22
    end
  end
  %   figure
  %   subplot(3,1,1),plot(sodiff(:,1),'.'); title('sum over pos-|neg| for phi<0');
  %   subplot(3,1,2),plot(sodiff(:,2),'.'); title('sum over pos-|neg| for phi>0');
  %   subplot(3,1,3),plot(sodiff(:,1)-sodiff(:,2),'.'); title('(sum for phi<0)-(sum for phi>0)');
  %   xlabel('channel');
  %   clear sodiff
  for k=1:length(evcode),fprintf('%3g [%s]\n',evcode(k),peaks(k).label); end
  poscode   = [12 22]; negcode   = [11 21];
  leftcode  = [21 22]; rightcode = [11 12];
  evcodes   = unique(evcode)';
  ev0 = {peaks(evcode==0).label}';  n0 = length(ev0);
  ev1 = {peaks(evcode==11).label}'; n1 = length(ev1);
  ev2 = {peaks(evcode==21).label}'; n2 = length(ev2);
  ev3 = {peaks(evcode==12).label}'; n3 = length(ev3);
  ev4 = {peaks(evcode==22).label}'; n4 = length(ev4);
  nev = [n0 n1 n2 n3 n4];
  disp([evcodes;nev(nev~=0)])
  % visualizer(ts_data_selection(data,'chanlabel',ev0));
  % visualizer(ts_data_selection(data,'chanlabel',ev1));
  % visualizer(ts_data_selection(data,'chanlabel',ev2));
  % visualizer(ts_data_selection(data,'chanlabel',ev3));
  % visualizer(ts_data_selection(data,'chanlabel',ev4));
  for k = 1:length(peaks)
    peaks(k).event_code = evcode(k);
  end
  poscode   = [12 22]; negcode   = [11 21];
  leftcode  = [21 22]; rightcode = [11 12];
  for k = 2:length(evcodes)
    highlight = find([peaks.event_code]==evcodes(k));
    if k == 2, figure; end
    subplot(2,2,k),dewar
    title(sprintf('Event %g',evcodes(k))); view([0 90])
  end
  figure
  subplot(1,2,1),bchan1 = find([peaks.event_code]==0); highlight=bchan1; dewar; view([0 90]);
  th = 80;  title('bad chans (slope method)');
  bchan2 = match_str({flip(1).sensor_info.label},{flip(1).sensor_info(per<th).label})';
  subplot(1,2,2),highlight = bchan2; dewar; view([0 90]);
  {peaks(bchan1(~ismember(bchan1,bchan2))).label} % bad in method 2 not 1
  {peaks(bchan2(~ismember(bchan2,bchan1))).label} % bad in method 1 not 2
  title('bad chans (mean polarity (flip) method)');
  bchans = unique([bchan1 bchan2]);  
end

%% II. Classifying polarities: use SO phase-locked averaging of smoothed, 
% % Hilbert-based gamma envelope.
% 
% % level-crossing analysis
% ngrid     = 1500;
% stepsize  = 25;          % # of points to shift each step
% npts      = 50;         % # of points in the window
% ch        = 4;
% x         = peaks(ch).PhiNodes*(180/pi);
% nx        = length(x);
% steps     = 1:stepsize:[ngrid-npts+1];
% nsteps    = length(steps);
% allcounts = {};
% 
% for type = 1:length(peaktypes)
%   pktype  = peaktypes{type};%'pospeak';
%   keep    = find(~(peaks(ch).pospeak_skip | peaks(ch).negpeak_skip));
%   pks     = peaks(ch).(pktype);
%   pks     = pks(keep);
%   npks    = length(pks);
%   count   = zeros(npks,nsteps);
%   for k = 1:npks
%     % select gamma=f(phi) around this peak
%     G        = peaks(ch).([pktype '_gnodes'])(:,keep(k)); % gamma = f(preset_phi)
%     xi       = linspace(PhiNodes(1),PhiNodes(end),ngrid);
%     Gi       = interp1(PhiNodes,G,xi,'spline');
%     L        = median(Gi)+std(Gi);%.9*max(Gi);
%     [ind,t0] = crossing(Gi,xi,L);
%     if k == 1, figure; plot(PhiNodes,G,'o',xi,Gi,'.-',xi(ind),Gi(ind),'g*'); axis tight; hline(L,'r'); title(pktype); xlabel('phi (rad)'); ylabel('gamma'); end
%     for j = 1:nsteps
%       this       = find((ind>=steps(j))&(ind<=(steps(j)+npts-1)));
%       count(k,j) = length(this);
%     end
%     clear G Gi L ind t0 this
%   end
%   allcounts{type} = count;
%   clear count pks npks
% end
% toc
% 
% tt = xi(steps+floor(stepsize/2))*(180/pi); % t(steps);
% figure; subplot(2,1,1),plot(tt,sum(allcounts{1},1),'b*-',tt,sum(allcounts{2},1),'ro-'); axis tight
% legend('pos','neg'); title('high-level crossing gamma counts'); vline(0,'k'); set(gca,'xlim',[-90 90]);
% tmp=sum(allcounts{1},1)-sum(allcounts{2},1); subplot(2,1,2)
% plot(tt,tmp,'k-'); axis tight; set(gca,'xlim',[-90 90]); vline(0,'k'); hline(0,'k');
% title('pos - neg'); ylabel('gamma difference'); xlabel('SO phase'); hold on
% [cidx,ctrs]=kmeans(tmp,2); plot(tt(cidx==1),tmp(cidx==1),'go',tt(cidx==2),tmp(cidx==2),'bo');
% 
% a = sum((tt(cidx==1)<0)); b = sum((tt(cidx==1)>0));
% c = sum((tt(cidx==2)<0)); d = sum((tt(cidx==2)>0));
% fprintf('cluster count for t<0 vs t>0:\n')
% fprintf('cluster 1 (diff @ %g): %g  %g\n',ctrs(1),a,b); fprintf('cluster 2 (diff @ %g): %g  %g\n',ctrs(2),c,d)
% [a-b c-d],[a-c b-d]

% figure
% h1 = find(diff(sodiff,[],2) > 0);  highlight=h1; subplot(2,2,1); dewar; view(0,90); title('(-) DOWN'); {flip(1).sensor_info(highlight).label}
% h2 = find(diff(sodiff,[],2) < 0);  highlight=h2; subplot(2,2,2); dewar; view(0,90); title('(+) DOWN'); {flip(1).sensor_info(highlight).label}
% h3 = find(diff(sodiff,[],2) == 0); highlight=h3; subplot(2,2,3); dewar; view(0,90); title('indeterminate'); {flip(1).sensor_info(highlight).label}
% 
% code = zeros(1,nchan);
% code(h1) = 1;
% code(h2) = -1;
% 
% % h1 [flipmat] -1
% % h2 [flipmat] +1
% intersect(h1,find(flip(1).matrix(1,:)==-1))
% intersect(h2,find(flip(1).matrix(1,:)==1))

%% quick look at what may be close to the flip-corrected matrix + events...
% 
% refchan = 1;
% flipvec = flip(1).matrix(:,refchan);
% 
% tmpdat  = init_dat;
% tmpdat.epochs.data = tmpdat.epochs.data .* repmat(flipvec,[1 length(tmpdat.epochs.time)]);
% 
% tmpevents = init_events;
% for k = 1:nchan
%   tmp = tmpevents(k).type;
%   if flipvec(k)
%     tmpevents(k).type(tmp==1) = 2;
%     tmpevents(k).type(tmp==2) = 1;
%   end
% end
% otherbadchans = [1 7]; % visually identified
% allbadchans   = [bchans otherbadchans];
% tmpevents(allbadchans)  = [];
% tmpdat                  = ts_data_selection(tmpdat,'badchans',allbadchans);
% events = tmpevents; save('tmpevents.mat','events'); clear events
% 
% visualizer(tmpdat); % load tmpevents.mat


%% eliminate inconsistent channels
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


%% flip data matrix

% OutFile = strrep(OutFile,'rejection','rejection_flipped');
% fprintf(fid,'Saving file: %s\n',OutFile);

% ...


%% define SO clusters + involved channels
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
showtrl = [183 250 459 653 672];
% 1:length(tc);%[1 11 44 53 64 77 81 93 101]);% showtrl = [37 39 40 44:47 55 65:67 76 81 91:97 99 101]; % []
plot_flag = 0;
normalize = 0;
winsize   = 2 ;                   % size of window to display centered at ref time
pad       = round(winsize*Fs/2);  % number of indices to display around tk
if ~isempty(showtrl)
  meanflipwave = zeros(2*pad+1,length(showtrl));
  sumflipwave  = zeros(1,length(showtrl)); 
  AIC          = zeros(length(showtrl),2);  
end

tic
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
  tmpcat = cat(1,xx(xs>=0,:),-xx(xs<0,:));
  tmpavg = mean(tmpcat,1);
  if plot_flag    
    figure(1); subplot(3,2,1),cla
    plot(tt,xx,'-','LineWidth',2); hold on; axis tight; set(gca,'ylim',ylim); set(gca,'xlim',xlim);
    xlabel('time (s)'); for k=1:length(ts),vline(ts(k),'y'); end  %   plot(ts,xs,'.','markersize',12,'Color','k','linewidth',2); 
    plot(tt,xx(ref,:),'k.-','LineWidth',4);
    vline(tk-cpad,'Color','k','linewidth',3); vline(tk+cpad,'Color','k','linewidth',3);
    vline(tk,'Color','b','linewidth',3); hline(0,'k');
    plot(ts(ref),xs(ref),'o','markersize',16,'linewidth',3,'Color','k');
    title(sprintf('cluster %g (%gs), ref=[%s]',c,tk,[allreflabels{:}])); 
    % flipped
    subplot(3,2,2);plot(tt,xx(xs>0,:),'b-',tt,-xx(xs<0,:),'r-'); title('flipped (red)'); hold on
    plot(tt,tmp,'k','LineWidth',4); set(gca,'ylim',ylim); set(gca,'xlim',xlim);
    plot(ts(xs<0),-xs(xs<0),'.',ts(xs>=0),xs(xs>=0),'.','markersize',16,'Color','k','linewidth',2); 
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
    options = statset('MaxIter',500);
    try
      f1 = gmdistribution.fit(ts',1,'Options',options); AIC(showtrl==c,1) = f1.AIC;
      f2 = gmdistribution.fit(ts',2,'Options',options); AIC(showtrl==c,2) = f2.AIC;
    end
    clear f1 f2    
  end  
end
toc
% save Workspace_s2_plot_clusters_AggHistogram_ClusterWindow0.6sec_filt0.1-4Hz_toi0-8855.14_grad1grad2.mat

%% use cluster to establish how to flip

% create state matrix S = [Sik] where
%   Sik = 1  if chan i peak in cluster k is > 0
%   Sik = -1 if chan i peak in cluster k is < 0
%   Sik = 0  if chan i is not involved in cluster k

% Nij is the number of clusters in which i & j are co-involved

ncluster  = length(tc);
S         = zeros(nchan,ncluster);
rij       = zeros(nchan,nchan);
Nij       = zeros(nchan,nchan);
T         = procdat.epochs.time;
sens      = procdat.sensor_info;
Fs        = procdat.sfreq;
normalize = 0;
winsize   = 2 ;                                                                 % size of window to display centered at ref time
pad       = round(winsize*Fs/2);                                                % number of indices to display around tk
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
  Nij(chans,chans)  = Nij(chans,chans) + 1;
  v         = zeros(1,length(xs));
  v(xs>=0)  =  1;
  v(xs< 0)  = -1;
  rij(chans,chans) = rij(chans,chans) + (v'*v);  
  
end

rij = rij ./ Nij;
% Ri  = 

% IF: ( 0 <= Ri <= 0.75 )
% ... [reject channel i]
% IF: ( 0.75 <= Ri <= 1 )
% ... [calc eik for each cluster k & reject if eik < 0.75]
%     [eik = f(pij) where pij = <.> over all trials]
% [reject any channels for which no trials remain]

% recalc pij after rejection
% select pij for arbitrary channel & use to flip the entire data matrix

% flip data matrix   - use combo of cluster cross-sections & flip matrix

% update event codes - pos/neg now corresponds to UP/DOWN

% save PEAKS         - these peaks will be the focus of subsequent analyses

% visualize flipped data matrix w/ SO clusters for selected intervals
% visualizer(procdat); % load interval-constrained peak clusters w/ new codes

%% analyze clusters

% find channel w detection nearest to the reference time for each cluster
%   => call that the given cluster's reference channel

% process clusters with thus defined reference channels:
% R(delay,distance), single-trial PLV (gamma), coherence, cross-correlogram

%% plots





















