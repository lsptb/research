function [clusters,err] = SO_cluster_corrcoef(data,detections,varargin);
% Purpose:
% ...
% Inputs:
% - gridscale, center-to-center distance [m] between ECoG macroelectrodes.
% 
% Modified by JSS on 10-Aug-2010
% - changed R-calc so that if a cluster contains both pos and neg peaks, R
% will be calculated using all involved channels and only one "condition"
% will be returned.
% Brain Volume must be given in cm^3

params = mmil_args2parms( varargin, ...
                         { 'cluster_RefChanType','EarliestDetection',{'EarliestDetection','HistPeak','Closest2HistPeak','MaxR'},...
                           'layout',[],[],...
                           'cluster_ClusterWindow',[],[],...
                           'corrcoef_MinChansPerCluster',20,[],...
                           'preproc',1,{0,1},...
                           'BrainVolume',1450,[],...
                           'iEEG_flag',0,{0,1},...
                           'badlabels',[],[],...
                           'SulciFactor',3,[],...
                           'gridscale',.01,[],...
                         }, ...
                         false );

t         = detections(1).tstart:1/detections(1).sfreq:detections(1).tstop;          % detections time vector
T         = data.epochs.time;                                                        % data time vector
Toffset   = nearest(t,T(1));
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
tc          = t(clusterIX);  

nchan   = data.num_sensors;
pos     = zeros(nchan,3); % (x,y,z) for each channel
if isempty(params.layout)
  trans = data.coor_trans.device2head;
  for k = 1:nchan
    loc         = data.sensor_info(k).loc;
    if any(strmatch('grad',data.sensor_info(k).typestring)) || ...
       any(strmatch('mag' ,data.sensor_info(k).typestring))
      loc       = trans*loc;
    end  
    pos(k,1:3)  = loc(1:3,4);
  end
  method = 'stereographic';%params.method; % gnomic, stereographic, ortographic, inverse, polar
  prj    = elproj(pos, method); % * [0 1; -1 0];
            % ELPROJ makes a azimuthal projection of a 3D electrode cloud
            %  on a plane tangent to the sphere fitted through the electrodes
            %  the projection is along the z-axis
  X = prj(:,1);   % x-coordinates
  Y = prj(:,2);   % y-coordinates  
elseif ischar(params.layout) && exist(params.layout,'file')
  [chNum,X,Y,Width,Height,Lbl,Rem] = textread(params.layout,'%f %f %f %f %f %q %q');
  for i=1:length(Lbl)
    if ~isempty(Rem{i})
      % this ensures that channel names with a space in them are also supported (i.e. Neuromag)
      Lbl{i} = [Lbl{i} ' ' Rem{i}];
    end
  end
  [sel,jnk] = match_str(Lbl,{data.sensor_info.label});
  X = X(sel);
  Y = Y(sel);
elseif isstruct(params.layout)
  X = lay.pos(:,1);
  Y = lay.pos(:,2);
else
  error('Layout not found.');
end

% cfg.layout =[];
% cfg.layout = params.layout;
% cfg.layout = prepare_layout(cfg); % create 2d layout from 3d positions
% distance2D = dist(cfg.layout.pos(:,1:2)');
distance2D = dist([X Y]');
              % distance(i,j) = sqrt( (x(i)-x(j)).^2 + (y(i) - y(j)).^2 );
% Angular distance
zshift     = .08; % shifts the Z center to almost the middle of the helmet
theta3D    = ts_BetweenSensorAngles(data.sensor_info,zshift);
              % theta(i,j) = real(acos(dot(Pi,Pj)/(norm(Pi)*norm(Pj)))*180/pi)
              %              where Pk = (x,y,z)-coord of k-th sensor
theta3D(isnan(theta3D)) = 0;
  % [theta] = nchan x nchan (same for distance2D)

% remove rejected chans from distance matrices
[sel1,sel2] = match_str({detections.label},{data.sensor_info.label});
distance2D  = distance2D(sel1,sel1);
theta3D     = theta3D(sel1,sel1);

% use Euclidean distance for iEEG and angular distance for M/EEG
if params.iEEG_flag 
  distance  = distance2D;
  % rescale [0,1] - arbitrary units
  distance  = distance / max(distance(:));
else
  distance  = theta3D;
end

% make sure distance matrices have the right number of channels
if size(distance,1) ~= nchan, error('size of distance matrix is incorrect.'); end

% Update info on clusters post-rejection
% create array of cluster numbers
cnum        = {detections.cluster_number};
cnum        = unique([cnum{:}]);

[a,cnum,b]  = intersect(cnum,clusterID);
clear a b

ncluster      = length(cnum);         % number of clusters with non-rejected channels
tc            = tc(cnum);             % update cluster reference times
InvolvedChans = InvolvedChans(cnum);  % update list of involved channels
cpad          = params.cluster_ClusterWindow / 2;
if params.preproc
  data = ts_preproc(data,'bpfilter','yes','bpfreq',[.1 4],'bandpass_detrend_flag',0,'blc','yes','blcwindow',[-inf inf]);
end
% split clusters into UP & DOWN (don't know which is UP, but one is)
clusters = [];
clusters(1).peaktype = 'pospeak';
clusters(2).peaktype = 'negpeak';
clusters(1).code = 1; % this will be pospeak
clusters(2).code = 2; % this will be negpeak
[clusters(1:2).cluster_number] = deal([]);

% get indices to bad channels (into detections and data)
if ~isempty(params.badlabels)
  [TF,badchans] = ismember(params.badlabels,{data.sensor_info.label});
  badchans      = badchans(badchans~=0);
else
  badchans      = [];
end
if ~isequal({data.sensor_info(badchans).label},{detections(badchans).label})
  error('Inconsistent channel indices.');
end

% Loop over clusters
warning('off','MATLAB:divideByZero');
err   = zeros(1,ncluster);
for k = 1:ncluster
  if Ninvolved(k) < params.corrcoef_MinChansPerCluster
    err(k) = 1;
    continue;
  end
  chans = InvolvedChans{k};
  chans = setdiff(chans,badchans);
  % ... this may be the place to remove bad channels...
  if max(cat(1,InvolvedChans{:})) > length(detections)
    err(k) = 1;
    error('Channel indices are inconsistent.');
%     [chans jnk] = match_str({detections.label},{orig_detections(chans).label});
  end
  if isempty(chans), err(k) = 1; continue; end
  tk    = tc(k);
  % find indices for each involved channel
  ix1 = []; ix2 = [];
  if ~isempty(detections(chans(1)).pospeak)
    s1 = arrayfun(@(x)x.pospeak(nearest(t(x.pospeak),tk)),detections(chans)); ix1 = abs(t(s1)-tk)<=cpad;
  end
  if ~isempty(detections(chans(1)).negpeak)
    s2 = arrayfun(@(x)x.negpeak(nearest(t(x.negpeak),tk)),detections(chans)); ix2 = abs(t(s2)-tk)<=cpad;
  end
  % should be either pospeak or negpeak for all chans
  CHANS   = chans;
  if      (~isempty(ix1) && all(ix1)) && (isempty(ix2) || ~all(ix2))
    codes = 1;
    IX    = s1(ix1);
  elseif  (~isempty(ix2) && all(ix2)) && (isempty(ix1) || ~all(ix1))
    codes = 2;
    IX   = s2(ix2);
  elseif (~isempty(ix1) && ~isempty(ix2)) && all((ix1+ix2)>=1)
    % split into pos and neg
    codes = [1 2];
    IX{1} = s1(ix1);
    IX{2} = s2(ix2);
    clear CHANS
    CHANS{1} = chans(ix1);
    CHANS{2} = chans(ix2);
  else
    err(k) = 1;
    continue;
  end
  for c  = 1:length(codes)
    code = codes(c);
    if iscell(IX),    ix    = IX{c};    else ix    = IX;    end
    if iscell(CHANS), chans = CHANS{c}; else chans = CHANS; end
    if length(chans) < 2, continue; end
    AbsIndex = ix;
    SelIndex = ix - Toffset + 1; 
%     SelIndex = cellfun(@(x)nearest(T,t(x)),num2cell(ix));                         % convert to data time indices
    if strcmp(params.cluster_RefChanType,'EarliestDetection')
      relref    = find(min(AbsIndex)==AbsIndex);                                    % peak index to 1st detection in the cluster
    elseif strcmp(params.cluster_RefChanType,'Closest2HistPeak')
      relref    = nearest(t(AbsIndex),tk);                                          % peak index for each chan nearest to tk in the data time vector
    elseif strcmp(params.cluster_RefChanType,'HistPeak')
      error('This option does not work since a distance cannot be defined.');
    elseif strcmp(params.cluster_RefChanType,'MaxR')
      % Calc & store (R,p) using all chans as ref then define ref as the
      % chan within 50ms of earliest det which produces the max R value
      % Loop over involved chans
      ts      = cellfun(@(x)T(x),num2cell(SelIndex));                                        % peak time for each chan from the data time vector
      try
        xs    = cellfun(@(x,y)data.epochs.data(x,y),num2cell(chans),num2cell(SelIndex)); % NOTE: chan k has peak at (ts(k),xs(k))
      catch
        xs    = cellfun(@(x,y)data.epochs.data(x,y),num2cell(chans),num2cell(SelIndex)'); % NOTE: chan k has peak at (ts(k),xs(k))
      end
      Nk      = length(chans);
      Rvals   = zeros(1,Nk);
      pvals   = zeros(1,Nk);
      for ch  = 1:Nk
        this  = chans(ch);
        D     = distance(this,chans);   % angular distance from reference
        tau   = abs(ts - ts(ch));            % delay wrt reference
        [tmpR,tmpP] = corrcoef([D' tau']);  % correlation coefficient b/w distance & delay wrt ref
        if length(tmpR)>1 && ~isnan(tmpR(1,2))
          Rvals(ch) = tmpR(1,2);
          pvals(ch) = tmpP(1,2);
        end
      end
      % chans within 50ms of the earliest detection
      potref = find((ts-min(ts)) < .05);
      % the chan within 50ms with the largest |R|
      relref = potref(abs(Rvals(potref))==max(abs(Rvals(potref))));
    end
    if length(relref) > 1, relref = relref(1); end
    ref       = chans(relref);
    reflabel  = detections(ref).label;
    evnt      = find([clusters.code]==code);                                          % UP or DOWN
    cnt       = length(clusters(evnt).cluster_number) + 1;                            % # of this type
    if cnt == 1
      clusters(evnt).sensor_info  = data.sensor_info;
      clusters(evnt).matfiles     = data.epochs.matfiles;
      clusters(evnt).duration     = data.epochs.duration;
      clusters(evnt).abs_toilim   = [detections(1).tstart detections(1).tstop];
      clusters(evnt).sel_toilim   = [data.epochs.time(1) data.epochs.time(end)];
      clusters(evnt).sfreq        = data.sfreq;
      if params.iEEG_flag
        clusters(evnt).deg2meters = params.gridscale / ((max(X(:))-min(X(:)))/sqrt(nchan));
        % Assumption: perfect nxn grid
      else
        clusters(evnt).deg2meters = (2*pi*((3*params.BrainVolume/(4*pi))^(1/3))/360)*(1/100)*params.SulciFactor;
        % X [deg]     = deg2meters*X      [m]
        % X [deg/sec] = deg2meters*X      [m/s] or [mm/ms]
        % X [deg/ms]  = deg2meters*X*1000 [m/s] or [mm/ms]
      end
      if isfield(data.epochs,'IntervalEndPoints'), clusters(evnt).IntervalEndPoints = data.epochs.IntervalEndPoints; end
    end
    ts  = cellfun(@(x)T(x),num2cell(SelIndex));                                        % peak time for each chan from the data time vector
    try
      xs = cellfun(@(x,y)data.epochs.data(x,y),num2cell(chans),num2cell(SelIndex));
    catch
      xs = cellfun(@(x,y)data.epochs.data(x,y),num2cell(chans),num2cell(SelIndex)'); % NOTE: chan k has peak at (ts(k),xs(k))
    end
    D   = distance(ref,chans);  % angular distance from reference
    tau = abs(ts - ts(relref));      % delay wrt reference
    if params.iEEG_flag
      x = D;
    else      
      r = (3*params.BrainVolume/(4*pi))^(1/3);
      x = D*(2*pi*r/360)*(1/100)*params.SulciFactor;
    end
    if tau(x~=0) == 0
      s = inf;
    else
      s = mean(x(tau~=0)./tau(tau~=0));
    end
    [tmp,P] = corrcoef([D' tau']);  % correlation coefficient b/w distance & delay wrt ref
    reftime = ts(relref);
    clusters(evnt).cluster_number(cnt)        = cnum(k);
    clusters(evnt).epochs(cnt).HistTime       = tk;
    clusters(evnt).epochs(cnt).RefChan        = reflabel;
    clusters(evnt).epochs(cnt).RefTime        = reftime;
    clusters(evnt).epochs(cnt).InvolvedChans  = {detections(chans).label};
    clusters(evnt).epochs(cnt).DetectionTimes = ts;
    clusters(evnt).epochs(cnt).DetectionAmps  = xs;
    clusters(evnt).epochs(cnt).Distance        = D;
    clusters(evnt).epochs(cnt).Delays         = tau;
    clusters(evnt).epochs(cnt).Speed          = s;
    if length(tmp)>1 && ~isnan(tmp(1,2))
      clusters(evnt).epochs(cnt).CorrCoef       = tmp(1,2);
    end
    clusters(evnt).epochs(cnt).AbsIndex       = AbsIndex;
    clusters(evnt).epochs(cnt).SelIndex       = SelIndex;
    if strcmp(params.cluster_RefChanType,'EarliestDetection')
      % run calculation once since all delays will be positive
      if length(tmp)>1 && ~isnan(tmp(1,2))
        clusters(evnt).epochs(cnt).R              = tmp(1,2);
        clusters(evnt).epochs(cnt).p              = P(1,2);
        clusters(evnt).epochs(cnt).N              = length(chans);
      end
    elseif strcmp(params.cluster_RefChanType,'Closest2HistPeak')
      % run calculation twice (once for posDelay & negDelay)
      [tmpR,tmpP] = corrcoef([D(tau>=0)' tau(tau>=0)']);
      if length(tmpR)>1 && ~isnan(tmpR(1,2))
        clusters(evnt).epochs(cnt).R_posDelay     = tmpR(1,2);
        clusters(evnt).epochs(cnt).p_posDelay     = tmpP(1,2);
        clusters(evnt).epochs(cnt).N_posDelay     = length(find(tau>=0));
      end
      [tmpR,tmpP] = corrcoef([D(tau<0)' tau(tau<0)']);
      if length(tmpR)>1 && ~isnan(tmpR(1,2))
        clusters(evnt).epochs(cnt).R_negDelay     = tmpR(1,2);
        clusters(evnt).epochs(cnt).p_negDelay     = tmpP(1,2);
        clusters(evnt).epochs(cnt).N_negDelay     = length(find(tau<0));
      end
    elseif strcmp(params.cluster_RefChanType,'MaxR')
      if length(tmp)>1 && ~isnan(tmp(1,2)) && ~isnan(P(1,2))
        if ~(tmp(1,2)==Rvals(relref) && P(1,2)==pvals(relref))
          error('problem with MaxR: vals from loop do not match ref val being stored');
        end
        clusters(evnt).epochs(cnt).RefTimeFromFirst = reftime - min(ts);
        clusters(evnt).epochs(cnt).R                = tmp(1,2);
        clusters(evnt).epochs(cnt).p                = P(1,2);
        clusters(evnt).epochs(cnt).N                = length(chans);
        clusters(evnt).epochs(cnt).R_allref         = Rvals;
        clusters(evnt).epochs(cnt).p_allref         = pvals;
        clear Rvals pvals
      end
    end
    % calc correlation for shuffled locations
    [tmpR,tmpP] = corrcoef([D(randperm(length(D)))' tau']);
    if length(tmpR)>1 && ~isnan(tmpR(1,2))
      clusters(evnt).epochs(cnt).R_RandomLoc    = tmpR(1,2);
      clusters(evnt).epochs(cnt).p_RandomLoc    = tmpP(1,2);
      clusters(evnt).epochs(cnt).N_RandomLoc    = length(find(tau>=0));
    end
  end
  clear IX
end
warning('on','MATLAB:divideByZero');
