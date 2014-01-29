addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/functions
addpath /space/mdkm1/2/kmdev/projects/jsherfey/sleep/scripts
%% 28-Jun-2010, JSS
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
tic
clear all; loadgrad = 1; loadeeg = 0; writelog = 0;
outpath   = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s1'; cd(outpath);
eventfile = 's1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat';
script_20100628_s1_MappingPolarity2CorticalState_setup
  % main variables:  data (grad), eeg, events, peaks
  % other variables: badlabels, fid, fiffiles, matfiles, subject
toc
% Data selection
testchans = {data.sensor_info.label}; %{'MEG 0143','MEG 0142','MEG 0213','MEG 0212'};
PhiNodes  = linspace(-pi/2,pi/2,40);%-pi/2:pi/48:pi/2;
tlim      = [1700 2700];%[500 2700];%[2000 2500];
sofreq    = [.1 4];
gmfreq    = [20 100];
GammaRes  = []; % gamma will be smoothed over a window = length(epoch)/GammRes
GammaEnv  = 0;
winsize   = 2.5;                                    % window size in seconds
% % note: window should be large enough to include peak and surrounding troughs
% select paired peaks
peaks     = select_peakpairs(peaks,1);
% select subset of peaks & data (test channels & time limits)
script_20100628_s1_MappingPolarity2CorticalState_data_selection
  % T & X         => all times & select chans
  % data & peaks  => select times & select chans
  % others: skipchans (no peaks) & pkt (time vector corresponding to peak indices)
toc
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

% for each channel
for ch = 1:nchan
  label    = testchans{ch};
  for type = 1:length(peaktypes)
    pktype = peaktypes{type};
    pks    = peaks(ch).(pktype);
    % bandpass filter for slow oscillation
    x = ts_freq_filt(X(ch,:)',Fs,sofreq,[0 0],'bandpass')';
    x = x(nearest(T,t(1)):nearest(T,t(end)));
    % bandpass filter for gamma activity    
    g = ts_freq_filt(X(ch,:)',Fs,gmfreq,[0 0],'bandpass')';
    g = ts_freq_filt(g',Fs,60,3,'notch')';
    g = g(nearest(T,t(1)):nearest(T,t(end)));

    % if negpeak => flip x
    if strcmp(pktype,'negpeak'), flip = -1; else flip = 1; end
    x = flip*x;

    % shift pks since length(t) < length(T) & pks comes from T
    pks      = pks - nearest(T,t(1)) + 1;
    % eliminate pks with insufficient data padding
    pks      = pks((pks-pad>=1)&(pks+pad<=ntime));
    % preallocation
    npk      = length(pks);
    xnodes   = zeros(nNode,npk);
    ynodes   = zeros(nNode,npk);
    gnodes   = zeros(nNode,npk);
    skip     = zeros(1,npk);
    GammaMinusYlag = zeros(1,npk);
    % loop over detection times {tk}
    for k = 1:npk
      tki = pks(k);
      tk  = t(tki);
      ii  = tki + (-pad:pad);
      xx  = x(ii);
      gg  = g(ii);
      if GammaEnv,           gg = abs(hilbert(gg)); end
      if ~isempty(GammaRes), gg = smooth(gg,ceil(length(ii)/GammaRes),'lowess'); end
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

      xL = flip*xL;
      yL = flip*yL;

      if debugflag
        % plot to verify monotonic phase; at least verify the crossings surround tk
        figure; 
        subplot(8,1,1),plot(tt,xx,'b'); axis tight; title('so'); hline(L,'k'); vline(tk,'k');
        subplot(8,1,2),plot(tt,xL,'b',tt(ix),xL(ix),'b*'); axis tight; title('xL: so shifted');hline(0,'k');vline(tk,'k');
        subplot(8,1,3),plot(tt,yL,'g',tt(ix),yL(ix),'g*'); axis tight; title('yL: so shifted scaled');hline(0,'k');vline(tk,'k');
        subplot(8,1,4),plot(tt,phi,'k',tt(ix),phi(ix),'k*'); axis tight; title('so phase');hline(0,'k');vline(tk,'k');
        subplot(8,1,5),plot(tt,gg,'r',tt(ix),gg(ix),'r*'); axis tight; title('smooth gamma filt');hline(0,'k');vline(tk,'k');
        subplot(8,1,6),plot(PhiNodes,xL(ix),'b*-'); axis tight; title('xL on phi');hline(0,'k');vline(median(PhiNodes),'k');
        subplot(8,1,7),plot(PhiNodes,yL(ix),'g*-'); axis tight; title('yL on phi');hline(0,'k');vline(median(PhiNodes),'k');
        subplot(8,1,8),plot(PhiNodes,gg(ix),'r*-'); axis tight; title('gamma on phi');hline(0,'k');vline(median(PhiNodes),'k');
        figure('Name','v = predefined phi values');
        [tmp,lags] = xcorr(yL(ix),gg(ix),'coeff'); dphi=PhiNodes(2)-PhiNodes(1);
        subplot(3,1,1),plot(lags*dphi,tmp); axis tight; title('xcorr(yL(v),gamma(v))');
        [tmp,lags] = xcorr(phi(ix),yL(ix),'coeff');
        subplot(3,1,2),plot(lags*dphi,tmp); axis tight; title('xcorr(phi(v),yL(v))');
        tmp1=lags(max(tmp)==tmp); vline(tmp1*dphi,'k');
        [tmp,lags] = xcorr(phi(ix),gg(ix),'coeff');
        subplot(3,1,3),plot(lags*dphi,tmp); axis tight; title('xcorr(phi(v),gamma(v))');
        tmp2=lags(max(tmp)==tmp); vline(tmp2*dphi,'k');
        xlabel(sprintf('gamma\\_lag - yL\\_lag = %g radians\n',(tmp2-tmp1)*dphi));
        pause; close
      end
      % get x-values at PhiNodes
      xnodes(:,k) = xL(ix);
      ynodes(:,k) = yL(ix);
      gnodes(:,k) = gg(ix);

      % cross-correlations
      dphi       = PhiNodes(2)-PhiNodes(1);
      [tmp,lags] = xcorr(phi(ix),yL(ix),'coeff'); tmp1=min(lags(max(tmp)==tmp));
      [tmp,lags] = xcorr(phi(ix),gg(ix),'coeff'); tmp2=min(lags(max(tmp)==tmp));
      GammaMinusYlag(k) = (tmp2-tmp1)*dphi;
  %     [tmp,lags] = xcorr(yL(ix),gg(ix),'coeff');
      clear tmp tmp1 tmp2
    end
    peaks(ch).([pktype '_skip'])           = skip;
    peaks(ch).([pktype '_GammaMinusYlag']) = GammaMinusYlag;
    peaks(ch).PhiNodes             = PhiNodes;
    peaks(ch).([pktype '_xnodes']) = xnodes;
    peaks(ch).([pktype '_ynodes']) = ynodes;
    peaks(ch).([pktype '_gnodes']) = gnodes;
    % calculate tentative average (no outlier rejection)
    peaks(ch).([pktype '_xnode_avg']) = sum(xnodes(:,~skip),2)/sum(~skip);
    peaks(ch).([pktype '_ynode_avg']) = sum(ynodes(:,~skip),2)/sum(~skip);
    peaks(ch).([pktype '_gnode_avg']) = sum(gnodes(:,~skip),2)/sum(~skip);
    clear xnodes ynodes gnodes xx gg tt xL yL phi x g pks ix GammaMinusYlag skip
  end % end loop over peak types
  toc
end
nkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
% set outliers to zero
% get nNode x N's from # of non-zero elements along 2nd dimension
% average over remaining peaks

posgamma_integral = arrayfun(@(x)(sum(x.pospeak_gnode_avg-min([x.pospeak_gnode_avg;x.pospeak_gnode_avg]))),peaks);
posgamma_integral(posgamma_integral==0) = [];
neggamma_integral = arrayfun(@(x)(sum(x.negpeak_gnode_avg-min([x.pospeak_gnode_avg;x.pospeak_gnode_avg]))),peaks);
neggamma_integral(neggamma_integral==0) = [];
%   figure; plot(1:nchan,posgamma_integral,'b.',1:nchan,neggamma_integral,'ro'); title('gamma integrals');
%   legend('pos','neg'); [posgamma_integral' neggamma_integral']
diffgam  = posgamma_integral-neggamma_integral; 
chans    = 1:nchan; 
thresh   = median(abs(diffgam))+2*std(abs(diffgam));
badchans = abs(diffgam)>(thresh);
  figure('Name','Gamma Integral: pos - neg'); subplot(2,1,1),plot(chans,diffgam,'.-',chans(badchans),diffgam(badchans),'ro');
  axis tight; hline(0,'k'); hline(-thresh,'r'); hline(thresh,'r');
  subplot(2,1,2),plot(chans(~badchans),diffgam(~badchans),'.-'); title('without outlier channels');
  axis tight; set(gca,'ylim',[-1 1]*(1.2*thresh)); hline(0,'k'); hline(-thresh,'r'); hline(thresh,'r');
  {peaks(badchans).label}

% nposkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
% nnegkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
% [nposkeep' nnegkeep'] % equivalent
% nkeep = arrayfun(@(x)(sum(~x.pospeak_skip)),peaks);
% figure; subplot(1,2,1),plot(nkeep,'.'); level1 = median(nkeep); hline(level1,'r');
% tmp = .5*std(nkeep); hline(level1+tmp,'g'); hline(level1-tmp,'g');
% subplot(1,2,2),try hist(nkeep,round(length(nkeep)/10)); end;
% HighFreqChans = [1:7 100:110]; LowFreqChans  = 45:55; % Examples

% two subplots: xnode & ynode; each w/ a pospeak/negpeak overlay
% each curve has nNode points
thresh = median(nkeep) - .5*std(nkeep);
nn = 10; startchan = 1; highlight = find(nkeep>thresh);
if 1%0
  nsel     = length(highlight);
  sodiff   = zeros(nchan,2);
  MaxPos   = zeros(nchan,1); % is the max neg or pos
  MaxLeft  = zeros(nchan,1); % is the max left or right
  evcode   = zeros(nchan,1);
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
      y21 = peaks(ch).negpeak_gnode_avg - off; y22 = peaks(ch).pospeak_gnode_avg - off;
      x   = peaks(ch).PhiNodes*(180/pi);
      tmp = y12-abs(y11);
      sodiff(ch,1) = sum(tmp(x<0));
      sodiff(ch,2) = sum(tmp(x>0));
      % get max vals
      ix           = find(abs(tmp)==max(abs(tmp)));
      MaxLeft(ch)  = PhiNodes(ix)<0;
      MaxPos(ch)   = tmp(ix)>0;
      evcode(ch,:) = str2num(sprintf('%g%g',[MaxLeft(ch)+1 MaxPos(ch)+1]));
      if MaxLeft(ch), str='left';       else str='right';       end
      if MaxPos(ch),  str=[str '-pos']; else str=[str '-neg'];  end
      str = [str ': ' num2str(evcode(ch))];
      subplot(4,nn,cnt),plot(x,y12,'b',x,y11,'r'); axis tight; 
      title(sprintf('%s:n=%g/%g',peaks(ch).label,sum(~peaks(ch).pospeak_skip),sum(~peaks(ch).negpeak_skip))); 
      vline(0,'k'); if cnt==1, legend('pos','neg'); end
      subplot(4,nn,nn+cnt),plot(x,y12-abs(y11),'k'); axis tight; hline(0,'k'); vline(0,'k'); title(str); axis off % title('SO: pos-|neg|)');
      subplot(4,nn,2*nn+cnt),plot(x,y22,'b',x,y21,'r'); axis tight; title('mean gamma'); vline(0,'k');
      subplot(4,nn,3*nn+cnt),plot(x,smooth(y22-y21,5,'lowess'),'k'); axis tight; title('gamma: pos-neg'); hline(0,'k'); vline(0,'k');
      xlabel('instantaneous phase of SO (Hz)'); % maybe smooth it
    end
  end
%   figure
%   subplot(3,1,1),plot(sodiff(:,1),'.'); title('sum over pos-|neg| for phi<0');
%   subplot(3,1,2),plot(sodiff(:,2),'.'); title('sum over pos-|neg| for phi>0');
%   subplot(3,1,3),plot(sodiff(:,1)-sodiff(:,2),'.'); title('(sum for phi<0)-(sum for phi>0)');
%   xlabel('channel');
  clear sodiff
end
for k=1:length(evcode),fprintf('%3g [%s]\n',evcode(k),peaks(k).label); end

poscode   = [12 22]; negcode   = [11 21];
leftcode  = [21 22]; rightcode = [11 12];
evcodes   = unique(evcode)';
ev0 = {peaks(evcode==0).label}';  n0 = length(ev0);
ev1 = {peaks(evcode==11).label}'; n1 = length(ev1);
ev2 = {peaks(evcode==21).label}'; n2 = length(ev2);
ev3 = {peaks(evcode==12).label}'; n3 = length(ev3);
ev4 = {peaks(evcode==22).label}'; n4 = length(ev4);
disp([evcodes;n0 n1 n2 n3 n4])
% visualizer(ts_data_selection(data,'chanlabel',ev0));
% visualizer(ts_data_selection(data,'chanlabel',ev1));
% visualizer(ts_data_selection(data,'chanlabel',ev2));
% visualizer(ts_data_selection(data,'chanlabel',ev3));
% visualizer(ts_data_selection(data,'chanlabel',ev4));

figure
highlight = find(nkeep>thresh); subplot(1,2,1); dewar; title('many detections');
highlight = find(nkeep<thresh); subplot(1,2,2); dewar; title('few detections');

% II. Classifying polarities: use SO phase-locked averaging of smoothed, 
% Hilbert-based gamma envelope.

% level-crossing analysis
stepsize  = 2;          % # of points to shift each step
npts      = 5;          % # of points in the window
ngrid     = 100;
ch        = 1;
x         = peaks(ch).PhiNodes*(180/pi);
nx        = length(x);
steps     = 1:stepsize:[nx-npts+1];
nsteps    = length(steps);
allcounts = {};

for type = 1:length(peaktypes)
  pktype  = peaktypes{type};%'pospeak';
  keep    = find(~(peaks(ch).pospeak_skip | peaks(ch).negpeak_skip));
  pks     = peaks(ch).(pktype);
  pks     = pks(keep);
  npks    = length(pks);
  count   = zeros(npks,nsteps);
  for k = 1:npks
    % select gamma=f(phi) around this peak
    G        = peaks(ch).([pktype '_gnodes'])(:,keep(k)); % gamma = f(preset_phi)
    xi       = linspace(PhiNodes(1),PhiNodes(end),ngrid);
    Gi       = interp1(PhiNodes,G,xi,'spline');
    L        = median(Gi)+std(Gi);%.9*max(Gi);
    [ind,t0] = crossing(Gi,xi,L);
    if k == 1, figure; plot(PhiNodes,G,'o',xi,Gi,'-',xi(ind),Gi(ind),'g*'); axis tight; hline(L,'r'); title(pktype); xlabel('phi (rad)'); ylabel('gamma'); end
    for j = 1:nsteps
      this       = find((ind>=steps(j))&(ind<=(steps(j)+npts-1)));
      count(k,j) = length(this);
    end
    clear G xi Gi L ind t0 this
  end
  allcounts{type} = count;
  clear count pks npks
end
toc

tt = PhiNodes(steps+floor(stepsize/2))*(180/pi); % t(steps);
figure; plot(tt,sum(allcounts{1},1),'b*-',tt,sum(allcounts{2},1),'ro-'); axis tight
legend('pos','neg'); title('high-level crossing gamma counts'); vline(0,'k');
set(gca,'xlim',[-90 90]);

clear tmp
fld   = 'gnodes';
ch    = 109;% 110; highlight(56); %LowFreqChans(8);
sel   = ~(peaks(ch).pospeak_skip | peaks(ch).negpeak_skip);
field = ['pospeak_' fld]; tmp(1) = ts_matrix2data(peaks(ch).(field)(:,sel),'time',PhiNodes); tmp(1).epochs.event_code=1;
field = ['negpeak_' fld]; tmp(2) = ts_matrix2data(peaks(ch).(field)(:,sel),'time',PhiNodes); tmp(2).epochs.event_code=2;
tmp = ts_combine_data(tmp);
% visualizer(tmp)
% take difference
tmp.epochs(1).data = tmp.epochs(1).data - tmp.epochs(2).data;
tmp.epochs         = tmp.epochs(1);
% visualizer(tmp);
% concatenate epochs to make continuous
% c=1;  tmpdat = (tmp.epochs(c).data(:)); %tmpdat(1:length(PhiNodes):length(tmpdat)) = max(tmpdat); 
%       tmp.epochs(c).data = tmpdat'; tmp.epochs(c).time = (0:length(tmpdat)-1)/40;%(0:length(tmpdat)-1)*(PhiNodes(2)-PhiNodes(1))*(180/pi);
% c=2;  tmpdat = (tmp.epochs(c).data(:)); %tmpdat(1:length(PhiNodes):length(tmpdat)) = min(tmpdat);
%       tmp.epochs(c).data = tmpdat'; tmp.epochs(c).time = (0:length(tmpdat)-1)/40;%(0:length(tmpdat)-1)*(PhiNodes(2)-PhiNodes(1))*(180/pi);
% visualizer(tmp);

% epoch diff
tmpdiff = peaks(ch).(['pospeak_' fld])(:,sel) - abs(peaks(ch).(['negpeak_' fld])(:,sel));
tmpdiff = ts_matrix2data(tmpdiff','time',PhiNodes,'continuous',1);
% visualizer(tmpdiff)

% avg
fld   = 'gnode_avg';
field = ['pospeak_' fld]; tmp(1) = ts_matrix2data(peaks(ch).(field),'time',PhiNodes,'continuous',1); tmp(1).epochs.event_code=1;
field = ['negpeak_' fld]; tmp(2) = ts_matrix2data(peaks(ch).(field),'time',PhiNodes,'continuous',1); tmp(2).epochs.event_code=2;
tmp = ts_combine_data(tmp);
% visualizer(tmp)

tmp1 = arrayfun(@(x)(mean(x.pospeak_GammaMinusYlag(~x.pospeak_skip))),peaks);
tmp2 = arrayfun(@(x)(mean(x.negpeak_GammaMinusYlag(~x.negpeak_skip))),peaks);
ii   = 1:length(tmp1);
figure('Name','Lag: gamma minus shifted half-wave');
subplot(3,1,1),plot(ii,tmp1,'bo'); axis tight; hline(0,'k'); title('pos')
subplot(3,1,2),plot(ii,tmp2,'r*'); axis tight; hline(0,'k'); title('neg')
subplot(3,1,3),plot(ii,tmp1-tmp2,'k.'); axis tight; hline(0,'k'); title('pos-neg'); xlabel('channel')

% plot average for one channel
      cnt = 1; nn = 1; figure
      off = 0;%min([peaks(ch).pospeak_xnode_avg;peaks(ch).pospeak_xnode_avg]);
      y11 = peaks(ch).negpeak_xnode_avg - off; y12 = peaks(ch).pospeak_xnode_avg - off;
      off = 0;%min([peaks(ch).pospeak_gnode_avg;peaks(ch).pospeak_gnode_avg]);
      y21 = peaks(ch).negpeak_gnode_avg - off; y22 = peaks(ch).pospeak_gnode_avg - off;
      x  = peaks(ch).PhiNodes*(180/pi);
      subplot(4,nn,cnt),plot(x,y12,'b',x,y11,'r'); axis tight; title(sprintf('%s:n=%g/%g',peaks(ch).label,sum(~peaks(ch).pospeak_skip),sum(~peaks(ch).negpeak_skip))); legend('pos','neg'); vline(0,'k');
      subplot(4,nn,nn+cnt),plot(x,y12-abs(y11),'k'); axis tight; title('SO: pos-|neg|)'); hline(0,'k'); vline(0,'k');
      subplot(4,nn,2*nn+cnt),plot(x,y22,'b',x,y21,'r'); axis tight; title('mean gamma'); legend('pos','neg'); vline(0,'k');
      subplot(4,nn,3*nn+cnt),plot(x,smooth(y22-y21,5,'lowess'),'k'); axis tight; title('gamma: pos-neg'); hline(0,'k'); vline(0,'k');
      xlabel('instantaneous phase of SO (Hz)'); % maybe smooth it
      fld  = 'gnodes';
      sel  = ~(peaks(ch).pospeak_skip | peaks(ch).negpeak_skip);
      Zavg = y12-abs(y11);
      Zdat = peaks(ch).(['pospeak_' fld])(:,sel) - abs(peaks(ch).(['negpeak_' fld])(:,sel));;
      ZL   = Zdat(x<0,:);
      ZR   = Zdat(x>0,:);
      ZLs  = sum(ZL,1);
      ZRs  = sum(ZR,1);
      n = 100; [ZLs(n),ZRs(n),ZLs(n)-ZRs(n)]
      
      Z = ts_matrix2data([ZLs;ZRs;ZLs-ZRs;ZL;ZR],'time',1:length(ZRs),'continuous',1);
%       visualizer(Z)

%%
% %% Test #2 => Epoch method #3 + Classification method #2.
% 
% % Load data
% clear all; loadgrad = 1; loadeeg = 0; writelog = 0;
% outpath   = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s1'; cd(outpath);
% eventfile = 's1_SO_init_peaks_filt0.01-4Hz_toi0-6350.11_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat';
% script_20100628_s1_MappingPolarity2CorticalState_setup
% 
% % I. Create epochs: use window with fixed size centered on paired detection times.
% peaks     = select_peakpairs(peaks,1);
% testchans = {'MEG 0143','MEG 0142','MEG 0213','MEG 0212'};
% tlim      = [500 2700];
% script_20100628_s1_MappingPolarity2CorticalState_data_selection
% 
% t         = data.epochs.time;
% ntime     = length(t);
% npeak     = cellfun(@length,{peaks.pospeak});
% nchan     = data.num_sensors;
% Fs        = data.sfreq;
% peaktypes = {'pospeak','negpeak'};
% 
% 
% 
% %% Test #3 => Epoch method #3 + Classification method #3.
% 
% 
