  Fs      = flipdat.sfreq;
  bpfreq  = [8 12];
  foi     = [8:12];
  NCO     = 5; width=NCO; %linspace(3,8,length(foi)); width=NCO; %[]; Nco     = 7; % width in freqanalysis
  sf      = foi / width;
%   sf      = [2 2 3 5 5 5 5 5 5]; width = foi ./ sf;
  tpad    = .7;
  tcut    = .5;
  st      = width ./ foi;
  
  t       = flipdat.epochs.time;
  pad     = round(tpad*flipdat.sfreq); pad = -pad:pad;
  nfft    = 2^(nextpow2(length(pad)));
  pad     = -nfft/2:(nfft/2-1);
  sel     = round(Fs*(tpad-tcut));     sel = sel:length(pad)-sel;
  nfreq   = length(foi);
  ntime   = length(sel);
  tic
    warning off
    % preallocation    
    for c = 1:length(clusters)
      for k = 1:length(clusters(c).epochs)
        clusters(c).epochs(k).TFR = complex(zeros(length(clusters(c).epochs(k).InvolvedChans),ntime,nfreq));
      end
    end
    % calculate wavelet spectra for each cluster
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
        xx.sfreq = round(xx.sfreq);
        tfr = ts_freqanalysis_fieldtrip(xx,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0,'verbose',0);
        clusters(c).epochs(k).TFR = tfr.timefreq.cmplx(:,sel,:);        
      end
%       clusters(c).parms.RiThreshold = TF_RiThreshold;
      clusters(c).parms.BandpassFc  = bpfreq;
      clusters(c).parms.MorletFreqs = foi;
      clusters(c).parms.MorletSF    = sf;
      clusters(c).parms.MorletST    = st;
      clusters(c).parms.MorletWidth = width;
      clusters(c).parms.DataPad     = tpad;
      clusters(c).parms.DataSel     = tcut;
    end
    warning on
toc
tic  
  tfr                 = tfr(1);
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
      if isempty(cix)
        tfr(c).sensor_info(i).badchan = 1;
        continue;
      end
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
TFR = ts_data_selection(TFR); % remove badchans


%% adjust for averaging on channel peaks
tic
% preallocation    
for c = 1:length(clusters)
  for k = 1:length(clusters(c).epochs)
    clusters(c).epochs(k).TFR2 = complex(zeros(length(clusters(c).epochs(k).InvolvedChans),ntime,nfreq));
  end
end
toc
Fs    = flipdat.sfreq;
for c = 1:length(clusters)
  for k = 1:length(clusters(c).epochs)
        [sel1,sel2] = match_str(clusters(c).epochs(k).InvolvedChans,{flipdat.sensor_info.label});
        for ch = 1:length(sel1)
          chan = sel1(ch);
          % find offset between HistTime and ch1 (Ref) for this trial
          HistTime    = clusters(c).epochs(k).HistTime;
          RefTime     = clusters(c).epochs(k).DetectionTimes(chan);
          OffsetTime  = RefTime - HistTime;
          Offset      = round(OffsetTime*Fs);
          MidPoint    = floor(ntime/2+1);
          RefPoint    = MidPoint + Offset;                
          % shift res to account for chan peak offset relative to HistTime
          if Offset < 0, Offset = abs(Offset); end
          wi                                    = MidPoint - Offset;
          clusters(c).epochs(k).TFR2(chan,MidPoint-wi+1:MidPoint,:)    = clusters(c).epochs(k).TFR(chan,RefPoint-wi+1:RefPoint,:);
          clusters(c).epochs(k).TFR2(chan,MidPoint+1:MidPoint+wi-2,:)  = clusters(c).epochs(k).TFR(chan,RefPoint+1:RefPoint+wi-2,:);
        end
  end
end
% re-average
toc
  tfr2                 = tfr(1);
  tfr2.timefreq        = tfr2.timefreq(1);
  nchan                = length(peaks);
  tfr2.sensor_info     = flipdat.sensor_info;
  tfr2.num_sensors     = flipdat.num_sensors;
  try tfr2.timefreq    = rmfield(tfr2.timefreq,{'power','cmplx'}); end
  tfr2.timefreq.time   = tfr2.timefreq.time - median(tfr2.timefreq.time);
  tfr2.timefreq.power  = zeros(nchan,ntime,nfreq);
  tfr2(2) = tfr2;
  % average over channels
  for c = 1:length(clusters)
    % create cell array fo TFR matrices for this cluster type
    alldat = arrayfun(@(x)x.TFR2,clusters(c).epochs,'uniformoutput',false);
    tfr2(c).timefreq.event_code = c;
    for i  = 1:nchan
      % find clusters involving this channel
      cix  = find(arrayfun(@(x)(ismember(peaks(i).label,x.InvolvedChans)),clusters(c).epochs));
      if isempty(cix)
        tfr2(c).sensor_info(i).badchan = 1;
        continue;
      end
      % get indices to the matrix in those clusters for this channel
      ind  = arrayfun(@(x)find(ismember(x.InvolvedChans,peaks(i).label)),clusters(c).epochs(cix));
      % get data from those matrices for this channel
      tmp = cellfun(@(x,y)x(y,:,:),alldat(cix),num2cell(ind),'uniformoutput',false);
      % concatenate data along fourth dimension
      tmp = cat(4,tmp{:});
      % calculate power from complex spectra
      tmp = abs(double(tmp)).^2;
      % average POWER along fourth dimension
      tfr2(c).timefreq.power(i,:,:) = mean(tmp,4);
      trials{c}{i}                 = cix;
      clear tmp
    end
  end
  TFR2 = ts_combine_data(tfr2);
  clear alldat
toc 
TFR2 = ts_combine_data(tfr2);
tt   = TFR2.timefreq(1).time;
TFR2 = ts_data_selection(TFR2); % remove badchans
TFR2 = ts_combine_conditions(TFR2,'combinations',{'1-2'},'neweventcodes',3);
blcwin = 'all';%[tt(1) + [.3 .7]*(tt(end)-tt(1))]; 
zlim = [-3 3];
tfz2 = ts_zscore(TFR2,'baselinetype','zscore','blcwindow',blcwin,'verbose',0);
tfz2 = SO_freqband_average(tfz2,[8 12]);
ts_ezplot(TFR2,'events',[1 2],'blc',1,'baselinetype','zscore','blcwindow',blcwin,'zlim',zlim,'verbose',0,'layout',layout,'chantype','grad');
ts_ezplot(TFR2,'events',3,'blc',1,'baselinetype','zscore','blcwindow',blcwin,'zlim',zlim,'verbose',0,'layout',layout,'chantype','grad','showlabels','yes');
ts_ezplot(tfz2,'events',[1 2],'zlim',[-3 3],'layout',layout,'chantype','grad');
ts_ezplot(tfz2,'events',3,'zlim',[-3 3],'layout',layout,'chantype','grad');



% distribution around t=0
tlim  = [.3 .7]; 
flim  = [8 12];
blcwin = 'all';
pvals = {}; N = [];
N(1)  = length(clusters(1).epochs);
N(2)  = length(clusters(2).epochs);
tstart= tic;
for c = 1:length(clusters)
  pvals{c} = zeros(N(c),length(clusters(c).sensor_info));
  fprintf('cond %g of %g (%g)\n',c,length(clusters),toc(tstart));
  for k = 1:N(c)
    tmpTFR = ts_data_selection(TFR2,'chanlabel',clusters(c).epochs(k).InvolvedChans,'toilim',tlim,'foilim',flim,'events',c);
    tmpTFR.timefreq = tmpTFR.timefreq(1);
    tmpTFR.timefreq.power = abs(double(clusters(c).epochs(k).TFR2)).^2;
    tmpTFR = ts_zscore(tmpTFR,'baselinetype','zscore','blcwindow',blcwin,'verbose',0);
    [sel1,sel2] = match_str({clusters(c).sensor_info.label},{tmpTFR.sensor_info.label});
    pvals{c}(k,sel1) = mean(mean(tmpTFR.timefreq.power,3),2);
  end
end
fprintf('Time elapsed: %g min\n',toc(tstart)/60);

% somehow examine the difference in power values between the two conditions

% plots
peaktimes = zeros(2,nchan);
badchans  = find(flipvec==0);
for c = 1:2
  clear tmp
  [tmp{1:nchan,1}] = deal([]); allchanX{c} = tmp;
  [tmp{1:nchan,1}] = deal([]); allchanN{c} = tmp;
  % rearrange & get histogram values for each channel
  nbins     = 50;
  nchan     = size(pvals{1},2);
  chanpvals = zeros(nchan,N(c));
  chanX     = zeros(nchan,nbins);
  chanN     = zeros(nchan,nbins);
  tmpbadchan= zeros(nchan,1);
  for ch = 1:nchan
    if ismember(ch,badchans) || sum(pvals{c}(:,ch)~=0)<nbins, tmpbadchan(ch)=1; continue; end
    chanpvals(ch,:) = pvals{c}(:,ch)';%cat(1,pvals{c}(:,ch),pvals{2}(:,ch))';
    tmp             = chanpvals(ch,:);
    tmp             = tmp(tmp~=0);
    [tmpN,tmpX]     = hist(tmp,nbins);
    chanX(ch,1:length(tmpX)) = tmpX;
    chanN(ch,1:length(tmpX)) = tmpN;
    tmpix           = find(tmpN==max(tmpN));
    peaktimes(c,ch) = chanX(tmpix(1));
    allchanX{c}{ch} = tmpX;
    allchanN{c}{ch} = tmpN;
  end
%   tmpdat = ts_matrix2avg(chanN,'time',chanX(1,:));
%   tmpdat.sensor_info = flipdat.sensor_info;
%   [tmpdat.sensor_info(find(flipvec==0 | tmpbadchan==1)).badchan] = deal(1);
%   tmpdat = ts_data_selection(tmpdat);
%   % ts_ezplot(ts_data_selection(tmpdat),'layout',layout,'chantype','grad');
%   goodgrad1 = find(ismember({flipdat.sensor_info.typestring},'grad1')' & flipvec~=0 & tmpbadchan~=1);
%   goodgrad2 = find(ismember({flipdat.sensor_info.typestring},'grad2')' & flipvec~=0 & tmpbadchan~=1);
%   for j = 1:2
%     if j==1, goodchans=goodgrad1; else goodchans=goodgrad2; end
%     seldat = ts_data_selection(tmpdat,'chantype',sprintf('grad%g',j));
%     seldat.averages.event_code = c;
%     ts_ezplot(seldat,'layout',layout,'title',sprintf('histogram of mean power (%g-%gHz, %g-%gs wrt SOpeak)',flim,tlim),'linestyle','.');
%     hsubplots = get(gcf,'Children');
%     for k = 1:length(hsubplots)
%       hplot = hsubplots(k);
%       tmp   = get(hplot,'Children');
%       hdata = tmp(end);
%       set(hdata,'XData',chanX(goodchans(k),:));
%     end
%   end  
end

ch = 1; 
meanvals = zeros(nchan,2);
figure
for m = 1:8
  if m==1, pp=1; end
  for n = 1:8
    subplot(8,8,pp);
    x1=allchanX{1}{ch}; n1=allchanN{1}{ch};
    x2=allchanX{2}{ch}; n2=allchanN{2}{ch};
    scatter(x1(n1~=0),n1(n1~=0),'r.'); title(flipdat.sensor_info(ch).label); hold on
    scatter(x2(n2~=0),n2(n2~=0),'b.'); axis tight
    ch = ch + 1; pp = pp + 1;
    warning off
    meanvals(ch,1) = mean(n1(n1~=0));
    meanvals(ch,2) = mean(n2(n2~=0));
    warning on
  end
end

figure; plot(1:nchan,meanvals(1:nchan,1),'r.',1:nchan,meanvals(1:nchan,2),'b.');
figure; plot(1:nchan,meanvals(1:nchan,1)-meanvals(1:nchan,2),'.'); axis tight; hline(0,'k');
