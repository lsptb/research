tic
addpath /space/mdeh1/10/halgdev/projects/jsherfey/sleep/scripts
outpath = '/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/s8';

% SO detection parameters
parms     = [];
parms.fc1 = .1;   % Hz (lower cut-off freq prior to detection)
parms.fc2 = 2;    % Hz (upper cut-off freq prior to detection)
parms.monothresh = -5; % increase to be stricter; rec: stop at fc1
parms.minzero = 300; % ms (minimum distance between zero-crossings)
parms.amp_stdfact = 0; % # std > or < mean to threshold

parms.monophase_flag  = 1;
parms.surround_flag   = 1;
parms.interdet_flag   = 0;
parms.zero_flag       = 1;    
parms.zeroplus_flag   = 1; 
parms.amp_flag        = 1;    
parms.TFrej_flag      = 0;    
parms.preproc_flag    = 1;

% % LOOP OVER SUBJECTS
subjects  = {'s5','s6','s7','s8'};
subj=4;%for subj  = 4%:length(subjects)
  fprintf('Subject %g of [%s]\n',subj,num2str(1:length(subjects)));
  subject = subjects{subj};
% % subject matfiles
% if subj == 1
%   matfiles = {...

  %% LOAD DATA AND FIND PEAKS
  grads = {'grad1','grad2'}; graphmarker = 'o+';
  gradtype=1;%for gradtype = 1:length(grads)
    fprintf('Gradiometer type %g of %g\n',gradtype,length(grads));
    matfiles = {...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_1_nb01_060808_' grads{gradtype} '.mat'] ...
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_2_nb01_060808_' grads{gradtype} '.mat'] ...  
    ['/home/halgdev/projects/jsherfey/sleep/' subject '/matfiles/SL_3_nb01_060808_' grads{gradtype} '.mat'] ...
    };
    data = SO_combine_matfiles(matfiles);
    toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8
    data    = ts_data_selection(data,'toilim',toilim);
    fprintf('Time = %g to %g sec\n',toilim);

    proc = ts_preproc(data,'bpfilter','yes','bpfreq',[.1 8]);
    visualizer(proc)
    
    %% time-frequency analysis
    tic
    alpha   = 1;
    bpfreq  = [5 75];    
    toi1    = [1050 1175];
    toi2    = [1050 1175];
    labels  = {'MEG0113','MEG0122','MEG0132','MEG0143','MEG0213','MEG0222','MEG0313','MEG0322','MEG0343','MEG0413','MEG0422','MEG0432','MEG0443','MEG0723'};
    % constant spectral resolution
    foi     = 10:4:70;
    sf      = 2;
%     % variable spectral resolution
%     foi     = [2:2:100 105:5:300];
%     sf      = [2*ones(1,50) 5*ones(1,40)];    
    % select data to display with TF results
    dat = ts_data_selection(data,'chanlabel',labels,'toilim',toi2);
    % filter data for display
    dat = ts_preproc(dat,'bpfilter','yes','bpfreq',[.1 10]);
    % select channels of interest
    proc    = ts_data_selection(data,'chanlabel',labels);
    % filter out undesired freqs before wavelet
    proc    = ts_preproc(proc,'bpfilter','yes','bpfreq',bpfreq);
    % remove data padding used for filtering
    proc    = ts_data_selection(proc,'toilim',toi1);
    % Morlet wavelet analysis
    tfdata  = ts_freqanalysis_fieldtrip(proc,'foi',foi,'sf',sf,'trials_flag',1,'save_flag',0);
    % remove edge artifacts produced by wavelet analysis
    tfdata  = ts_data_selection(tfdata,'toilim',toi2);
    % calculate z-scores wrt entire recording
%     zdata   = ts_zscore(tfdata,'baselinetype','relchange','blcwindow','all','verbose',0);
    % correct for 1/f^2 scaling of power spectrum
    zdata = tfdata;
    freqs = tfdata.timefreq.frequencies;
    fcorr = repmat(permute(freqs.^(alpha),[1 3 2]),[zdata.num_sensors length(zdata.timefreq.time) 1]);
    zdata.timefreq.power = zdata.timefreq.power .* fcorr;
    clear fcorr
    % average over a frequency band
    freqband= [20 60];
    tfband  = SO_freqband_average(zdata,freqband);
    % normalize band-average and scale wrt the processed time series
    % note: this is necessary for overlaying band-avg & time-domain data
    for ch = 1:tfband.num_sensors
      % de-mean
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) - mean(tfband.averages.data(ch,:));
      % normalize
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) / max(abs(tfband.averages.data(ch,:)));
      % rescale
      tfband.averages.data(ch,:) = tfband.averages.data(ch,:) * max(abs(dat.epochs(1).data(ch,:)));
    end
    % add band-average to the time-domain data structure
    dat.epochs(2)      = dat.epochs(1);
    dat.epochs(2).data = tfband.averages.data;
    % change the event code for the tfband
    dat.epochs(2).event_code = 2;
    toc
    % visualize raw time series, band-average, and TFR power
    visualizer(dat,zdata);
    
    visualizer(ts_data_selection(data,'chanlabel',labels,'toilim',toi2),...
               ts_zscore(tfdata,'baselinetype','absolute','blcwindow','all','verbose',0));
    
    % multiplot power
    ts_ezplot(zdata,'zlim',[-3 3],'showlabels','yes');
    
    % overlay smoothed band-avgs for all channels
    tmp = tfband.averages.data;
    for k = 1:size(tmp,1)
      tmp(k,:) = smooth(tmp(k,:),10);
    end
    figure; plot(tfband.averages.time,tmp); axis tight

    % overlay raw data, band-avg, and filtered signal
    filtdat = dat;
    filtdat(2) = ts_data_selection(proc,'toilim',toi2);
    filtdat(2).epochs.event_code = 3;
    filtdat    = ts_combine_data(filtdat);
    visualizer(dat);
    
    % determine alpha for freq-corr if given TFR over long period
%     % note: this only works if tfdata is very long
%     ch   = 8;
%     fix  = 1:(length(freqs)-3);
%     x    = log10(freqs(fix));
%     y    = log10(squeeze(mean(tfdata.timefreq.power(ch,:,fix),2)));
%     P    = polyfit(x',y,1);
%     yfit = P(1)*x+P(2);
%     figure; plot(x,y,'-',x,yfit,'--')
%     alpha = -P(1)
    
%% UP/DOWN state detection (Mulkovski et al, Cerebral Cortex, 2007)

data    = SO_combine_matfiles(matfiles);
toilim  = [600 1350]; % 12.5 minutes during NREM 3/4 for s8
data    = ts_data_selection(data,'toilim',toilim);

labels  = {'MEG0113','MEG0143','MEG0213','MEG0222','MEG0322','MEG0343'};
% bandpass filter 20-100Hz
dat     = ts_data_selection(data,'toilim',[1000 1300],'chanlabel',labels);
proc    = ts_preproc(dat,'bpfilter','yes','bpfreq',[20 100]);
dat     = ts_data_selection(dat,'toilim',[1025 1250]);
proc    = ts_data_selection(proc,'toilim',[1025 1250]);

% calculate standard deviation in a running 5- or 10-ms window (=RMS)
StdDevWindow = 10; % ms
SmoothWindow = 50; % ms

Fs = proc.sfreq;
dt = 1/Fs;

% convert windows from ms to indices and force odd
StdDevWindow = round(StdDevWindow/(1000*dt));
SmoothWindow = round(SmoothWindow/(1000*dt));
if ~mod(StdDevWindow,2),StdDevWindow=StdDevWindow+1; end
if ~mod(SmoothWindow,2),SmoothWindow=SmoothWindow+1; end
begntime  = length(proc.epochs.time);
endntime  = begntime - 2*floor(StdDevWindow/2);
beg_ind   = floor(StdDevWindow/2)+1;
end_ind   = begntime - floor(StdDevWindow/2);

tic
% initialize std struct and calculate std in moving window
stdproc = proc;
stdproc.epochs.time = proc.epochs.time(beg_ind:end_ind);
stdproc.epochs.data = zeros(proc.num_sensors,endntime);
for k = 1:endntime
  stdproc.epochs.data(:,k) = std(proc.epochs.data(:,k:k+StdDevWindow-1),0,2);
end
toc
% % very slow
% lo_ind = num2cell(1:end_ind);
% hi_ind = num2cell(beg_ind:begntime);
% stdproc.epochs.data(1,:) = cellfun(@(a,b)std(proc.epochs.data(1,a:b)),lo_ind,hi_ind)

% linear smoothing with a 50-ms window
METHOD = 'lowess';
SPAN   = SmoothWindow;
smoothstd = stdproc;
smoothstd.epochs.data = zeros(size(stdproc.epochs.data));
tic
for k = 1:stdproc.num_sensors
  smoothstd.epochs.data(k,:) = smooth(stdproc.epochs.data(k,:),SPAN,METHOD);
  toc
end


%%    
% % %     
% % %     
% % % 
% % %     reflabels   = {data.sensor_info.label}; %{'MEG0633','MEG1043','MEG0713','MEG0723','MEG0743','MEG0733'};
% % %     nref        = length(reflabels);
% % %     nchan       = data.num_sensors;
% % %     fprintf('Finding slow oscillations in %g reference channels\n',nref);
% % %     [chans,jnk] = match_str({data.sensor_info.label},reflabels);
% % %     refdata     = ts_data_selection(data,'channels',chans);
% % %     peaks   = SO_peaks(refdata,parms);
% % %     t       = refdata.epochs.time;
% % %     events  = [];
% % %     for k   = 1:length(peaks)
% % %       events(k).label = peaks(k).label;
% % %       events(k).time  = [t([peaks(k).pospeak peaks(k).negpeak])];
% % %       events(k).type  = [1*ones(1,length(peaks(k).pospeak)) 2*ones(1,length(peaks(k).negpeak))];
% % %     end
% % %     outfile = sprintf('%s/s8_SOpeaks_filt%g-%gHz_toi%g-%g_%s_ref%s.mat',outpath,parms.fc1,parms.fc2,toilim,grads{gradtype},'all');%[reflabels{:}]);
% % %     save(outfile,'events','peaks');
% % %     toc
% % %     data = ts_preproc(data,'bpfilter','yes','bpfreq',[parms.fc1 parms.fc2]);
% % %     toc
% % %     
% % %     
% % %     
% % %     
% % %     
% % %     %% EPOCH CONTINUOUS DATA BASED ON PEAKS IN REFERENCE CHANNEL
% % %     % select the peak indices that fall within the selected data range
% % % %     [sel1,sel2] = match_str({peaks.label},{data.sensor_info.label});
% % % %     peaks = peaks(sel1);
% % % %     for k = 1:length(peaks)
% % % %       t = peaks(k).tstart:1/peaks(k).sfreq:peaks(k).tstop;
% % % %       peaks(k).pospeak = peaks(k).pospeak(t(peaks(k).pospeak)>=toilim(1) & t(peaks(k).pospeak)<=toilim(2));
% % % %       peaks(k).negpeak = peaks(k).negpeak(t(peaks(k).negpeak)>=toilim(1) & t(peaks(k).negpeak)<=toilim(2));
% % % %       peaks(k).tstart  = toilim(1);
% % % %       peaks(k).tstop   = toilim(2);
% % % %     end
% % % %     tic
% % % %     % refchan-based epoching
% % % %     for refchan = 1:nref
% % % %       fprintf('Epoching %g gradiometers using reference %g of %g\n',nchan,refchan,nref);
% % % %       reflabel = reflabels{refchan};
% % % %       % get the peaks that will be used to epoch every channel
% % % %       clear pospeaks negpeaks sodata
% % % %       pos  = peaks(refchan).pospeak;
% % % %       neg  = peaks(refchan).negpeak;
% % % %       % copy the reference channel peaks for the other channels
% % % %       [pospeaks{1:nchan}] = deal(pos); % use pospeaks in reference chan for all chans
% % % %       [negpeaks{1:nchan}] = deal(neg); % use negpeaks in reference chan for all chans
% % % %       % epoching parameters
% % % %       epochpad = 1000; % ms
% % % %       % pospeaks - data
% % % %       npeaks  = cellfun(@length,pospeaks);
% % % %       posdata = SO_epochs(data,pospeaks,epochpad);
% % % %       % negpeaks - data
% % % %       npeaks  = cellfun(@length,negpeaks);
% % % %       negdata = SO_epochs(data,negpeaks,epochpad);
% % % %       save(fullfile(outpath,sprintf('%s_SO-refchan-%s_toi%g-%gsec_pos_epoch_data_%s.mat',subject,reflabel,toilim,grads{gradtype})),'posdata','pospeaks');
% % % %       save(fullfile(outpath,sprintf('%s_SO-refchan-%s_toi%g-%gsec_neg_epoch_data_%s.mat',subject,reflabel,toilim,grads{gradtype})),'negdata','negpeaks');
% % % %       clear posdata negdata pospeaks negpeaks
% % % %       toc
% % % %     end
% % %     %% CALCULATE LATENCIES BETWEEN PEAKS IN DIFFERENCE CHANNELS
% % %     load(outfile)
% % %     testflag     = 1;
% % %     plotflag     = 1;
% % %     flipflag     = 0; % 0-no, 1-avg, 2-trial
% % %     load s8_ReferenceSO_600-1350sec_allrefchan_grad1.mat
% % %     theta        = ts_BetweenSensorAngles(ReferenceSO(1).sensor_info);
% % %     threshfactor = .5;
% % %     reflabels    = {events.label};
% % %     nref         = length(reflabels);
% % %     nchan        = length(events);%data.num_sensors;
% % %     clear data
% % %     % loop over reference channels
% % %     ReferenceSO = [];
% % %     tic
% % %     for refchan = 3%1:50%length(reflabels)
% % %       reflabel = reflabels{refchan};
% % %       fprintf('Calculating SO latencies wrt reference %g of %g (%s)\n',refchan,length(reflabels),reflabel);
% % %       infile1  = fullfile(outpath,sprintf('%s_SO-refchan-%s_toi%g-%gsec_neg_epoch_data_%s.mat',subject,reflabel,toilim,grads{gradtype})); % negdata
% % %       infile2  = fullfile(outpath,sprintf('%s_SO-refchan-%s_toi%g-%gsec_pos_epoch_data_%s.mat',subject,reflabel,toilim,grads{gradtype})); % posdata
% % %       load(infile1); if flipflag==1, negavg = SO_average(negdata,length(peaks(refchan).negpeak)*ones(1,negdata.num_sensors)); end
% % %       load(infile2); if flipflag==1, posavg = SO_average(posdata,length(peaks(refchan).pospeak)*ones(1,posdata.num_sensors)); end
% % %       % combine negative and positive SOs
% % %       thisdata = negdata;  % thisdata = posdata;
% % %       thisdata.epochs.num_trials(2) = posdata.epochs.num_trials;
% % %       thisdata.epochs.data          = cat(3,thisdata.epochs.data,posdata.epochs.data);
% % %       thisdata.epochs.matfiles      = {thisdata.epochs.matfiles{:} infile1 infile2};
% % %       refindex = strmatch(reflabel,{thisdata.sensor_info.label},'exact');  
% % %       ntrial   = sum(thisdata.epochs.num_trials);
% % %       ReferenceSO(refchan).reflabel        = reflabel;
% % %       ReferenceSO(refchan).sensor_info     = thisdata.sensor_info;    
% % %       ReferenceSO(refchan).rawmatfiles     = thisdata.epochs.matfiles;
% % %       ReferenceSO(refchan).epochmatfiles   = {infile1 infile2};
% % %       ReferenceSO(refchan).peaktypes       = {'negpeak','pospeak'};
% % %       ReferenceSO(refchan).num_trials      = [negdata.epochs.num_trials posdata.epochs.num_trials];
% % %       clear posdata negdata
% % %       if plotflag
% % %         nrow = ceil(sqrt(ntrial));
% % %         ncol = ceil(sqrt(ntrial));
% % %         cnt  = 0;
% % %         figure(refchan); set(gcf,'Name',sprintf('%s (%g of %g)',reflabel,refchan,nref));
% % %       end      
% % %       % loop of SO cycles
% % %       for trl  = 1:ntrial
% % %         negreference = trl <= ReferenceSO(refchan).num_trials(1);
% % % %         fprintf('trial %g of %g\n',trl,ntrial);
% % %         % get amplitude of peak in refchan
% % %         refA   = thisdata.epochs.data(refindex,nearest(thisdata.epochs.time,0),trl);
% % %         % get amplitude threshold based on this cycle in the refchan
% % %         threshold = threshfactor*abs(refA);%max(abs([thisdata.epochs.data(refindex,:,trl)]));
% % %         % preallocation
% % %         keep = zeros(1,nchan);
% % %         flip = zeros(1,nchan);
% % %         pkA  = zeros(1,nchan);
% % %         pkt  = zeros(1,nchan);
% % %         grapholor = repmat('b',1,nchan);
% % %         ReferenceSO(refchan).cycle(trl).reftime   = 0;
% % %         ReferenceSO(refchan).cycle(trl).refpeak   = 0;
% % %         ReferenceSO(refchan).cycle(trl).keepchan  = keep;
% % %         ReferenceSO(refchan).cycle(trl).flipchan  = flip;
% % %         ReferenceSO(refchan).cycle(trl).amplitude = pkA; 
% % %         ReferenceSO(refchan).cycle(trl).latency   = pkt;
% % %         if refchan == 1 && trl == 1, [ReferenceSO(2:nref)] = ReferenceSO; end
% % %         % loop over channels
% % %         for ch = 1:nchan
% % %           if trl==1, ReferenceSO(refchan).sensor_info(ch).angle = theta(refindex,ch); end
% % %           t  = thisdata.epochs.time;
% % %           y  = thisdata.epochs.data(ch,:,trl);
% % %           if flipflag == 1
% % %             % flip based on polarity of averages
% % %             if negreference
% % %               % negpeak refchan => flip if peak(t=0/max) > 0
% % %               if negavg.averages.data(ch,nearest(negavg.averages.time,0)) > 0
% % %                 y = -y;
% % %                 thisdata.epochs.data(ch,:,trl) = y;
% % %                 flip(ch) = 1;
% % %               end
% % %             else
% % %               % pospeak refchan => flip if peak(t=0/max) < 0
% % %               if posavg.averages.data(ch,nearest(posavg.averages.time,0)) < 0
% % %                 y = -y;
% % %                 thisdata.epochs.data(ch,:,trl) = y;
% % %                 flip(ch) = 1;
% % %               end
% % %             end  
% % %             % we are only interested in peaks with same polarity now
% % %           end
% % %           dt = diff(t);
% % %           dy = diff(y);
% % %           m  = dy./dt;
% % %           I1 = (m(1:end-1) < 0) & (m(2:end) > 0); % neg to pos slope
% % %           I2 = (m(1:end-1) > 0) & (m(2:end) < 0); % pos to neg slope
% % %           I  = find(I1 | I2); % peak times
% % %           
% % %           % eliminate peaks of opposite polarity
% % %           if negreference
% % %             % eliminate positive peaks
% % %             I(y(I)>0) = [];
% % %           else
% % %             % eliminate negative peaks
% % %             I(y(I)<0) = [];
% % %           end
% % %           if isempty(I)
% % %             A        = 0;
% % %             pkt(ch)  = 0;
% % %             pkA(ch)  = 0;
% % %             keep(ch) = 0;
% % %           else
% % %             % select remaining peak that is closest to t=0
% % %             I0  = I(nearest(t(I),0)); % index to peak closest to t=0
% % %             A   = y(I0); % amplitude of peak closest to t=0
% % %             pkt(ch) = t(I0);
% % %             pkA(ch) = A;
% % % 
% % %             % amplitude threshold
% % %             if abs(A) > threshold
% % %               keep(ch) = 1;
% % %             end
% % %           end
% % %             % testing
% % %             if testflag
% % %                 plot(t,y,'b',t(I),y(I),'k.',t(I0),A,'go','MarkerSize',8);
% % %                 ylims = max([abs(y) abs(refA)]);
% % %                 ylims = 1.1*[-ylims ylims];
% % %                 axis([-1 1 ylims])
% % %                 hline(0,'k'); vline(0,'k'); 
% % %                 h=refline(0,refA); set(h,'Color','b');
% % %                 h=refline(0,threshold); set(h,'Color','r','LineStyle','--');
% % %                 h=refline(0,-threshold); set(h,'Color','r','LineStyle','--');
% % %                 title(sprintf('%g/%s',ch,thisdata.sensor_info(ch).label)); 
% % %                 if abs(A)<threshold
% % %                   ylims = get(gca,'ylim');
% % %                   xlims = get(gca,'xlim');
% % %                   h = refline(diff(ylims)/diff(xlims)); set(h,'Color','r');
% % %                   h = refline(-diff(ylims)/diff(xlims)); set(h,'Color','r');
% % %                 else
% % %                   if t(I0) < 0, xpts=[t(I0) 0]; else xpts=[0 t(I0)]; end
% % %                   ypts = [A A]; h=line(xpts,ypts); set(h,'LineWidth',2);
% % %                   ypts = [0 0]; h=line(xpts,ypts); set(h,'LineWidth',2);
% % %                   if A<0, ypts=[A 0]; else ypts=[0 A]; end
% % %                   xpts = [t(I0) t(I0)]; h=line(xpts,ypts); set(h,'LineWidth',2);
% % %                   xpts = [0 0];         h=line(xpts,ypts); set(h,'LineWidth',2);
% % %                   if (refA < 0 && A > 0) || (refA > 0 && A < 0)
% % %                     title(sprintf('%g/%s (FLIP!)',ch,thisdata.sensor_info(ch).label)); 
% % %                   end
% % %                 end
% % %                 pause
% % %             end
% % %           if flipflag == 2
% % %             % find peak closest to t=0
% % %             I0 = I(nearest(t(I),0)); % index to peak closest to t=0
% % %             A  = y(I0); % amplitude of peak closest to t=0
% % %             
% % %             % keep if abs(peak) @ I0 > threshold
% % %             if abs(A) > threshold
% % %               keep(ch) = 1;
% % %               pkt(ch)  = t(I0);
% % %               % graphcolor = red if peak before ref; else blue
% % %               if t(I0) < 0, graphcolor(ch) = 'r'; end
% % %               % flip if peak polarity is opposite of the reference
% % %               if (refA < 0 && A > 0) || (refA > 0 && A < 0)
% % %                 flip(ch) = 1;
% % %                 y = -y;
% % %                 A = -A;
% % %               end
% % %             end
% % %             thisdata.epochs.data(ch,:,trl) = y;
% % %             pkA(ch) = A;
% % %           elseif flipflag == 1 % flip with average
% % %             
% % %           elseif flipflag == 0 % do not flip
% % %             
% % %           else
% % %             error('Illegal flipflag');
% % %           end
% % %         end % end loop over channels
% % %         % record info
% % %         ReferenceSO(refchan).cycle(trl).reftime   = 0;
% % %         ReferenceSO(refchan).cycle(trl).refpeak   = refA;
% % %         ReferenceSO(refchan).cycle(trl).keepchan  = keep;
% % %         ReferenceSO(refchan).cycle(trl).flipchan  = flip;
% % %         ReferenceSO(refchan).cycle(trl).latency   = pkt;
% % %         ReferenceSO(refchan).cycle(trl).amplitude = pkA;
% % %         % calculate correlations and such
% % %         angle = [ReferenceSO(refchan).sensor_info(keep==1).angle]';
% % %         delay = pkt(keep==1)';
% % %         % sort by angular distance
% % %         [angle,I] = sort(angle);
% % %         delay = delay(I);
% % %         % remove nans
% % %         rmix = isnan(angle);
% % %         angle(rmix) = [];
% % %         delay(rmix) = [];
% % %         
% % %         pos_ind       = delay > 0;
% % %         pos_skipflag  = 0;
% % %         if sum(pos_ind) < 2, pos_skipflag = 1; end
% % %         neg_ind       = delay < 0;
% % %         neg_skipflag  = 0;
% % %         if sum(neg_ind) < 2, neg_skipflag = 1; end
% % %         
% % %         if ~pos_skipflag
% % %           [c,r]         = polyfit(angle(pos_ind),delay(pos_ind),1);
% % %           tfit_posdelay = c(1)*angle(pos_ind) + c(2);
% % %           Err_posdelay  = sqrt(sum(abs(tfit_posdelay-delay(pos_ind))).^2/length(tfit_posdelay));        
% % %           [R,pval]      = corrcoef(angle(pos_ind),delay(pos_ind));
% % %           R_pos         = R(1,2);
% % %           pval_pos      = pval(1,2);
% % %         end
% % %         if ~neg_skipflag
% % %           [c,r] = polyfit(angle(neg_ind),delay(neg_ind),1);
% % %           tfit_negdelay = c(1)*angle(neg_ind) + c(2);
% % %           Err_negdelay  = sqrt(sum(abs(tfit_negdelay-delay(neg_ind))).^2/length(tfit_negdelay));
% % %           [R,pval]      = corrcoef(angle(neg_ind),delay(neg_ind));
% % %           R_neg    = R(1,2);
% % %           pval_neg = pval(1,2);    
% % %         end
% % %         
% % %         if ~(pos_skipflag && neg_skipflag)
% % %           [R,pval] = corrcoef(angle,delay);
% % %           R        = R(1,2);
% % %           pval     = pval(1,2);
% % %         end
% % %         % plot
% % %         if plotflag
% % %           cnt = cnt + 1;
% % %           subplot(nrow,ncol,cnt), 
% % %           if ~pos_skipflag, plot(angle(pos_ind),delay(pos_ind),['b' graphmarker(gradtype)],angle(pos_ind),tfit_posdelay,'k-'); hold on; end   
% % %           if ~neg_skipflag, plot(angle(neg_ind),delay(neg_ind),['r' graphmarker(gradtype)],angle(neg_ind),tfit_negdelay,'k-'); hold on; end
% % %           axis tight; set(gca,'ylim',[-1 1]); set(gca,'xlim',[0 200]);%[min(distances) max(distances)])
% % %           title(sprintf('(%g) R = %g, %g',trl,(R_pos+R_neg)/2,(pval_pos<.05) && (pval_neg<.05)));
% % %           vline(0,'k'); hline(0,'k');  %axis off
% % %         end
% % %         % record corr info
% % %         if ~(pos_skipflag && neg_skipflag), ReferenceSO(refchan).cycle(trl).corrcoef = R; end
% % %         if ~pos_skipflag, ReferenceSO(refchan).cycle(trl).posdelay_corrcoef = R_pos; end
% % %         ReferenceSO(refchan).cycle(trl).N_pos    = sum(pos_ind);
% % %         if ~neg_skipflag, ReferenceSO(refchan).cycle(trl).negdelay_corrcoef = R_neg; end
% % %         ReferenceSO(refchan).cycle(trl).N_neg    = sum(neg_ind);
% % %         ReferenceSO(refchan).cycle(trl).angle    = angle;
% % %         ReferenceSO(refchan).cycle(trl).delay    = delay;
% % % 
% % %       end % end loop over trials
% % %       if plotflag
% % %         xlabel('angular distance (deg)')
% % %         ylabel('latency (sec)')  
% % %         pause
% % %       end
% % %       toc
% % %     end % end loop over reference channels
% % % %     outfile = sprintf('%s_ReferenceSO_%g-%gsec_%s_%s.mat',subject,toilim,'allrefchan1-50',grads{gradtype});%[reflabels{:}]);
% % % %     save(fullfile(outpath,outfile),'ReferenceSO');
% % % 
% % %     %% PLOT RESULTS
% % %     plotflag     = 1;
% % %     % loop over reference channels
% % %     if plotflag
% % %       for refchan = 3%1:nref
% % %         ntrial = 49;%length(ReferenceSO(refchan).cycle);
% % %         for page = 1:ceil(length(ReferenceSO(refchan).cycle)/ntrial)
% % %             nrow = ceil(sqrt(ntrial));
% % %             ncol = ceil(sqrt(ntrial));
% % %             cnt  = 0;
% % %             reflabel = reflabels{refchan};
% % %             figure(refchan); set(gcf,'Name',sprintf('%s (%g of %g)',reflabel,refchan,nref));
% % %           for trl  = [1:ntrial] + ntrial*(page-1)
% % %             fprintf('%g of %g\n',trl,ntrial)
% % %             angle = ReferenceSO(refchan).cycle(trl).angle;
% % %             delay = ReferenceSO(refchan).cycle(trl).delay;
% % %             R     = ReferenceSO(refchan).cycle(trl).corrcoef;
% % %             R_pos = ReferenceSO(refchan).cycle(trl).posdelay_corrcoef;
% % %             R_neg = ReferenceSO(refchan).cycle(trl).negdelay_corrcoef;
% % %             pos_ind = delay > 0;
% % %             neg_ind = delay < 0;
% % %             % calculate best fit line and least-square error
% % %             [c,r] = polyfit(angle(pos_ind),delay(pos_ind),1);
% % %             tfit_posdelay = c(1)*angle(pos_ind) + c(2);
% % %             Err_posdelay  = sqrt(sum(abs(tfit_posdelay-delay(pos_ind))).^2/length(tfit_posdelay));
% % %             [c,r] = polyfit(angle(neg_ind),delay(neg_ind),1);
% % %             tfit_negdelay = c(1)*angle(neg_ind) + c(2);
% % %             Err_negdelay  = sqrt(sum(abs(tfit_negdelay-delay(neg_ind))).^2/length(tfit_negdelay));
% % %             % plot
% % %               cnt = cnt + 1;
% % %               subplot(nrow,ncol,cnt), 
% % %               plot(angle(pos_ind),delay(pos_ind),'b.',angle(pos_ind),tfit_posdelay,'k-'); hold on;       
% % %               plot(angle(neg_ind),delay(neg_ind),'r.',angle(neg_ind),tfit_negdelay,'k-'); hold on;       
% % %               axis tight; set(gca,'ylim',[-.5 .5]); set(gca,'xlim',[0 200]);%[min(distances) max(distances)])
% % %               title(sprintf('R:%g/%g',R_pos,R_neg));%title(sprintf('(%g) R = %g',trl,(R_pos+R_neg)/2));
% % %               vline(0,'k'); hline(0,'k');  %axis off
% % %           end % end loop over trials
% % %             xlabel('angular distance (deg)')
% % %             ylabel('latency (sec)')  
% % %             pause
% % %         end
% % %       end % end loop over reference channels
% % %     end
% % %     plotflag = 1;
% % %     if plotflag
% % %       Nref = length(ReferenceSO);
% % %       nrow = ceil(sqrt(Nref));
% % %       ncol = ceil(sqrt(Nref));
% % %       Nthresh = 10; % only plot R if N > Nthresh = # of channels in calc
% % %       figure(refchan)
% % %       for refchan = 1:Nref
% % %         subplot(nrow,ncol,refchan)
% % %         Npos    = [ReferenceSO(refchan).cycle.N_pos];
% % %         Npos(arrayfun(@(x)isempty(x.posdelay_corrcoef),ReferenceSO(refchan).cycle)) = [];
% % %         posplot = [ReferenceSO(refchan).cycle.posdelay_corrcoef];
% % %         posplot = posplot(Npos > Nthresh);
% % %         Nneg    = [ReferenceSO(refchan).cycle.N_neg];
% % %         Nneg(arrayfun(@(x)isempty(x.negdelay_corrcoef),ReferenceSO(refchan).cycle)) = [];
% % %         negplot = [ReferenceSO(refchan).cycle.negdelay_corrcoef];
% % %         negplot = negplot(Nneg > Nthresh);
% % %         plot(posplot,'b.'); axis tight; set(gca,'ylim',[-1 1]); hold on; hline(0,'k'); 
% % %         plot(negplot,'r.')
% % %         hline(.3,'Color','k','LineWidth',2); 
% % %         hline(-.3,'Color','k','LineWidth',2); 
% % %         if refchan==1
% % %           title(sprintf('%s/%s,ref=%s',subject,grads{gradtype},reflabels{refchan}))
% % %         else
% % %           title(sprintf('ref=%s',reflabels{refchan}))
% % %         end
% % %       signum(refchan) = sum(negplot>.3)+sum(negplot<-.3)+sum(posplot>.3)+sum(posplot<-.3);
% % %       totnum(refchan) = length(negplot) + length(posplot);        
% % %       end
% % %       xlabel('SO cycle')
% % %       ylabel('R(distance,delay)')
% % %     end
% % %     
% % %     plotflag = 1;
% % %     if plotflag
% % %       percent  = 100*(signum./totnum);
% % %       pthresh = 50;
% % %       refplots = reflabels(percent>pthresh);
% % %       refind   = find(percent>pthresh);
% % %       refrates = []; refcnt = 1;
% % %       avgnegrefdelay = []; avgposrefdelay = [];
% % %       for refchan = refind
% % %         fprintf('ref %g of %g\n',refcnt,length(refind))
% % %           ntrial = length(ReferenceSO(refchan).cycle);
% % %           showtrial = zeros(1,ntrial);
% % % %           for trl = 1:ntrial
% % % %             if ~isempty(ReferenceSO(refchan).cycle(trl).posdelay_corrcoef) && ReferenceSO(refchan).cycle(trl).posdelay_corrcoef>.3,showtrial(trl)=1; end
% % % %             if ~isempty(ReferenceSO(refchan).cycle(trl).negdelay_corrcoef) && ReferenceSO(refchan).cycle(trl).negdelay_corrcoef<-.3,showtrial(trl)=1; end
% % % %           end
% % % %           showtrial([2 7 13 26 36 39 49 53 109 147 155 159 194 200 243 249 294 303 330 375 387 404 433 438 455 465])=1;
% % %           showtrial(1:end)=1;
% % %           nshow = sum(showtrial);
% % %             nrow = ceil(sqrt(nshow));
% % %             ncol = ceil(sqrt(nshow));
% % %             cnt  = 0;
% % % %             figure(refchan); set(gcf,'Name',sprintf('%s (%g of %g)',reflabel,refchan,nref));
% % %           reflabel = reflabels{refchan};
% % %           rates = []; avgposdelays = []; avgnegdelays=[];
% % %           for trl  = 1:ntrial
% % %             if ~showtrial(trl), continue; end
% % %             angle = ReferenceSO(refchan).cycle(trl).angle;
% % %             delay = ReferenceSO(refchan).cycle(trl).delay;
% % %             R     = ReferenceSO(refchan).cycle(trl).corrcoef;
% % %             R_pos = ReferenceSO(refchan).cycle(trl).posdelay_corrcoef;
% % %             R_neg = ReferenceSO(refchan).cycle(trl).negdelay_corrcoef;
% % %             pos_ind = delay > 0;
% % %             neg_ind = delay < 0;
% % %             % calculate best fit line and least-square error
% % %             [c,r] = polyfit(angle(pos_ind),delay(pos_ind),1);
% % %             tfit_posdelay = c(1)*angle(pos_ind) + c(2);
% % %             [c,r] = polyfit(angle(neg_ind),delay(neg_ind),1);
% % %             tfit_negdelay = c(1)*angle(neg_ind) + c(2);
% % %             % plot
% % % %               if (isempty(R_pos) || R_pos<.3) && (isempty(R_neg) || R_neg>-.3), continue; end
% % % %               cnt = cnt + 1;
% % % %               subplot(nrow,ncol,cnt), 
% % % %               if R_pos > .3, plot(angle(pos_ind),delay(pos_ind),'b.',angle(pos_ind),tfit_posdelay,'k-'); hold on; end
% % % %               if R_neg < -.3, plot(angle(neg_ind),delay(neg_ind),'r.',angle(neg_ind),tfit_negdelay,'k-'); hold on; end
% % % %               axis tight; set(gca,'ylim',[-.5 .5]); set(gca,'xlim',[0 200]);%[min(distances) max(distances)])
% % % %               title(num2str(trl));%title(sprintf('(%g) R = %g',trl,(R_pos+R_neg)/2));
% % % %               vline(0,'k'); hline(0,'k');  hold off; %axis off
% % % %               drawnow
% % % %               pause
% % %               rates(trl) = mean([mean(delay(pos_ind)./angle(pos_ind)),mean(-delay(neg_ind)./angle(neg_ind))]);
% % %               avgposdelays(trl) = mean(delay(pos_ind));
% % %               avgnegdelays(trl) = mean(delay(neg_ind));
% % %           end % end loop over trials
% % %           refrates(refcnt) = mean(rates); 
% % %           avgposrefdelay(refcnt) = mean(avgposdelays);
% % %           avgnegrefdelay(refcnt) = mean(avgnegdelays);
% % %           refcnt = refcnt + 1;
% % % %           xlabel('angular distance (deg)')
% % % %           ylabel('latency (sec)')  
% % % %           pause
% % %       end % end loop over reference channels      
% % %     end
% % %     [avgposrefdelay' avgnegrefdelay']
% % %     refrates'
% % %     mean([avgposrefdelay(~isnan(avgposrefdelay)) -avgnegrefdelay(~isnan(avgnegrefdelay))]) % sec
% % %     reflabels(refind)
% % %     
% % % %   end % end loop over gradiometer types
% % % % end % end loop over subjects
% % % 
% % % 
% % % % Correct delays for refchan offset
% % % % delay = ReferenceSO(refchan).cycle(trl).delay; 
% % % % delay = delay - delay(refindex);
% % % % 
% % % % nkeepchan = arrayfun(@(x)length(find(x.keepchan)),ReferenceSO(1).cycle);
% % % % nposdelay = arrayfun(@(x)length(find(x.latency(x.keepchan==1)>0)),ReferenceSO(1).cycle);
% % % % nposdelay = arrayfun(@(x)length(find(x.latency(x.keepchan==1)<0)),ReferenceSO(1).cycle);
% % % % arrayfun(@(x)length(find(x.latency(x.keepchan==1)==0)),ReferenceSO(1).cycle)

