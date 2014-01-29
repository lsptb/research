load s8_SO_cluster_results_Ref-MEG0143MEG0213MEG0723MEG1323MEG1423MEG2122-negpeak_filt0.01-4Hz_toi600-1350_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_18-Jun-2010.mat
% load s8_SO_cluster_results_Ref-MEG0143MEG0213MEG0723MEG1323MEG1423MEG2122-negpeak_filt0.01-4Hz_toi600-1350_grad1grad2_refall_smooth0.05sec_zero2zero0.25-1sec_mindist0.5sec_derivpks_pairedpeaks_18-Jun-2010.mat
avgdat = flip_dat;

% visualizer: avgdata => view cross-correlation
ref = 1;
avgdat.epochs(1).data = SOreference(ref).clusters(1).avgdata;
avgdat.epochs(1).time = SOreference(ref).clusters(1).epochtime;
avgdat.epochs(2) = avgdat.epochs(1);
avgdat.epochs(2).data = SOreference(ref).clusters(2).avgdata;
evnts = [1 2];
visualizer(ts_data_selection(avgdat,'events',evnts))

% look at paired peak power & compare to continuous filtered signal

ylim   = [0 1E-19];
foilim = 30 50];
for k = 1:length(SOreference)
  nr  = 10;
  nc  = 10;
  nneg= length(SOreference(k).clusters(1).PSD_mean);
  npos= length(SOreference(k).clusters(2).PSD_mean);
  nmax= min(nneg,npos);
  foi = SOreference(k).clusters(1).epochs(1).PSDfreq;
  fix = foi>=foilim(1) & foi<=foilim(2);
  P0  = zeros(size(SOreference(k).clusters(1).epochs(1).PSDfreq(fix)));
  nsmooth = round(length(fix)/10);
  if k == 1 && ~exist('hfig1','var'), hfig1 = figure; else figure(hfig1); end
  set(gcf,'Name',SOreference(k).reflabel);
  cnt = 0; Pm1=P0; Pm2=P0;
  for j = 1:nmax
    cnt = cnt + 1;
    f1  = SOreference(k).clusters(1).epochs(j).PSDfreq(fix);
    p1  = smooth(SOreference(k).clusters(1).epochs(j).PSD(fix),nsmooth)';
    f2  = SOreference(k).clusters(2).epochs(j).PSDfreq(fix);
    p2  = smooth(SOreference(k).clusters(2).epochs(j).PSD(fix),nsmooth)';
    subplot(nr,nc,cnt)
    plot(f1,p1,'b-',f2,p2,'r-'); axis tight, vline(30,'k');vline(50,'k');%set(gca,'ylim',ylim);
    ntk = t(SOreference(k).peaks(strmatch(SOreference(k).reflabel,{SOreference(k).peaks.label})).negpeak(j));
    ptk = t(SOreference(k).peaks(strmatch(SOreference(k).reflabel,{SOreference(k).peaks.label})).pospeak(j));
    title(sprintf('%g: t=%g/%g',j,ntk,ptk));
    if cnt==1,legend('neg','pos'); else axis off; end
    Pm1 = Pm1 + p1/(nr*nc);
    Pm2 = Pm2 + p2/(nr*nc);
    PSDmu1(cnt) = mean(p1);
    PSDmu2(cnt) = mean(p2);
    if mod(j+1,nr*nc)==1 && j+1<=nmax
      if k == 1 && ~exist('hfig2','var'), hfig2 = figure; else figure(hfig2); end
      subplot(4,1,1),plot(f1,Pm1,'b-o',f2,Pm2,'r-o'); axis tight, %set(gca,'ylim',ylim);
      title(sprintf('%s: mean power (n=%g)',SOreference(k).reflabel,nr*nc)); legend('neg','pos')
      epochix = j-nr*nc+1:j;
      nepoch1 = length(epochix);
      PSD1    = SOreference(k).clusters(1).PSD_mean(epochix);
      nepoch2 = length(epochix);
      PSD2    = SOreference(k).clusters(2).PSD_mean(epochix);
      if parms.Coh_flag
        Coh1    = SOreference(k).clusters(1).Coh_mean(epochix);
        Coh2    = SOreference(k).clusters(2).Coh_mean(epochix);          
        subplot(4,1,2),plot(1:nepoch1,Coh1,'b-o',1:nepoch2,Coh2,'r-o'); title(sprintf('Coherence (%g-%gHz)',parms.Coh_foilim)); axis tight, set(gca,'ylim',[0 .4]);
      end
      subplot(4,1,3),plot(1:nepoch1,PSD1,'b-o',1:nepoch2,PSD2,'r-o'); title(sprintf('Power (%g-%gHz)',parms.PSD_foilim)),xlabel('epoch #'); axis tight, %set(gca,'ylim',ylim);
      subplot(4,1,4),plot(1:nepoch1,PSDmu1,'b-o',1:nepoch2,PSDmu2,'r-o'); title(sprintf('Power (%g-%gHz)',foilim)),xlabel('epoch #'); axis tight, %set(gca,'ylim',ylim);
      pause
      cnt = 0; Pm1=P0; Pm2=P0;
      figure(hfig1)
    end      
  end
end

