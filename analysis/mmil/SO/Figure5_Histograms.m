tic
SubjID = 's1';
cd(sprintf('/space/mdkm1/2/kmdev/projects/jsherfey/sleep/meg/%s/matfiles',SubjID));
load proc_epoch_data_ICA; data = epoch_data; clear epoch_data
load SO_detections;                       alldet  = detections;
load SO_clustered_detections;             cstdet  = detections;
load SO_clustered_consistent_detections;  condet  = detections;
load SO_clusters_noflip;
clusters_Near_all = clusters_Near_all(1);
clusters_Nmax_all = clusters_Nmax_all(1);
load SO_clusters_consistency;
clusters_Cear_all = clusters_Cear_all(1);
clusters_Cmax_all = clusters_Cmax_all(1);

%% Histograms
figure('name','Histograms (maxR, all trials)');
% cluster size
Ccstsiz   = cellfun(@length,{clusters_Cmax_all.epochs.InvolvedChans});
Ncstsiz   = cellfun(@length,{clusters_Nmax_all.epochs.InvolvedChans});
[Nc,Xc]   = hist(Ccstsiz,25);
[Nn,Xn]   = hist(Ncstsiz,25);
subplot(2,3,1); plot(Xc,Nc,'b-',Xn,Nn,'r-'); legend('Consistent','No flip')
xlabel('cluster size (# involved channels)'); ylabel('count'); axis tight; vline(40,'k'); vline(20,'k');

% max delay
Cmaxdel   = cellfun(@max,{clusters_Cmax_all.epochs.Delays});
Nmaxdel   = cellfun(@max,{clusters_Nmax_all.epochs.Delays});
[Nc,Xc]   = hist(Cmaxdel,25);
[Nn,Xn]   = hist(Nmaxdel,25);
subplot(2,3,2); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('max delay (s)');

% mean delay
Cmeandel  = cellfun(@mean,{clusters_Cmax_all.epochs.Delays});
Nmeandel  = cellfun(@mean,{clusters_Nmax_all.epochs.Delays});
[Nc,Xc]   = hist(Cmeandel,25);
[Nn,Xn]   = hist(Nmeandel,25);
subplot(2,3,3); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('mean delay (s)');

% inter-detection-interval
Cidi  = diff([clusters_Cmax_all.epochs.RefTime]);
Nidi  = diff([clusters_Nmax_all.epochs.RefTime]);
[Nc,Xc]   = hist(Cidi(Cidi<200),200);
[Nn,Xn]   = hist(Nidi(Nidi<200),200);
subplot(2,3,4); plot(Xc,Nc,'b-',Xn,Nn,'r-');
xlabel('inter-detection interval (s)'); set(gca,'xlim',[0 20]);

% absmax amplitude
Cmaxamp   = cellfun(@(x)max(abs(x)),{clusters_Cmax_all.epochs.DetectionAmps});
Nmaxamp   = cellfun(@(x)max(abs(x)),{clusters_Nmax_all.epochs.DetectionAmps});
[Nc,Xc]   = hist(Cmaxamp,25);
[Nn,Xn]   = hist(Nmaxamp,25);
subplot(2,3,5); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('absmax amplitude (fT/cm2)');

% mean amplitude
Cmeanamp  = cellfun(@mean,{clusters_Cmax_all.epochs.DetectionAmps});
Nmeanamp  = cellfun(@mean,{clusters_Nmax_all.epochs.DetectionAmps});
[Nc,Xc]   = hist(Cmeanamp,25);
[Nn,Xn]   = hist(Nmeanamp,25);
subplot(2,3,6); plot(Xc,Nc,'b-',Xn,Nn,'r-'); axis tight;
xlabel('mean amplitude (fT/cm2)');


%% (% pos clustered detections) vs (% neg clustered detections)
nchan = length(alldet);
Nalln = cellfun(@length,{alldet.negpeak});
Nallp = cellfun(@length,{alldet.pospeak});
Ncstn = cellfun(@length,{cstdet.negpeak});
Ncstp = cellfun(@length,{cstdet.pospeak});
fracn = Ncstn ./ Nalln;
fracp = Ncstp ./ Nallp;
figure;
subplot(2,1,1),plot(fracp,fracn,'.'); axis([0 1 0 1]); lsline; axis square
ylabel('neg peaks'); title('% detections in clusters (no flipping)');
Ncstn = cellfun(@length,{condet.negpeak});
Ncstp = cellfun(@length,{condet.pospeak});
fracn = Ncstn ./ Nalln;
fracp = Ncstp ./ Nallp;
subplot(2,1,2),plot(fracp,fracn,'.'); axis([0 1 0 1]); lsline; axis square
xlabel('pos peaks'); ylabel('neg peaks'); title('% detections in clusters (consistency analysis)');
