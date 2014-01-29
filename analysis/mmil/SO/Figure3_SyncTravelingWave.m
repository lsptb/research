% Figure 3 (21-Aug-2010, JSS)
% Last significant modification: 18-Oct-2010, JSS
%     - epoch_data was recalc with all files and "grad" was added to filename
cd ~jsherfey/projects/sleep3/
load s1/matfiles/proc_grad_epoch_data_ICA.mat
load s1/matfiles/SO_clusters_consistency.mat
params    = SO_params(1);

synctrial = 6498;
travtrial = 3096;
clim      = [0 .2];

offset = 2408.52; % 0 for original data; 2408.5 for new

% sync example
cst = clusters_Cmax_all(1).epochs(synctrial);
sync = ts_data_selection(epoch_data,'toilim',cst.HistTime+[-5 5],'chanlabel',cst.InvolvedChans);
sync = ts_preproc(sync,'bpfilter','yes','bpfreq',[.01 4]);

tpad = 1; badchans = [43 50];
sync = ts_data_selection(sync,'toilim',cst.HistTime+[-tpad tpad],'badchans',badchans);
SYNC = sync;
% for k = 1:sync.num_sensors
%   SYNC.epochs.data(k,:) = SYNC.epochs.data(k,:) / max(abs(SYNC.epochs.data(k,round(length(SYNC.epochs.time)/2)+[-20:20])));
% end
figure('color','w')
subplot(3,4,2); imagesc(SYNC.epochs.time,1:SYNC.num_sensors,SYNC.epochs.data);
set(gca,'clim',[-.05E-7 .05E-7]); set(gca,'xlim',cst.HistTime+[-tpad tpad]);
title('synchronous wave');ylabel('channel');
subplot(3,4,6); plot(sync.epochs.time,SYNC.epochs.data); axis tight
set(gca,'xlim',cst.HistTime+[-tpad tpad]);
xlabel('time (s)'); ylabel('amp (fT/cm)'); set(gca,'ylim',[-.13E-7 .13E-7])

refix  = round(strmatch(cst.RefChan,{epoch_data.sensor_info.label})/2);
sel    = setdiff(1:length(cst.InvolvedChans),badchans);
ch     = match_str({epoch_data.sensor_info.label},cst.InvolvedChans(sel));
deldat = epoch_data;
deldat.epochs.data = nan(deldat.num_sensors,1);
deldat.epochs.data(ch,:) = cst.Delays(sel);
deldat.epochs.time = 0;
deldat.averages = deldat.epochs;
deldat = rmfield(deldat,'epochs');
tmpdat = ts_data_selection(deldat,'chantype','grad1'); x1 = tmpdat.averages.data;
[highlight,jnk] = match_str({tmpdat.sensor_info.label},cst.InvolvedChans);
% subplot(3,4,1); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'hlmarker','*','hlmarkersize',4,'electrodes','on','style','straight','newfig',0); %colorbar
% title('delay maps'); ylabel('grad1 (dB/dx)'); set(gca,'clim',clim);
tmpdat = ts_data_selection(deldat,'chantype','grad2'); x2 = tmpdat.averages.data;
[highlight,jnk] = match_str({tmpdat.sensor_info.label},cst.InvolvedChans);
% subplot(3,4,5); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'hlmarker','*','hlmarkersize',4,'electrodes','on','style','straight','newfig',0);% colorbar
% ylabel('grad2 (dB/dy)');  set(gca,'clim',clim);
% combine grad1 & grad2 delay maps
x  = nan(102,1);
ix = ~isnan(x1) & isnan(x2);  x(ix) = x1(ix);
ix = isnan(x1) & ~isnan(x2);  x(ix) = x2(ix);
ix = ~isnan(x1) & ~isnan(x2); x(ix) = mean([x1(ix) x2(ix)],2);
ex = find(isnan(x));
m  = true(102,1);
m(ex) = 0;
x(ex) = 0;
tmpdat.averages.data = x;
% subplot(3,4,1); ts_ezplot(tmpdat,'mask',m,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',refix,'hlmarker','*','hlmarkersize',4,'hlcolor','r','electrodes','on','style','straight','newfig',0); %colorbar
% title('delay map'); ylabel('grads'); set(gca,'clim',clim);

tmpref  = tmpdat.sensor_info(refix).label;
tmpdat  = ts_data_selection(tmpdat,'badchans',ex);
refix   = strmatch(tmpref,{tmpdat.sensor_info.label});
subplot(3,4,1); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',refix,'hlmarker','*','hlmarkersize',12,'hlcolor','r','electrodes','on','style','straight','newfig',0); %colorbar
title('delay map'); ylabel('grads'); set(gca,'clim',clim); axis equal


[highlight,jnk] = match_str({epoch_data.sensor_info.label},cst.InvolvedChans);
subplot(3,4,5); dewar;view(45,30);%view(2); % view(-90,30) %view(0,90); %view(2);
tmpix = strmatch(tmpref,{epoch_data.sensor_info.label});
hold on; plot3(epoch_data.sensor_info(tmpix).loc(1,4),epoch_data.sensor_info(tmpix).loc(2,4),epoch_data.sensor_info(tmpix).loc(3,4),'*','MarkerSize',12,'Color','r');
axis equal; axis off; title('involved grads');
subplot(3,4,10); plot(cst.Delays,cst.Theta3D,'.'); axis tight; 
title(sprintf('R=%g (p=%g)',cst.R,cst.p)); xlabel('delay (s)'); ylabel('angle (deg)');
xlim=get(gca,'xlim'); set(gca,'xlim',[xlim(1) .2]); set(gca,'ylim',[0 180]);lsline

% traveling wave
cst = clusters_Cmax_all(1).epochs(travtrial);
trav = ts_data_selection(epoch_data,'toilim',cst.HistTime+[-5 5],'chanlabel',cst.InvolvedChans);
trav = ts_preproc(trav,'bpfilter','yes','bpfreq',[.01 4]);

tpad = 1; badchans = 44;
trav = ts_data_selection(trav,'toilim',cst.HistTime+[-tpad tpad],'badchans',badchans);
TRAV = trav;
% for k = 1:trav.num_sensors
%   TRAV.epochs.data(k,:) = TRAV.epochs.data(k,:) / max(abs(TRAV.epochs.data(k,round(length(TRAV.epochs.time)/2)+[-20:20])));
% end
subplot(3,4,3); imagesc(TRAV.epochs.time,1:TRAV.num_sensors,TRAV.epochs.data);
set(gca,'clim',[-.05E-7 .05E-7]); set(gca,'xlim',cst.HistTime+[-tpad tpad]);
title('traveling wave');
subplot(3,4,7); plot(trav.epochs.time,TRAV.epochs.data); axis tight
set(gca,'xlim',cst.HistTime+[-tpad tpad]);
xlabel('time (s)'); set(gca,'ylim',[-.13E-7 .13E-7])

refix  = round(strmatch(cst.RefChan,{epoch_data.sensor_info.label})/2);
sel    = setdiff(1:length(cst.InvolvedChans),badchans);
ch     = match_str({epoch_data.sensor_info.label},cst.InvolvedChans(sel));
deldat = epoch_data;
deldat.epochs.data = nan(deldat.num_sensors,1);
deldat.epochs.data(ch,:) = cst.Delays(sel);
deldat.epochs.time = 0;
deldat.averages = deldat.epochs;
deldat = rmfield(deldat,'epochs');
tmpdat = ts_data_selection(deldat,'chantype','grad1'); x1 = tmpdat.averages.data;
[highlight,jnk] = match_str({tmpdat.sensor_info.label},cst.InvolvedChans);
% subplot(3,4,4); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'hlmarker','*','hlmarkersize',4,'electrodes','on','style','straight','newfig',0); %colorbar
% title('delay maps'); ylabel('grad1 (dB/dx)'); set(gca,'clim',clim);
tmpdat = ts_data_selection(deldat,'chantype','grad2'); x2 = tmpdat.averages.data;
[highlight,jnk] = match_str({tmpdat.sensor_info.label},cst.InvolvedChans);
% subplot(3,4,8); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',highlight,'hlmarker','*','hlmarkersize',4,'electrodes','on','style','straight','newfig',0);% colorbar
% ylabel('grad2 (dB/dy)');  set(gca,'clim',clim);
% combine grad1 & grad2 delay maps
x  = nan(102,1);
ix = ~isnan(x1) & isnan(x2);  x(ix) = x1(ix);
ix = isnan(x1) & ~isnan(x2);  x(ix) = x2(ix);
ix = ~isnan(x1) & ~isnan(x2); x(ix) = mean([x1(ix) x2(ix)],2);
ex = find(isnan(x));
m  = true(102,1);
m(ex) = 0;
x(ex) = 0;
tmpdat.averages.data = x;
% subplot(3,4,4); ts_ezplot(tmpdat,'mask',m,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',refix,'hlmarker','*','hlmarkersize',4,'hlcolor','r','electrodes','on','style','straight','newfig',0); %colorbar
% title('delay map'); ylabel('grads'); set(gca,'clim',clim);

tmpref  = tmpdat.sensor_info(refix).label;
tmpdat  = ts_data_selection(tmpdat,'badchans',ex);
refix   = strmatch(tmpref,{tmpdat.sensor_info.label});
subplot(3,4,4); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',refix,'hlmarker','*','hlmarkersize',12,'hlcolor','r','electrodes','on','style','straight','newfig',0); %colorbar
title('delay map'); ylabel('grads'); set(gca,'clim',clim); axis equal

% x  = zeros(102,1);
% ix = x1~=0 & x2==0; x(ix) = x1(ix);
% ix = x1==0 & x2~=0; x(ix) = x2(ix);
% ix = x1~=0 & x2~=0; x(ix) = mean([x1(ix) x2(ix)],2);
% tmpdat.averages.data = x;
% subplot(3,4,4); ts_ezplot(tmpdat,'topoplot',1,'layout',params.layout,'toprows',1,'topcols',1,'highlight',refix,'hlmarker','*','hlmarkersize',4,'hlcolor','r','electrodes','on','style','straight','newfig',0); %colorbar
% title('delay map'); ylabel('grads'); set(gca,'clim',clim);

[highlight,jnk] = match_str({epoch_data.sensor_info.label},cst.InvolvedChans);
subplot(3,4,8); dewar;view(45,30);%view(2); %view(90,30)% view(0,90); %view(2);
tmpix = strmatch(tmpref,{epoch_data.sensor_info.label});
hold on; plot3(epoch_data.sensor_info(tmpix).loc(1,4),epoch_data.sensor_info(tmpix).loc(2,4),epoch_data.sensor_info(tmpix).loc(3,4),'*','MarkerSize',12,'Color','r');
axis equal; axis off; title('involved grads');
subplot(3,4,11); plot(cst.Delays,cst.Theta3D,'.'); axis tight;
title(sprintf('R=%g (p=%g)',cst.R,cst.p)); xlabel('delay (s)');
xlim=get(gca,'xlim'); set(gca,'xlim',[xlim(1) .2]); set(gca,'ylim',[0 180]); lsline
text(.01,160,sprintf('speed=%gm/s',cst.Speed));


