% load concatenated raw MEG data for one subject
% load corresponding SO paired peaks

toilim   = [500 1200];
gammafreq= [20 60];
sofreq   = [.2 4];
peaktypes= {'pospeak','negpeak'};
winsize  = 3; 					% duration of window centered on peaks

nchan = data.num_sensors;
Fs    = data.sfreq;
pad   = ceil(winsize*Fs/2);%pad 	 = round(winsize*Fs/2) - 1;
Ngrid = round(winsize*Fs); 			% # pts in fit
T     = data.epochs.time;
tix 	= find(T>=toilim(1)&T<=toilim(2));
tix 	= tix(tix-pad>=1 & tix+pad<=length(T));
t     = T(tix);

% slow oscillation bandpass filter
dat 	 = ts_preproc(data,'bpfilter','yes','bpfreq',sofreq,'bandpass_detrend_flag',0);

% loop over channels
for ch = 1:nchan
% select data & peaks bounded by toilim
x   = dat.epochs.data(ch,tix);
gam = ts_freq_filt(data.epochs.data(ch,tix(1)-pad:tix(end)+pad)',Fs,'bandpass',gammafreq,[0 0])';
gam = gam(pad+1:end-pad); 			% remove gamma filter padding
% gam & x should have the same length

% remove peaks outside toilim (necessary to determine # epochs for preallocation)
tmp1 = T(peaks(ch).pospeak)>=toilim(1)&T(peaks(ch).pospeak)<=toilim(2);
tmp2 = T(peaks(ch).negpeak)>=toilim(1)&T(peaks(ch).negpeak)<=toilim(2);
tmp  = tmp1 | tmp2; % keep the same epochs for pospeak & negpeak
peaks(ch).pospeak = peaks(ch).pospeak(tmp);
peaks(ch).negpeak = peaks(ch).negpeak(tmp);

% preallocate interpolated epoch data
nepochs    = length(peaks(ch).pospeak); 	% note: # pospeak = # negpeak
normdata   = zeros(2,Ngrid,nepochs);		% this will be the data we compare
normgamma  = zeros(2,Ngrid,nepochs);
avgdata    = zeros(2,Ngrid);			% this will be the data we compare
avggamma   = zeros(2,Ngrid);

% loop over peaktypes
for pktype = 1:2
peaktype   = peaktypes{pktype};

pks   = peaks(ch).(peaktype);

% flip x if negpeak
if strcmp(peaktype,'negpeak'), x = -x; end

% loop over detections
for k = 1:length(pks)
% select one peak & define window around detection
pk    = pks(k);
tk    = T(pk);
ix    = nearest(t,tk); 			% peak index wrt toilim subset
wn    = ix + [-pad:pad]; 		% window around peak at t=tk
xx    = x(wn); 				% data in window
tt    = t(wn); 				% time in window
gam   = gam(wn);			% gamma in window
SO(k).tk = tk; 				% store peak time
% figure, subplot(1,2,1),plot(tt,xx); title(['peak at t=' num2str(tk)]); subplot(1,2,2),plot(tt,gam)

% find max/min in window & define level
ymax  = max(xx);
ymin  = min(xx);
level = (ymax-ymin)/2;

% find indices to level-crossings surrounding the peak
[crs,tcrs] = crossing(xx,level); 		% indices to all level-crossings in window
mid   = nearest(tt,tk); 		% index to peak in window
beg   = max(crs(crs<mid)); 		% index to crossing before peak closest to peak
end   = min(crs(crs>mid)); 		% index to crossing after peak closest to peak
tbeg  = tt(beg); 			% time of level-crossing before peak
tend  = tt(end); 			% time of level-crossing after peak
%tbeg  = tcrs(crs==beg); 
%tend  = tcrs(crs==end); 

% select subset between t=tstart & t=tstop
% figure, subplot(1,2,1),plot(tt,xx); title(['level=' num2str(level)]); hline(ymin,'k'); hline(ymax,'k'); hline(level,'r'); subplot(1,2,2),plot(tt,gam);
tt    = tt(beg:end);
xx    = xx(beg:end);
gam   = gam(beg:end) - (max(gam)-min(gam))/2;
SO(k).epoch_t0 = tt(1);
SO(k).epoch_tf = tt(end);
SO(k).epoch_Fs = 1/(tt(2)-tt(1));
SO(k).epoch_x  = xx;

% note: xx(1)=xx(end)=level
% subtract level so that xx(1)=xx(end)=0
xx    = xx - xx(1);
% figure, plot(tt,xx); title(['epoch=[tstart,tstop], shifted by ' num2str(level)]);

% undo flip if negpeak
if strcmp(peaktype,'negpeak'), xx = -xx; end

% define artificial time vector t = 0 to 1 sec for this epoch
%Nk    = length(xx); 		% # time points b/w tstart & tstop
%tt    = linspace(0,1,Nk); 	% artificial time vector
% NOTE: tt ~[0,1] may not be necessary using TriScatteredInterp(); maybe shift tt to start at tt=0

% Interpolation: method 1
%figure, subplot(1,2,1),plot(tt,xx); title('artificial time'); subplot(1,2,2),plot(tt,gam)
%ttk 	 = linspace(0,1,Ngrid);
%[Tk,Gk] = griddata(tt,gam,ttk);
%[Tk,Xk] = griddata(tt,xx, ttk);
% (Tk,Xk) defines this epoch which will be used to compare activity around different peaks

% Interpolation: method 2
%Tk 	 = linspace(0,1,Ngrid);
%Gk 	 = interp1(tt,gam,Tk,'linear');
%Xk 	 = interp1(tt,xx,Tk,'linear');

% Interpolation: method 3
G 	= TriScatteredInterp(tt,gam,'linear');
F	= TriScatteredInterp(tt,xx,'linear');
Tk	= linspace(0,1,Ngrid);
Gk 	= G(Tk);
Xk 	= F(Tk);

SO(k).fit_tf = Tk(end);
SO(k).fit_Fs = 1/(Tk(2)-Tk(1));
SO(k).fit_x  = Xk;
SO(k).fit_t0 = Tk(1);
% figure, subplot(1,2,1),plot(Tk,Xk); title('interpolated epoch'); subplot(1,2,2),plot(Tk,Gk);

normgamma(pktype,:,k) = Gk;
normdata(pktype,:,k)  = Xk;

clear Xk Gk Tk tt xx Nk mid beg end tbeg tend crs level ymax ymin
end % end loop over detections

avggamma(pktype,:) = mean(normgamma,3);
avgdata(pktype,:)  = mean(normdata, 3);

end % end loop over peak types

alldata(ch).label 	= data.sensor_info(ch).label;
alldata(ch).cycles 	= SO;
alldata(ch).normdata	= normdata;
alldata(ch).normgamma   = normgamma;
alldata(ch).avgdata 	= avgdata;
alldata(ch).avggamma 	= avggamma;
alldata(ch).num_grid 	= Ngrid;
alldata(ch).avgtime 	= SO(1).fit_t0:1/SO(1).fit_Fs:SO(1).fit_tf;
clear normdata normgamma SO avggamma avgdata

tmpdat1 = alldata(ch).avgdata(1,:);
tmpdat2 = alldata(ch).avgdata(2,:);
tmpgam1 = alldata(ch).avggamma(1,:);
tmpgam2 = alldata(ch).avggamma(2,:);
tmptime = alldata(ch).avgtime;
figure('Name',alldata(ch).label)
subplot(2,1,1),plot(tmptime,tmpdat1,'b',tmptime,tmpdat2,'r'); title('mean SO'); legend(peaktypes{1},peaktypes{2});
subplot(2,1,2),plot(tmptime,tmpgam1,'b',tmptime,tmpgam2,'r'); title('mean gamma');
clear tmpdat1 tmpdat2 tmpgam1 tmpgam2 tmptime
pause

end % end loop over channels
