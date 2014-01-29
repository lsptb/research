function O = IPSPsSectionsTA(x,yF,yIC,Size,Notch,method,FreqRange,Bins,SmoothWindow,PlotsYN,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec)

warning off

Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));
options.WindowStyle = 'normal';

%% ANALYSIS 1
% filter and power spec LFP
dF1 = yF.*1e6;
if ~isempty(Notch), dF1 = Bsfft(dF1,Fs,Notch(1),Notch(2)); end
OF = PowerSpecTA(x,yF,FreqRange,Bins,'Normalized',[]);
switch method
    case {'pwelch','pw'};
        [PxxE,fE] = pwelch(detrend(yF),Bins,[],Bins,Fs);
    case {'multitaper','mtm'};
        [PxxE,fE] = pmtm(detrend(yF),SmoothWindow,Bins,Fs);
end
FR = find(FreqRange(1)<=f & f<=FreqRange(end));
PeakPowerE = max(PxxE(FR)); %#ok<FNDSB>
PxxE = PxxE./PeakPowerE;

% filter IC trace
dIC = Bpfft(yIC*1e3,Fs,0,400); dIC_orig = dIC; O.IC_orig = dIC_orig;
dIC = dIC - offset_voltage;
t = x - x(1); t_orig = x - x(1);
FR2(:,1) = FreqRange;

I = tonic_injected_current; % display this on the figure legend, only if there are no sections
O.tonic_injected_current = I;
O.sections_label_num = sections_label_num;
O.sections_start_sec = sections_start_sec;
O.sections_length_sec = sections_length_sec;


% filter IPSPs
if strcmp('Y',PlotsYN)
    hf = figure; maximize(gcf), subplot(3,1,[1 2]), plot(t,dIC,'-k'), axis tight, set(gcf,'Toolbar','figure')
    assignin('base','t2',t), assignin('base','dIC',dIC), assignin('base','Fs',Fs)
    TA_SlideDisplay_NEW('plot','t2','dIC',[],{'downsample' 'dIC' []},[1 10],{'sgolaySubFun' 'dIC' []},[3 101]);
    uicontrol('Position',[50 75 200 40],'String','Done','Callback','uiresume(gcbf)');
    uiwait(gcf); close(hf);
    dIC = evalin('base','OUT.out1'); O.IC_golay = dIC; xxtt = evalin('base','xt');
    t_factor = round(length(t)/length(xxtt));
    if t_factor>1, dIC_orig = downsample(dIC_orig,t_factor); dF1 = downsample(dF1,t_factor); t = downsample(t,t_factor); Fs = Fs/t_factor; end
else oo1.out1 = sgolayfilt(dIC,2,41);
    assignin('base','OUT',oo1)
    dIC = evalin('base','OUT.out1'); O.IC_golay = dIC;
end

if ~isempty(sections_start_sec)
    for k = 1:length(sections_start_sec)
        dIC{k} = dIC1((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
        dF{k} = dF1((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
    end
else dIC{1} = dIC1;
    dF{1} = dF1;
end


% find IPSPs, discounting those too small or those with spikes
[st,sti] = findpeaks(dIC{k});
[tr,tri] = findpeaks(dIC{k}*-1);
tr = tr.*-1;
if tri(1)<sti(1), tr(1) = []; tri(1) = []; end
if length(tri)>length(sti)
    tri(length(sti)+1:end) = []; tr(length(st)+1:end) = [];
end
if length(tri)<length(sti)
    sti(length(tri)+1:end) = []; st(length(tr)+1:end) = [];
end
stHeight = st - tr;
c = 1; sp = 0;
for k = 1:length(stHeight)
    if stHeight(k)>Size && sp==0;
        if stHeight(k)>20
            ST(end) = []; STi(end) = []; TR(end) = []; TRi(end) = [];
            if c>1; c = c-1; end
            sp = 1;
        else ST(c) = st(k); STi(c) = sti(k); TR(c) = tr(k); TRi(c) = tri(k);
            i_tr(c) = tri(k);
            c = c+1;
        end
    elseif stHeight(k)>Size && sp==1;
        sp = 0;
    end
end

% find IPSP widths, discounting compound IPSPs
SThh = (ST-TR).*0.33; SThh = SThh+TR;
c = 1; O.skipped = 0;
O.IPSPcompound(:,1) = nan(size(round(-0.004*Fs):round(0.15*Fs)));
for k = 1:length(STi)-2
    try
        wetemp = crossing(dIC{k}(TRi(k):end),[],SThh(k));
    catch ME1
        if strcmp(ME1.identifier,'MATLAB:badsubscript')
            O.skipped = O.skipped+1;
            wetemp = [];
        end
    end
    if ~isempty(wetemp)
        wendi(k) = wetemp(1);
        wendi(k) = wendi(k) + TRi(k);
        drev = flipud(dIC{k});
        TRirev(k) = length(drev) - TRi(k);
        wstemp = crossing(drev(TRirev(k):end),[],SThh(k));
        wstarti(k) = wstemp(1);
        wstarti(k) = wstarti(k) + TRirev(k);
        wstarti(k) = length(dIC{k}) - wstarti(k);
        if wendi(k)<STi(k+1) % separate out compound
            ST2(c) = ST(k); STi2(c) = STi(k); TR2(c) = TR(k); TRi2(c) = TRi(k);
            STHH(c) = SThh(k);
            WEndi(c) = wendi(k);
            WStarti(c) = wstarti(k);
            i_TR(c) = i_tr(k);
            c = c+1;
        else
            if STi(k)+round(0.15*Fs) < length(dIC_orig) && STi(k)-round(0.004*Fs) > 1
                O.IPSPcompound(:,end+1) = dIC_orig(STi(k)-round(0.004*Fs):STi(k)+round(0.15*Fs)) - dIC_orig(STi(k));
            end
        end
    else STHH(k) = NaN; %SThh(k) = NaN;
    end
end
wendi(wendi==0) = [];
wstarti(wstarti==0) = [];
STHH(isnan(STHH)) = [];
c = 1; ci = 1;
for k = 1:length(i_TR) % continue to separate out compound
    CurrTri = find(tri==i_TR(k));
    if (CurrTri+1)<=length(sti)
        NextTri = sti(CurrTri+1);
    end
    if NextTri<WEndi(k) && isempty(find(i_TR==NextTri,1))
        % do nothing
    else ST3(c) = ST2(k); STi3(c) = STi2(k); TR3(c) = TR2(k); TRi3(c) = TRi2(k);
        STHH2(c) = STHH(k);
        WEndi2(c) = WEndi(k);
        WStarti2(c) = WStarti(k);
        if TRi3(c)<round(30*Fs)
            ci = ci + 1;
        end
        c = c+1;
    end
end
for k = 2:length(TRi)
    IEIms(k-1) = ((TRi(k) - TRi(k-1))./Fs).*1000;
end
for k = 2:length(TRi3)
    IEImsUni(k-1) = ((TRi3(k) - TRi3(k-1))./Fs).*1000;
end

% height, width, intervals & number
O.n = length(TR3);
O.n30sec = ci;
O.nTotal = length(TR);
O.Amplitude = ST - TR;
O.AmplitudeUni = ST3 - TR3;
O.Widthms = ((wendi - wstarti)/Fs)*1e3;
O.WidthmsUni = ((WEndi2 - WStarti2)/Fs)*1e3;
O.IEIms = IEIms;
O.IEImsUni = IEImsUni;

% average IPSP
dFA = Bpfft(dF,Fs,FR2(1),FR2(2));
dFA = dFA - mean(dFA);
for k = 1:length(STi3)
    if round(0.005*Fs)<STi3(k) && STi3(k)<(length(dIC{k})-round(0.15*Fs))
        O.IPSPsUni(:,k) = dIC_orig(STi3(k)-round(0.005*Fs):STi3(k)+round(0.15*Fs));
        O.IPSPsZeroedUni(:,k) = dIC_orig(STi3(k)-round(0.004*Fs):STi3(k)+round(0.15*Fs)) - dIC_orig(STi3(k));
        O.AngLFPUni(:,k) = dFA(STi3(k)-round(0.004*Fs):STi3(k)+round(0.15*Fs));
        O.AngUni(:,k) = (angle(hilbert(O.AngLFPUni(:,k))));%*180)/pi;
        O.AngLocUni(k) = O.AngUni(round(0.004*Fs)+1,k);
    end
end
for k = 1:length(STi)
    if round(0.01*Fs)<STi(k) && STi(k)<(length(dIC{k})-round(0.15*Fs))
        O.IPSPs(:,k) = dIC_orig(STi(k)-round(0.004*Fs):STi(k)+round(0.15*Fs));
        O.IPSPsZeroed(:,k) = dIC_orig(STi(k)-round(0.004*Fs):STi(k)+round(0.15*Fs)) - dIC_orig(STi(k));
        O.AngLFP(:,k) = dFA(STi(k)-round(0.004*Fs):STi(k)+round(0.15*Fs));
        O.Ang(:,k) = (angle(hilbert(O.AngLFP(:,k))));%*180)/pi;
        O.AngLoc(k) = O.Ang(round(0.004*Fs)+1,k);
    end
end
O.IPSPmeanUni = mean(O.IPSPsZeroedUni,2);
O.IPSPmean = mean(O.IPSPsZeroed,2);

% IPSP slope (on)
wh = 0;
while wh<=50
    try
    % IPSP slope (on) Uni
    O.IPSPmeanUniGOLAY = [];
    [~,msti] = findpeaks(O.IPSPmeanUni,'NPeaks',1);
    [mtr,mtri] = findpeaks(O.IPSPmeanUni*(-1),'NPeaks',1);
    mtr = mtr*(-1);
    bt = mtr*0.2;
    btctemp = crossing(O.IPSPmeanUni(msti:end),[],bt);
    btc = btctemp(1); btc = (btc + msti)/5;
    tp = mtr*0.8;
    tpctemp = crossing(flipud(O.IPSPmeanUni(1:mtri)),[],tp);
    tpc = tpctemp(1); tpc = (mtri - tpc)/5;
    btc = [btc bt]; tpc = [tpc tp];
    coord = tpc - btc;
    O.IPSPSlopeOnUni = coord(2)/coord(1);
    % IPSP slope (off) Uni
    [moff,moffi] = findpeaks(O.IPSPmeanUni(mtri:end),'NPeaks',1);
    mtr2 = mtr - moff;
    tp = (mtr2*0.8) + moff;
    bt = (mtr2*0.2) + moff;
    btctemp = crossing(O.IPSPmeanUni(mtri:end),[],bt);
    btc = btctemp(1); btc = (btc + mtri)/5;
    tpctemp = crossing(O.IPSPmeanUni(mtri:end),[],tp);
    tpc = tpctemp(1); tpc = (tpc + mtri)/5;
    btc = [btc bt]; tpc = [tpc tp];
    coord = tpc - btc;
    O.IPSPSlopeOffUni = coord(2)/coord(1);
    O.IPSPAreaOffUni = sum(O.IPSPmeanUni(mtri:(moffi + mtri)));
    O.IPSPTimemsOffUni = (moffi + mtri - msti)/5;
    catch ME
        if strcmp({'MATLAB:badsubscript'},ME.identifier)
            O.IPSPmeanUni = sgolayfilt(O.IPSPmeanUni,2,5);
            wh = wh+1;
        else wh = 51; %#ok<NASGU>
            rethrow(ME)
        end
    end
    if ~exist('ME','var')
        disp(['golay filtered IPSP: ' num2str(wh) ' iterations'])
        O.GolayIteration = wh;
        if wh>0
            O.IPSPmeanUniGOLAY = O.IPSPmeanUni;
            O.IPSPmeanUni = mean(O.IPSPsZeroedUni,2);
        end
        wh = 51;
    end
    clear ME
end

%% PLOTS 1
% IPSP trace

if strcmp('Y',PlotsYN)
    
    % show all IPSPs
    figure, maximize(gcf), subplot(3,3,[1 2]), hold on
    plot(t,dIC_orig,'-c'), plot(t,dIC{k},'-k')
    scatter(WStarti2./Fs,STHH2,[],'r','Marker','.')
    scatter(WEndi2./Fs,STHH2,[],'r','Marker','.')
    axis tight, xlabel('Time (s)'), ylabel('Membrane Potential (mV)')
    % plot average IPSP
    subplot(3,3,[3 6 9]), hold on
    t_ls = linspace(-4,150,size(O.IPSPsZeroedUni,1));
    plot(t_ls,O.IPSPsZeroedUni,'Color',[0.8 0.8 0.8])
    if isempty(O.IPSPmeanUniGOLAY)
        plot(t_ls,O.IPSPmeanUni,'-k','LineWidth',1), axis tight, xlim([-5 80])
    else
        plot(t_ls,O.IPSPmeanUniGOLAY,':k','LineWidth',1)
        plot(t_ls,O.IPSPmeanUni,'-k','LineWidth',1), axis tight, xlim([-5 80])
    end
    gy = get(gca,'YLim'); if gy(2)>15; ylim([gy(1) 15]), end
    xlabel('Time (ms)'), ylabel({'Potential (mV)'})
    
end

%% ANALYSIS 2
% power spectra of IPSPs
if strcmp('Y',PlotsYN)
    ans1 = inputdlg('Choose a section to cross-correlate','',1,{'10 15'},options);
else ans1{1} = num2str([t(end-10*Fs) t(end)]);
end
ct = str2num(ans1{1}); %#ok<ST2NM>
if length(round(ct(1)*Fs):round(ct(2)*Fs))<Bins/2
    Bins = round(length(round(ct(1)*Fs):round(ct(2)*Fs))/10);
end
OIC = PowerSpecTA(x,yIC*1e3,FreqRange,Bins,'Normalized',[]);
switch method
    case {'pwelch','pw'};
        [PxxE2,fE2] = pwelch(detrend(yIC(round(ct(1)*Fs):round(ct(2)*Fs))*1e3),Bins/2,[],Bins/2,Fs);
    case {'multitaper','mtm'};
        [PxxE2,fE2] = pmtm(detrend(yIC(round(ct(1)*Fs):round(ct(2)*Fs))*1e3),Window,Bins/2,Fs);
    otherwise; error('UNKNOWN METHOD SPECIFIED.  Must be either pwelch or multitaper')
end
FR = round(FreqRange(1)*(Bins/Fs)):round(FreqRange(end)*(Bins/Fs));
PeakPowerE = max(PxxE(FR));
PxxE2 = PxxE2./PeakPowerE;
maxF = max(dF(round(ct(1)*Fs):round(ct(2)*Fs)));
minF = min(dF(round(ct(1)*Fs):round(ct(2)*Fs)));
maxIC = max(dIC{k}(round(ct(1)*Fs):round(ct(2)*Fs)));
minIC = min(dIC{k}(round(ct(1)*Fs):round(ct(2)*Fs)));

%% cross correlations
[co,lags] = xcov(dF(round(ct(1)*Fs):round(ct(2)*Fs)),dIC{k}(round(ct(1)*Fs):round(ct(2)*Fs)));
[coAIC] = xcov(dIC{k}(round(ct(1)*Fs):round(ct(2)*Fs)));
[coALFP] = xcov(dF(round(ct(1)*Fs):round(ct(2)*Fs)));
lagst = lags./Fs;
[~,loc] = findpeaks(co*-1);
st = find(lags==0,1,'first');
stn = st;
loctemp = loc;
locnear = loctemp(findnearest(stn,loctemp));
locnear2 = locnear;
O.XCorrDelay = ((stn-locnear2)/Fs)*1000;

co = (co-min(co))./(max(co)-min(co));
coAIC = (coAIC-min(coAIC))./(max(coAIC)-min(coAIC));
coALFP = (coALFP-min(coALFP))./(max(coALFP)-min(coALFP));

%% PLOTS 2
% plot histogram of IEI
if strcmp('Y',PlotsYN)
    
    subplot(3,3,4), hold on
    bar(0:5:300,histc(IEImsUni,0:5:300),1,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]), axis tight
    xlabel('IEI (ms)'), ylabel('Counts'), set(gca,'Box','on')

    % plot histogram of amplitude
    subplot(3,3,5), hold on
    bar(0:0.2:20,histc(O.AmplitudeUni,0:0.2:20),1,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]), axis tight, xlim([Size Size+15])
    xlabel('Amplitude (mV)'), ylabel('Counts'), set(gca,'Box','on')

    % plot histogram of width
    subplot(3,3,7), hold on
    bar(0:1:60,histc(O.WidthmsUni,0:1:60),1,'FaceColor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]), axis tight
    xlabel('Width (ms)'), ylabel('Counts'), set(gca,'Box','on')

    % plot LFP & IPSP traces
    figure, maximize(gcf), subplot(3,3,[1 2]), hold on, plot(t_orig,dF,':c'), plot(t,dF,'-k'), axis tight
    title('LFP'), ylabel({'Membrane';'Potential (mV)'}), xlim(ct), ylim([minF maxF]), set(gca,'Box','on')
    subplot(3,3,[4 5]), hold on, plot(t,dIC_orig,'-k'), plot(t,dIC{k},':b'), axis tight, xlim(ct), ylim([minIC maxIC])
    title('IPSPs'), xlabel('Time (s)'), ylabel({'Membrane';'Potential (mV)'}), set(gca,'Box','on')

    % plot correlations
    subplot(3,3,9), plot(lagst,co,'-k'), title('LFP:IC X-correlation')
    ylim([0 1]), xlim([-0.5 0.5]), line([0 0],[0 1],'LineStyle',':'), set(gca,'Box','on')
    text(350,100,['Delay: ' num2str(O.XCorrDelay) 'ms'],'FontSize',10,'Units','pixels')
    line([O.XCorrDelay/(-1000) O.XCorrDelay/(-1000)],[0 1],'LineStyle',':','Color','r');
    subplot(3,3,7), plot(lagst,coAIC,'-k'), title('IPSP autocorrelation')
    ylim([0 1]), xlim([-0.5 0.5]), line([0 0],[0 1],'LineStyle',':'), set(gca,'Box','on')
    cu = coAIC;
    [py,px] = findpeaks(cu(round(end/2)+1:end),'NPeaks',1);
    txt1 = sprintf('%s %0.3g','Power: ',py);
    text(290,-30,txt1,'FontSize',10,'Units','pixels')
    text(290,-50,['Freq: ' num2str(round(Fs/px)) 'Hz'],'FontSize',10,'Units','pixels')
    O.ICOscFreq = Fs/px;
    O.ICAutoCorrPower = py;
    subplot(3,3,8), plot(lagst,coALFP,'-k'), title('LFP autocorrelation')
    ylim([0 1]), xlim([-0.5 0.5]), line([0 0],[0 1],'LineStyle',':'), set(gca,'Box','on')
    cu = coALFP;
    [py,px] = findpeaks(cu(round(end/2)+1:end),'NPeaks',1);
    txt1 = sprintf('%s %0.3g','Power: ',py);
    text(290,-30,txt1,'FontSize',10,'Units','pixels')
    text(290,-50,['Freq: ' num2str(round(Fs/px)) 'Hz'],'FontSize',10,'Units','pixels')
    O.LFPOscFreq = Fs/px;
    O.LFPAutoCorrPower = py;

    % plot power spectra
    subplot(3,3,3), hold on, area(OF.f,OF.Pxx,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]), 
    area(fE,PxxE,'FaceColor',[0.4 0.4 0.4])
    axis tight, xlim([2 100]), ylim([0 1]), set(gca,'Box','on')
    title('LFP Spectra'), ylabel({'Power';'({\mu}V^2/Hz.KHz)'})
    subplot(3,3,6), hold on, area(OIC.f,OIC.Pxx,'FaceColor',[0.9 0.9 0.9],'EdgeColor',[0.9 0.9 0.9]), 
    area(fE2,PxxE2,'FaceColor',[0.4 0.4 0.4])
    axis tight, xlim([2 100]), ylim([0 1]), set(gca,'Box','on')
    title('IPSP Spectra'), xlabel('Frequency (Hz)'), ylabel({'Power';'({\mu}V^2/Hz.KHz)'})

end

% phase angle historgram
if strcmp('Y',PlotsYN)
        pn = 2;
    [~,~,~,~,Ang_peak{1},~,Distrib_kurt{1},Peak_Density{1},DataE_st{1},DataE_en{1}] = HistFitTA2(25,0.5,O.AngLocUni);
    P = 0;
    while P==0
        ans3 = inputdlg({'Do another phase plot with tightly filtered LFP? If yes, what freq range? If no, delete contents and press ''Ok'' (not ''Cancel'')'},...
            'phase plot changes',1,{'20 30'});
        if ~isempty(ans3{1})
            FR2(1:2,pn) = str2num(ans3{1}); %#ok<ST2NM>
            dF2 = Bpfft(dF,Fs,FR2(1,pn),FR2(2,pn));
            dF2 = dF2 - mean(dF2);
            for k = 1:length(STi3)
                if round(0.004*Fs)<STi3(k) && STi3(k)<(length(dIC{k})-round(0.15*Fs))
                    O.AngLFPUni(:,k,pn) = dF2(STi3(k)-round(0.004*Fs):STi3(k)+round(0.15*Fs));
                    O.AngUni(:,k,pn) = (angle(hilbert(O.AngLFPUni(:,k,pn)))); %*180)/pi;
                    O.AngLocUni(pn,k) = O.AngUni(round(0.004*Fs)+1,k,pn);
                end
            end
            [~,~,~,~,Ang_peak{pn},~,Distrib_kurt{pn},Peak_Density{pn},DataE_st{pn},DataE_en{pn}] = HistFitTA2(25,0.5,O.AngLocUni(pn,:));
            pn = pn+1;
        else P = 1;
        end
    end
else
    set(0,'DefaultFigureVisible','off')
    [~,~,~,~,Ang_peak{1},~,Distrib_kurt{1},Peak_Density{1},DataE_st{1},DataE_en{1}] = HistFitTA2(25,0.5,O.AngLocUni);
    FR2 = [1,400;10,80;20,30;30,60;2,15]';
    for pn = 2:5
        dF2 = Bpfft(dF,Fs,FR2(1,pn),FR2(2,pn));
        dF2 = dF2 - mean(dF2);
        for k = 1:length(STi3)
            if round(0.004*Fs)<STi3(k) && STi3(k)<(length(dIC{k})-round(0.15*Fs))
                O.AngLFPUni(:,k,pn) = dF2(STi3(k)-round(0.004*Fs):STi3(k)+round(0.15*Fs));
                O.AngUni(:,k,pn) = (angle(hilbert(O.AngLFPUni(:,k,pn)))); %*180)/pi;
                O.AngLocUni(pn,k) = O.AngUni(round(0.004*Fs)+1,k,pn);
            end
        end
        [~,~,~,~,Ang_peak{pn},~,Distrib_kurt{pn},Peak_Density{pn},DataE_st{pn},DataE_en{pn}] = HistFitTA2(25,0.5,O.AngLocUni(pn,:));
    end
    close all
    set(0,'DefaultFigureVisible','on')
end

O.OIC = OIC;
O.OF = OF;
O.PhaseAngleHistFreqRanges = FR2;
O.Fs = Fs;

warning on


end

