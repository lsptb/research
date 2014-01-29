function O = XCorrTA(x,y1,y2,bpFiltParms,Notch,NormAbs)

warning('off')

% Initial bits n bobs
Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));

% Actions
d1 = y1.*1e6;
d2 = y2.*1e6;
if ~isempty(bpFiltParms);
    d1 = Bpfft(d1,Fs,bpFiltParms(1),bpFiltParms(2));
    d2 = Bpfft(d2,Fs,bpFiltParms(1),bpFiltParms(2));
end
if ~isempty(Notch), d1 = Bsfft(d1,Fs,Notch(1),Notch(2)); end
if ~isempty(Notch), d2 = Bsfft(d2,Fs,Notch(1),Notch(2)); end

[c,lags] = xcov(detrend(d1),detrend(d2));
lagst = lags./Fs;

switch NormAbs
    case 'Normalized'
        c = (c-min(c))./(max(c)-min(c)); %c = cellfun(@(x) x./max(x),c,'Uni',0);
end

cu{1} = c;
[py,px] = findpeaks(cu{1}(round(end/2)+1:end),'NPeaks',1);
O.OscFreq{1} = Fs/px;
O.XCorrPower{1} = py;
for k = 1:2;
    cu{k+1} = envelope(cu{k},'cubic');
    [py,px] = findpeaks(cu{k+1}(round(end/2)+1:end),'NPeaks',1);
    O.OscFreq{k+1} = Fs/px;
    O.XCorrPower{k+1} = py;
end

O.AutoCorr = cu;
O.AutoCorrLagsms = lagst;
O.raw1 = y1;
O.raw2 = y2;
O.Fs = Fs;
O.AutoCorrLags_msPerPoint = (lags(2)-lags(1))/Fs;

warning('on')

end

