function O = SpikeFieldCohTA(x,yF,yIC,bpFiltParms,Notch,varargin)

warning off %#ok<WNOFF>

Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));

if nargin==7
    section_start_sec = varargin{1};
    section_length_sec = varargin{2};
end

yIC = detrend(yIC)*1e3;
if ~isempty(bpFiltParms), yF = Bpfft(yF,Fs,bpFiltParms(1),bpFiltParms(2)); end
if ~isempty(Notch), yF = Bsfft(yF,Fs,Notch(1),Notch(2)); end
yF = yF - mean(yF);

if nargin==7
    for k = 1:length(section_start_sec)
        y1{k} = yIC(section_start_sec(k)*Fs:(section_start_sec(k)+section_length_sec(k))*Fs);
        y2{k} = yF(section_start_sec(k)*Fs:(section_start_sec(k)+section_length_sec(k))*Fs);
    end
else y1{1} = yIC; y2{1} = yF;
end

for k = 1:length(y1)
    thresh(k) = -(sqrt(mean(y1{k}.^2))) + (std(y1{k})*4); thresh(2) = 9.5;
    [~,loc] = findpeaks(y1{k},'MinPeakHeight',thresh(k));
    if isempty(loc)
        fo{k} = NaN; so{k} = NaN; Ang{k} = NaN; AngLoc{k} = NaN;
    else
        num_spikes(k) = length(loc);
        for t = 1:length(loc);
            if loc(t)-250>0 && loc(t)+249<length(y2{k})
                fo{k}(:,t) = y2{k}(loc(t)-250:loc(t)+249);
                so{k}(:,t) = y1{k}(loc(t)-250:loc(t)+249);
                Ang{k}(:,t) = angle(hilbert(fo{k}(:,t))); %*180)/pi;
                AngLoc{k}(t) = Ang{k}(end/2,t);
            else AngLoc{k}(t) = NaN;
                idxnan(k) = t;
            end
        end
        if exist('idxnan','var'), fo{k}(:,idxnan) = NaN; so{k}(:,idxnan) = NaN; Ang{k}(:,idxnan) = NaN; end
        ISI{k} = diff(loc)/Fs;
        rates(k) = 1/nanmedian(ISI{k});
    end
end

[~,~,~,~,Ang_peak,~,Distrib_kurt,Peak_Density,DataE_st,DataE_en] = cellfun(@(x) HistFitTA2(25,0.5,x),AngLoc,'Uni',0);
SpikeFieldCoh_RoseAll_TA(Ang_peak,DataE_st,DataE_en,rates,num_spikes)

O.field = y2;
O.spike = y1;
O.field_sections = fo;
O.spike_sections = so;
O.Ang_sections = Ang;
O.AngLoc = AngLoc;
O.Distrib_kurt = Distrib_kurt;
O.Peak_Density = Peak_Density;
O.ISI = ISI;
O.SpikeRate_persec = rates;
O.num_spikes = num_spikes;
O.Fs = Fs;

warning on %#ok<WNON>

end

