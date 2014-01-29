function O = PowerSpecWobbleTA(x,y,FreqRange,Notch,mtmBandWidth,tonic_injected_current,varargin) % (removed NormAbs, removed GolayFiltParms & replaced with [1 501])
%
% (Please note, this is optimized for frequencies above 2Hz. A few lines
% may need changing if you are working with very low frequencies).
% If there are multiple sections, varargin is as follows:
% varargin = {[sections_label_num] [sections_start_sec] [sections_length_sec]};
%
% If using MultiFunTA, and you don't know the current injection:
% Ce = {zeros(length(Cf),1)};
%

%% Initial bits n bobs

warning('off')

for k = 1:length(varargin)
    if iscell(varargin{k})
        if isnumeric(varargin{k}{1})
            varargin{k} = cell2num(varargin{k});
        end
    end
end

Fs = round(1/(x(2)-x(1)));
d1 = y.*1e3;
O.raw = y;
if ~isempty(Notch), d1 = Bsfft(d1,Fs,Notch(1),Notch(2)); end
d1_filt = sgolayfilt(d1,1,501); % sgolayfilt(d1,GolayFiltParms(1),GolayFiltParms(end));
I = tonic_injected_current; % display this on the figure legend, only if there are no sections
O.tonic_injected_current = I;
O.FreqRange = FreqRange;
if ~isempty(varargin)
    sections_label = varargin{1};
    sections_start_sec = varargin{2};
    sections_length_sec = varargin{3};
    
    if length(sections_start_sec) < length(sections_label)
        if isnan(sections_start_sec),
            sections_start_sec = x(2);
        end
    end

    if length(sections_length_sec) < length(sections_label)
        sections_length_sec = ones(size(sections_label))*sections_length_sec;
    end
    if sections_start_sec(1) <= x(1), sections_start_sec(1) = x(2); end
    if length(sections_start_sec) < length(sections_label)
        for k = 1:length(sections_length_sec)-1
            sections_start_sec(end+1) = sections_start_sec(k) + sections_length_sec(k);
        end
    end
    if isnumeric(sections_label)
        O.sections_label = sections_label;
    else O.sections_label = cell2char(sections_label');
    end
    O.sections_start_sec = sections_start_sec;
    O.sections_length_sec = sections_length_sec;
    for k = 1:length(sections_start_sec)
        d{k} = d1((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
        df{k} = d1_filt((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
        xt{k} = x((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
    end
else d{1} = d1;
    df{1} = d1_filt;
    xt{1} = x;
end

ScaleUp = 4;

%%  find the sections of trace to analyse and analyse them

for k = 1:length(d)
    
    % find spikes and separate out the sections between them for fft later:
    [~,pkl] = findpeaks(d{k},'MinPeakHeight',mean(d{k})+15);
    df2 = sgolayfilt(df{k},1,501); % sgolayfilt(df{k},GolayFiltParms(1),GolayFiltParms(end));
    c1 = 1;
    if ~isempty(pkl) % if there are any spikes
        pkl(end+1) = length(d{k});
        for kp = 1:length(pkl)-1
            dtemp = df2(pkl(kp)+(0.02*Fs):pkl(kp+1)-(0.02*Fs)); % looking for the end of the AHP
            if length(dtemp)/Fs >= 0.5 % if it's a long enough section to justify spectral analysis
                [~,pkl2] = findpeaks(dtemp,'NPeaks',1);
                if ~isempty(pkl2) % if the AHP appears to end before the end of that section
                    d_sect(1:length(d{k}(pkl(kp)+(0.02*Fs)+pkl2:pkl(kp+1)-(0.02*Fs))),c1) = d{k}(pkl(kp)+(0.02*Fs)+pkl2:pkl(kp+1)-(0.02*Fs));
                    d_sect_filt(1:length(d{k}(pkl(kp)+(0.02*Fs)+pkl2:pkl(kp+1)-(0.02*Fs))),c1) = df{k}(pkl(kp)+(0.02*Fs)+pkl2:pkl(kp+1)-(0.02*Fs));
                    marker(c1) = length(d{k}(pkl(kp)+(0.02*Fs)+pkl2:pkl(kp+1)-(0.02*Fs)))+1;
                    c1 = c1+1;
                end
            end
        end
        % the section before the first spike:
        dtemp = df2(1:pkl(1)-(0.02*Fs));
        if length(dtemp)/Fs >= 0.5
            d_sect(1:length(df{k}(1:pkl(1)-(0.02*Fs))),c1) = d{k}(1:pkl(1)-(0.02*Fs));
            d_sect_filt(1:length(df{k}(1:pkl(1)-(0.02*Fs))),c1) = df{k}(1:pkl(1)-(0.02*Fs));
            marker(c1) = length(df{k}(1:pkl(1)-(0.02*Fs)))+1;
        end
    % for if there are no spikes
    else d_sect = d{k};
        d_sect_filt = df{k};
    end
    
    if ~exist('marker','var'), marker = length(d{k}); end
    
    if exist('d_sect','var')
        
        for km = 1:length(marker)
            if marker(km) < size(d_sect,1)
                d_sect(marker(km):end,km) = NaN;
                d_sect_filt(marker(km):end,km) = NaN;
            end
        end
        
        if ~isempty(varargin)
            if isnumeric(sections_label)
                NAME = genvarname(['sect_' num2str(sections_label(k))]);
            else NAME = genvarname(['sect_' sections_label(k)]);
            end
            O.data_sections.(NAME) = d_sect;
        else O.data_sections = d_sect; %NAME = genvarname(['sect_' num2str(k)]);
        end
        

        % ANALYSIS:
        for ks = 1:size(d_sect,2)
            [Pxx,f] = pmtm(detrend(d_sect(1:marker(ks)-1,ks)),mtmBandWidth,Fs*ScaleUp,Fs);
            FR = find(FreqRange(1)<=f & f<=FreqRange(end));
            PeakPower = max(Pxx(FR));
            AreaPower = trapz(f(FR),Pxx(FR));
            %maxPP = max(PeakPower);
            %maxAP = max(AreaPower);
            %if strcmp('Normalized',NormAbs)
            %    Pxx = Pxx./maxPP;
            %    AreaPower = AreaPower./maxAP;
            %end
            [mn,n] = max(Pxx(FR));
            FRup = FR(1)+n + 2*ScaleUp;
            if FR(1)+n > 4*ScaleUp % THIS WOULD NEED ALTERING IF LOOKING AT VERY LOW FREQ IN THE FUTURE!
                FRdown = FR(1)+n - 2*ScaleUp;
            else FRdown = FR(1)+n;
            end
            [~,n2] = findpeaks(Pxx(FRdown:FRup),'MinPeakHeight',mn*0.85); % for finding the mean freq when peaks are rough
            if ~isempty(n2), n = nanmean(n2); end
            OscFreq = roundn(f(FRdown(1)+round(n)),-1);

            if ~isempty(varargin)
                O.AreaPower.(NAME)(ks) = AreaPower;
                O.PeakPower.(NAME)(ks) = PeakPower;
                O.OscFreq.(NAME)(ks) = OscFreq;
                O.Pxx.(NAME)(:,ks) = Pxx;
                O.f.(NAME)(:,ks) = f;
            else O.AreaPower(ks) = AreaPower;
                O.PeakPower(ks) = PeakPower;
                O.OscFreq(ks) = OscFreq;
                O.Pxx(:,ks) = Pxx;
                O.f(:,ks) = f;
            end
            
        end

        d_sect_length = 0;
        for ks = 1:length(marker)
            d_sect_length = d_sect_length + marker(ks)-1;
        end
        O.total_length_sections_sec(k) = d_sect_length/Fs;
        
        if ~isempty(varargin)
            O.Pxx_mean.(NAME) = nanmean(O.Pxx.(NAME),2);
            O.f_mean.(NAME) = O.f.(NAME)(:,1);
            O.AreaPower_mean.(NAME) = trapz(f(FR),O.Pxx_mean.(NAME)(FR));
            O.PeakPower_mean.(NAME) = max(O.Pxx_mean.(NAME)(FR));
            [mn,n] = max(O.Pxx_mean.(NAME)(FR));
            FRup = FR(1)+n + 2*ScaleUp;
            if FR(1)+n > 4*ScaleUp
                FRdown = FR(1)+n - 2*ScaleUp;
            else FRdown = FR(1)+n;
            end
            [~,n2] = findpeaks(O.Pxx_mean.(NAME)(FRdown:FRup),'MinPeakHeight',mn*0.85);
            if ~isempty(n2), n = nanmean(n2); end
            O.OscFreq_mean.(NAME) = roundn(O.f_mean.(NAME)(FRdown(1)-1+round(n)),-1);
            O.data.(NAME) = d{k};
            O.t.(NAME) = xt{k};
        else O.Pxx_mean = nanmean(O.Pxx,2);
            O.f_mean = O.f(:,1);
            O.AreaPower_mean = trapz(f(FR),O.Pxx_mean(FR));
            O.PeakPower_mean = max(O.Pxx_mean(FR));
            [mn,n] = max(O.Pxx_mean(FR));
            FRup = FR(1)+n + 2*ScaleUp;
            if FR(1)+n > 4*ScaleUp
                FRdown = FR(1)+n - 2*ScaleUp;
            else FRdown = FR(1)+n;
            end
            [~,n2] = findpeaks(O.Pxx_mean(FRdown:FRup),'MinPeakHeight',mn*0.85);
            if ~isempty(n2), n = nanmean(n2); end
            O.OscFreq_mean = roundn(O.f_mean(FRdown(1)-1+round(n)),-1);
            O.data = d{k};
            O.t = xt{k};
        end
        O.Fs = Fs;
        O.Pxx_HzPerBin = 1/ScaleUp;
        if d_sect_length/Fs<2
            O.LessThan2Secs = 'True';
        else O.LessThan2Secs = 'False';
        end
        
    else O.LessThan2Secs = 'True';
        
    end

    clear marker d_sect d_sect_filt
    
end

warning('on')

end

