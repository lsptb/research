function O = CharDepolTonicSpikesTA(x,yF,yIC,bpFiltParms,Notch,PlotsYN,offset_voltage,tonic_injected_current,varargin)
% This function is part of a set of scripts for characterizing cells:
% CharHyperpolStepTA.m
% CharDepolStepTA.m
% CharDepolTonicSpikeTA.m
% CharIhTA.m
%
% If there are multiple sections, varargin is as follows:
% varargin = {[sections_label_num] [sections_start_sec] [sections_length_sec]};
% If the section is to be taken as a whole (i.e. no sections_start_sec,
% etc), leave the last three inputs as empty matrices, but still give the
% tonic injected current.
%
% If using MultiFunTA, and you don't know the current injection:
% Ce = {zeros(length(Cf),1)};
%

warning off %#ok<WNOFF>

for k = 1:length(varargin)
    if iscell(varargin{k}), varargin{k} = cell2num(varargin{k}); end
end

Fs = round(1/(x(2)-x(1)));
I = tonic_injected_current;
O.tonic_injected_current = I;
O.offset_voltage = offset_voltage;
if ~isempty(varargin)
    sections_label_num = varargin{1};
    sections_start_sec = varargin{2};
    sections_length_sec = varargin{3};
    O.sections_label_num = sections_label_num;
    O.sections_start_sec = sections_start_sec;
    O.sections_length_sec = sections_length_sec;
end
if ~isnan(offset_voltage), yIC = (yIC*1e3) - offset_voltage; else yIC = (yIC*1e3); end
if ~isempty(bpFiltParms), yF = Bpfft(yF,Fs,bpFiltParms(1),bpFiltParms(2)); end
if ~isempty(Notch), yF = Bsfft(yF,Fs,Notch(1),Notch(2)); end
yF = yF - mean(yF);
if exist('sections_start_sec','var')
    
    if length(sections_start_sec) < length(sections_label_num)
        if isnan(sections_start_sec),
            sections_start_sec = x(2);
        end
    end

    if length(sections_length_sec) < length(sections_label_num)
        sections_length_sec = ones(size(sections_label_num))*sections_length_sec;
    end
    if sections_start_sec(1) <= x(1), sections_start_sec(1) = x(2); end
    if length(sections_start_sec) < length(sections_label_num)
        for k = 1:length(sections_length_sec)-1
            sections_start_sec(end+1) = sections_start_sec(k) + sections_length_sec(k);
        end
    end
    for k = 1:length(sections_start_sec)
        y1{k} = yIC((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
        y2{k} = yF((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
    end
else y1{1} = yIC;
    y2{1} = yF;
end

for k = 1:length(y1)
    thresh(k) = -15; %-(sqrt(mean(y1{k}.^2))) + (std(y1{k})*4);
    [~,loc] = findpeaks(y1{k},'MinPeakHeight',thresh(k));
    if isempty(loc)
        fo{k}(1) = NaN; so{k}(1) = NaN; Ang{k}(1) = NaN; AngLoc{k}(1) = NaN;
        ISI{k}(1) = NaN; rates(k) = 0; num_spikes(k) = 0;
    else
        num_spikes(k) = length(loc);
        for t = 1:length(loc);
            if loc(t)-200>0 && loc(t)+599<length(y2{k})
                fo{k}(:,t) = y2{k}(loc(t)-200:loc(t)+599);
                so{k}(:,t) = y1{k}(loc(t)-200:loc(t)+599);
                Ang{k}(:,t) = angle(hilbert(fo{k}(:,t))); %*180)/pi;
                AngLoc{k}(t) = Ang{k}(200,t);
            else AngLoc{k}(t) = NaN;
                fo{k}(1:800,t) = deal(NaN);
                so{k}(1:800,t) = deal(NaN);
                Ang{k}(1:800,t) = deal(NaN);
            end
        end
        ISI{k} = diff(loc)/Fs;
        rates(k) = 1/nanmedian(ISI{k});
    end
    if isempty(loc), LOC{k} = NaN; else LOC{k} = loc; end
end

av_spike = cellfun(@(x) nanmean(x,2),so,'Uni',0);
amp_m = cellfun(@(x) nanmax(x,[],1),so,'Uni',0);

O.Time = x - x(1);
O.Fs = Fs;
O.yIC = yIC;
O.yF = yF;

if ~all(cell2mat(cellfun(@isnan,AngLoc,'Uni',0)))
    
    for k = 1:length(AngLoc)
        if isnan(AngLoc{k})
            [Ang_peak{k},~,Distrib_kurt{k},Peak_Density{k},DataE_st{k},DataE_en{k}] = deal(NaN);
        else
            [~,~,~,~,Ang_peak{k},~,Distrib_kurt{k},Peak_Density{k},DataE_st{k},DataE_en{k}] = HistFitTA2(25,0.5,AngLoc{k},PlotsYN);
        end
    end
    if length(Ang_peak)>1, SpikeFieldCoh_RoseAll_TA(Ang_peak,DataE_st,DataE_en,rates,num_spikes), end
    
    for k = 1:length(AngLoc)
    if ~isnan(AngLoc{k})
        for kso = 1:size(so{k},2)
            [~,SpikeSt_i] = max(diff(diff(so{k}(:,kso))));
            SpikeSt_pk = so{k}(SpikeSt_i,kso);
            Spike_amp_mV{k}(kso) = (amp_m{k}(kso)-SpikeSt_pk);
            Spike_hh = ((amp_m{k}(kso)-SpikeSt_pk)/2)+SpikeSt_pk;
            search_sect = so{k}(SpikeSt_i:end,kso);
            cross_i = crossing(search_sect,1:length(search_sect),Spike_hh);
            Spike_width_ms{k}(kso) = (cross_i(2)-cross_i(1))/Fs;
            [~,SpikeEn_i] = crossing(search_sect,1:length(search_sect),search_sect(1));
            search_sect2 = y1{k}(LOC{k}(kso)-200+SpikeSt_i+SpikeEn_i(2)-2:end);
            if Fs/2<length(search_sect2)
                ahp_base = min(search_sect2(1:Fs/2));
            else ahp_base = min(search_sect2(1:end));
            end
            baseline = nanmean(y1{k}(LOC{k}(kso)-200:end));
            ahp_hh = ((baseline-ahp_base)/2)+ahp_base;
            cross_i_ahp = crossing(downsample(search_sect2,10),1:length(downsample(search_sect2,10)),ahp_hh);
            if length(cross_i_ahp)>1
                ahp_halfWidth{k}(kso) = ((cross_i_ahp(2)-cross_i_ahp(1))*10)/Fs;
            else ahp_halfWidth{k}(kso) = NaN;
            end
        end
    else Spike_amp_mV{k} = NaN;
        Spike_width_ms{k} = NaN;
        ahp_halfWidth{k} = NaN;
    end
    end
    
    O.Time_Sections = x(1:800) - x(1);
    O.SpikeRate_persec_AllSect = rates;
    O.Num_Spikes_AllSect = num_spikes;
    if exist('sections_start_sec','var')
        for k = 1:length(sections_start_sec)
            if isnumeric(sections_label_num)
                NAME = genvarname(['sect_' num2str(sections_label_num(k))]);
            else NAME = genvarname(['sect_' sections_label_num(k)]);
            end
            O.Field.(NAME) = y2{k};
            O.Spike.(NAME) = y1{k};
            O.Field_Sections.(NAME) = fo{k};
            O.Spike_Sections.(NAME) = so{k};
            O.Spike_ave.(NAME) = av_spike{k};
            O.Spike_Amplitude.(NAME) = Spike_amp_mV{k};
            O.Spike_Width.(NAME) = Spike_width_ms{k};
            O.SpikeAHP_HalfWidth_ms.(NAME) = ahp_halfWidth{k};
            O.Ang_Sections.(NAME) = Ang{k};
            O.Ang_spike_in_Field.(NAME) = AngLoc{k};
            O.Ang_Distrib_kurt.(NAME) = Distrib_kurt{k};
            O.Ang_Peak_of_Density.(NAME) = Peak_Density{k};
            O.ISI.(NAME) = ISI{k};
            O.SpikeRate_persec.(NAME) = rates(k);
            O.Num_Spikes.(NAME) = num_spikes(k);
        end
    else O.Field = y2{1};
        O.Spike = y1{1};
        O.Field_Sections = fo{1};
        O.Spike_Sections = so{1};
        O.Spike_ave = av_spike{1};
        O.Spike_Amplitude = Spike_amp_mV{k};
        O.Spike_Width = Spike_width_ms{k};
        O.SpikeAHP_HalfWidth_ms = ahp_halfWidth{k};
        O.Ang_Sections = Ang{1};
        O.Ang_spike_in_Field = AngLoc{1};
        O.Ang_Distrib_kurt = Distrib_kurt{1};
        O.Ang_Peak_of_Density = Peak_Density{1};
        O.ISI = ISI{1};
        O.SpikeRate_persec = rates;
        O.Num_Spikes = num_spikes;
    end    
    
else O.NoSpikeSections = cell2num(cellfun(@isnan,AngLoc,'Uni',0));

end

warning on %#ok<WNON>

end

