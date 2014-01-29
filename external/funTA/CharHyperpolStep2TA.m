function O = CharHyperpolStep2TA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec,step_length)
% This function is part of a set of scripts for characterizing cells:
% CharHyperpolStepTA.m
% CharDepolStepTA.m
% CharDepolTonicSpikeTA.m
% CharIhTA.m
%
% PLEASE MAKE SURE THERE IS NO SPIKE BEFORE THE START OF THE FIRST STEP!!
%

warning('off') %#ok<WNOFF>

% disp('ASSUMING STEP LENGTH OF 0.4 SECS. IF INCORRECT, PLEASE ABORT AND ALTER IN ''CharDepolStepTA.m'' LINE 9.')
% step_length = 0.4;

Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));
O.Fs = Fs;
I = tonic_injected_current;
if isempty(baseline_length_sec), baseline_length_sec = 1/Fs; end
O.tonic_injected_current = I;
O.sections_label_num = sections_label_num;
O.sections_start_sec = sections_start_sec;
O.sections_length_sec = sections_length_sec;
O.baseline_start_sec = baseline_start_sec;
O.baseline_length_sec = baseline_length_sec;
O.y_orig = y;
y = (y*1e3)-offset_voltage;

for k = 1:length(sections_start_sec)
    sl{k} = round((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
end
bl = round((baseline_start_sec-x(1))*Fs:(baseline_start_sec+baseline_length_sec-x(1))*Fs);
O.Baseline_mV = mean(y(bl));
bl_factor = std(diff(y(bl)))*3;
[~,mloc_on] = findpeaks(diff(y),'MinPeakHeight',bl_factor,'NPeaks',1);
mloc_off = round(mloc_on + (step_length*Fs) - 5); % the 'minus five' takes it back to the stimulus artifact peak in the positive direction

% each {k} below is a section
y_sections = cellfun(@(x) y(x),sl,'Uni',0);
rr = cellfun(@(x) rem(length(x),Fs),y_sections,'Uni',0);
rl = cellfun(@(x) floor(length(x)/Fs),y_sections,'Uni',0);
for k = 1:length(y_sections), y_sections{k}(end-rr{k}+1:end) = []; end
y_sect_mat = cellfun(@(x,y) reshape(x,Fs,y),y_sections,rl,'Uni',0);
y_mean = cellfun(@(x) mean(x,2),y_sect_mat,'Uni',0);



% UNFINISHED - WORK FROM HERE.

% are there peaks in here?  If there are, spikes analysis.  If there
% aren't, find the peak in the first half.
for kys = 1:length(y_sect_mat)
    for k = 1:size(y_sect_mat{kys},2)
        [spike_height{kys}{k},spike_rloc{kys}{k}] = findpeaks(y_sect_mat{kys}(mloc_on:mloc_off,k),'MinPeakHeight',-40);
        if isempty(spike_height{kys}{k})
            NoSpike_TF{kys}(k) = 1;
            %Amplitudes{kys}{k} = NaN;
        else NoSpike_TF{kys}(k) = 0;
            %Amplitudes{kys}{k} = (spike_rloc{kys}{k}); % need subtract the baseline - but I haven't written this into the code yet
        end
    end
end
y_mean = cellfun(@(x,y) mean(x(:,logical(y)),2),y_sect_mat,NoSpike_TF,'Uni',0);
[~,IhPk] = cellfun(@(x) max(x(mloc_on:((round(mloc_off-mloc_on)/3))+mloc_on)),y_mean,'Uni',0);
IhPk = cellfun(@(x) x+mloc_on,IhPk,'Uni',0);
IhSta1 = cellfun(@(x) x-20,IhPk,'Uni',0);
IhSta2 = cellfun(@(x) x+80,IhPk,'Uni',0);

O.Peak_mV = cellfun(@(x,y,z) mean(x(y:z)),y_mean,IhSta1,IhSta2,'Uni',0);
VDiff = cellfun(@(x) (x-O.Baseline_mV),O.Peak_mV,'Uni',0);
I_ForDiff = num2cell(sections_label_num);
O.Resistance_Mohms = cellfun(@(x,y) x/y,VDiff,I_ForDiff,'Un',0);

cmap = [0.7,0,0;0.7,0,0;0,0.7,0;0,0.7,0];
if I<0, cmap2 = autumn(length(y_sections)+3); elseif I>0, cmap2 = winter(length(y_sections)+3); elseif I==0, cmap2 = bone(length(y_sections)+3); end
cmap2 = flipud(cmap2); cmap2(1:3,:) = []; cmap2 = num2cell(cmap2,2);
figure, hold on, cellfun(@(x,y) plot(1/Fs:1/Fs:length(x)/Fs,x,'Color',y),y_mean,cmap2','Uni',0);
legend(cellfun(@(x) num2str(x),(num2cell(sections_label_num)),'Uni',0))
for k = 1:length(IhSta1)
    scatter([bl(1)/Fs bl(end)/Fs IhSta1{k}/Fs IhSta2{k}/Fs],...
        [O.Baseline_mV O.Baseline_mV O.Peak_mV{k} O.Peak_mV{k}],[],cmap,'filled');
    for k2 = 1:length(spike_rloc{k}), scatter((spike_rloc{k}{k2}+mloc_on)/Fs,spike_height{k}{k2},[],'*g'), end
end
axis tight, set(gca,'Box','on'), xlabel('Time (s)'), ylabel('Membrane Potential (mV)'), title('Plot to show accuracy in getting event points')

% further analysis for any subsections that contained spikes
for kys = 1:length(spike_height)
    for k = 1:length(spike_height{kys})
        if NoSpike_TF{kys}(k)==0
            O.Spikes_Num{kys}(k) = length(spike_height{kys}{k});
            O.Spikes_ISIs{kys}{k} = (diff(spike_rloc{kys}{k}))/Fs;
            O.Spikes_ISI_median{kys}(k) = median(O.Spikes_ISIs{kys}{k});
            O.Spikes_InstFreq{kys}{k} = 1./((diff(spike_rloc{kys}{k}))/Fs);
            %O.Spikes_Amplitudes{kys}{k} = Amplitudes;
            %O.Spikes_Amplitude_median{kys}{k} = median(O.Spikes_Amplitudes{kys}{k});
        else O.Spikes_Num{kys}(k) = NaN;
            O.Spikes_InstFreq{kys}{k} = NaN;
        end
    end
end

O.Fs = Fs;

warning('on') %#ok<WNON>
    
end

