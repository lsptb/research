function O = CharDepolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec,step_length,plot_flag)
% This function is part of a set of scripts for characterizing cells:
% CharHyperpolStepTA.m
% CharDepolStepTA.m
% CharDepolTonicSpikeTA.m
% CharIhTA.m
%
% PLEASE MAKE SURE THERE IS NO SPIKE BEFORE THE START OF THE FIRST STEP!!
%
% CHANGELOG:
% 20140223 - added plot_flag to suppress plotting; default=1
% 20140331 - added changes to mloc_on/off from Tallie's version emailed on
%            03-Mar-2014; and assume fac=0 (missing SectionInputsTA)

if nargin<10, plot_flag=1; end

warning('off') %#ok<WNOFF>

inp{1} = sections_label_num;
inp{2} = sections_start_sec;
inp{3} = sections_length_sec;
inp{4} = baseline_start_sec;
inp{5} = baseline_length_sec;
inp{6} = step_length;
for k = 1:length(inp)
    if iscell(inp{k})
        if isnumeric(inp{k}{1})
            inp{k} = cell2num(inp{k});
        end
    end
end

Fs = round(1/(x(2)-x(1)));
O.Fs = Fs;
I = tonic_injected_current;
if isempty(baseline_length_sec), baseline_length_sec = 1/Fs; end
O.tonic_injected_current = I;
O.offset_voltage = offset_voltage;
O.sections_label_num = sections_label_num;
O.sections_start_sec = sections_start_sec;
O.sections_length_sec = sections_length_sec;
O.baseline_start_sec = baseline_start_sec;
O.baseline_length_sec = baseline_length_sec;
O.y_orig = y;
if isnan(baseline_length_sec), baseline_length_sec = 1/Fs; end
if isnan(baseline_start_sec), baseline_start_sec = x(2); end
if ~isnan(offset_voltage)
    y = (y*1e3) - offset_voltage;
else y = (y*1e3);
end

if length(sections_start_sec) < length(sections_label_num)
    if isnan(sections_start_sec), sections_start_sec = x(2); end
end

if length(sections_length_sec) < length(sections_label_num)
    sections_length_sec = ones(size(sections_label_num))*sections_length_sec;
end
    if length(sections_start_sec) < length(sections_label_num)
        for k = 1:length(sections_length_sec)-1
            sections_start_sec(end+1) = sections_start_sec(k) + sections_length_sec(k);
        end
    end
if sections_start_sec(1) <= x(1) || isnan(sections_start_sec(1))
    sections_start_sec(1) = x(2);
end
for k = 1:length(sections_start_sec)
    sl{k} = round((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
end
bl = round((baseline_start_sec-x(1))*Fs:(baseline_start_sec+baseline_length_sec-x(1))*Fs);
O.Baseline_mV = mean(y(bl));
if length(bl)<50
    bl_factor = std(diff(y(bl:bl+50)))*7;
else bl_factor = std(diff(y(bl)))*3;
end
% [~,mloc_on] = findpeaks(diff(y)*(-1),'MinPeakHeight',bl_factor,'NPeaks',1);
% mloc_off = round(mloc_on + (step_length*Fs) - 5); % the 'minus five' takes it back to the stimulus artifact peak in the positive direction
for k = 1:length(sections_start_sec)
    [~,mloc_on{k}] = findpeaks(diff(y(sl{k}))*(-1),'MinPeakHeight',bl_factor,'NPeaks',1);
    mloc_off{k} = round(mloc_on{k} + (step_length*Fs) - round(5*Fs/5000)); % the 'minus five' takes it back to the stimulus artifact peak in the positive direction
end

% each {k} below is a section
y_sections = cellfun(@(x) y(x),sl,'Uni',0);
rr = cellfun(@(x) rem(length(x),Fs),y_sections,'Uni',0);
rl = cellfun(@(x) floor(length(x)/Fs),y_sections,'Uni',0);
for k = 1:length(y_sections), y_sections{k}(end-rr{k}+1:end) = []; end
y_sect_mat = cellfun(@(x,y) reshape(x,Fs,y),y_sections,rl,'Uni',0);

% are there peaks in here?  If there are, spikes analysis.
for kys = 1:length(y_sect_mat)
    for k = 1:size(y_sect_mat{kys},2)
      [spike_height{kys}{k},spike_rloc{kys}{k}] = findpeaks(y_sect_mat{kys}(mloc_on{kys}:mloc_off{kys},k),'MinPeakHeight',O.Baseline_mV+40);
      %[spike_height{kys}{k},spike_rloc{kys}{k}] = findpeaks(y_sect_mat{kys}(mloc_on:mloc_off,k),'MinPeakHeight',O.Baseline_mV+40);
        if isempty(spike_height{kys}{k})
            NoSpike_TF{kys}(k) = 1;
        else NoSpike_TF{kys}(k) = 0;
        end
    end
end

%O.VoltOff = cellfun(@(x,y) mean(mean(x(mloc_off-200:mloc_off-100,logical(y)),2)),y_sect_mat,NoSpike_TF,'Uni',0);
O.VoltOff = cellfun(@(x,y,z) mean(mean(x(z-200:z-100,logical(y)),2)),y_sect_mat,NoSpike_TF,mloc_off,'Uni',0);

y_mean = cellfun(@(x,y) mean(x(:,logical(y)),2),y_sect_mat,NoSpike_TF,'Uni',0);
y_mean2 = cellfun(@(x,y) x(:,logical(~y)),y_sect_mat,NoSpike_TF,'Uni',0);
y_mean2(cellfun(@any,NoSpike_TF)) = [];

VDiff = cellfun(@(x) (x-O.Baseline_mV),O.VoltOff,'Uni',0);
I_ForDiff = num2cell(sections_label_num);
O.Resistance_Mohms = cellfun(@(x,y) x/y,VDiff,I_ForDiff,'Un',0);

cmap = [0.7,0,0;0.7,0,0;0,0.7,0;0,0.7,0];
if I<0
    cmap2 = autumn(length(y_sections)+3);
elseif I>0
    cmap2 = winter(length(y_sections)+3);
elseif I==0 || isnan(I)
    cmap2 = bone(length(y_sections)+3);
end
if plot_flag
  cmap2 = flipud(cmap2); cmap2(1:3,:) = []; cmap2 = num2cell(cmap2,2);
  figure, hold on, cellfun(@(x,y) plot(1/Fs:1/Fs:length(x)/Fs,x,'Color',y),y_mean,cmap2','Uni',0);
  cellfun(@(x,y) plot(1/Fs:1/Fs:length(x)/Fs,x(:,1),'-k'),y_mean2,'Uni',0);
  legend(cellfun(@(x) num2str(x),(num2cell(sections_label_num)),'Uni',0))
  gx = get(gca,'XLim'); gy = get(gca,'YLim');
  GX = gx(1)+(range(gx)/20); GY = gy(1)+(range(gy)/20);
  if isnan(offset_voltage), text(GX,GY,'(offset voltage unknown)'), end
%   for k = 1:length(O.VoltOff)
%       scatter([bl(1)/Fs bl(end)/Fs (mloc_off-200)/Fs (mloc_off-100)/Fs],...
%           [O.Baseline_mV O.Baseline_mV O.VoltOff{k} O.VoltOff{k}],[],cmap,'filled');
%       for k2 = 1:length(spike_rloc{k}), scatter((spike_rloc{k}{k2}+mloc_on)/Fs,spike_height{k}{k2},[],'*g'), end
%   end
  for k = 1:length(O.VoltOff)
      if 0%fac==1;
      scatter([(bl(1)/Fs)-sections_start_sec(1)+x(1) (bl(end)/Fs)-sections_start_sec(1)+x(1) (mloc_off{k}-200)/Fs (mloc_off{k}-100)/Fs],...
          [O.Baseline_mV O.Baseline_mV O.VoltOff{k} O.VoltOff{k}],[],cmap,'filled');
      else scatter([0 0 (mloc_off{k}-200)/Fs (mloc_off{k}-100)/Fs],...
          [O.Baseline_mV O.Baseline_mV O.VoltOff{k} O.VoltOff{k}],[],cmap,'filled');
      end
      for k2 = 1:length(spike_rloc{k}), scatter((spike_rloc{k}{k2}+mloc_on{k})/Fs,spike_height{k}{k2},[],'*g'), end
  end  
  axis tight, set(gca,'Box','on'), xlabel('Time (s)'), ylabel('Membrane Potential (mV)'), title('Plot to show accuracy in getting event points')
  xlim([0 numel(y_mean{1})/Fs]);
end

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
for k = 1:length(y_mean)
    O.y_NoSpike_sect(1:length(y_mean{k}),:) = y_mean{k};
end
O.y_NoSpike_sect(O.y_NoSpike_sect==0) = NaN;
for k = 1:length(y_mean2)
    O.y_Spiking_sect(1:size(y_mean2{k},1),1:size(y_mean2{k},2),k) = y_mean2{k};
end
O.y_Spiking_sect(O.y_NoSpike_sect==0) = NaN;

%O.y_sect_mat = y_sect_mat;

warning('on') %#ok<WNON>
    
end

