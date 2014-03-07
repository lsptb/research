function O = CharHyperpolStepTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec,plot_flag)
% This function is part of a set of scripts for characterizing cells:
% CharHyperpolStepTA.m
% CharDepolStepTA.m
% CharDepolTonicSpikeTA.m

% CHANGELOG
% - 20140223 - adjusted window on moving avg for event detection (100=>10)
% - added plot_flag to suppress plotting; default=1
if nargin<10, plot_flag=1; end
movavgwinsize = 10;

warning('off')

inp{1} = sections_label_num;
inp{2} = sections_start_sec;
inp{3} = sections_length_sec;
inp{4} = baseline_start_sec;
inp{5} = baseline_length_sec;
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
if isnan(baseline_length_sec), baseline_length_sec = 1/Fs; end
if isnan(baseline_start_sec), baseline_start_sec = x(2); end
O.tonic_injected_current = I;
O.offset_voltage = offset_voltage;
O.sections_label_num = sections_label_num;
O.sections_start_sec = sections_start_sec;
O.sections_length_sec = sections_length_sec;
O.baseline_start_sec = baseline_start_sec;
O.baseline_length_sec = baseline_length_sec;
O.y_orig = y;
if ~isnan(offset_voltage), y = (y*1e3) - offset_voltage; else y = (y*1e3); end

if length(sections_start_sec) < length(sections_label_num)
    if isnan(sections_start_sec),
        sections_start_sec = x(2);
    end
end
if length(sections_length_sec) < length(sections_label_num)
    sections_length_sec = ones(size(sections_label_num))*sections_length_sec;
end
if length(sections_start_sec) < length(sections_label_num)
    for k = 1:length(sections_length_sec)-1
        sections_start_sec(end+1) = sections_start_sec(k) + sections_length_sec(k);
    end
end
if sections_start_sec(1) <= x(1), sections_start_sec(1) = x(2); end
for k = 1:length(sections_start_sec)
    sl{k} = round((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
end

% each {k} below is a section
y_sections = cellfun(@(x) y(x),sl,'Uni',0);
rr = cellfun(@(x) rem(length(x),Fs),y_sections,'Uni',0);
rl = cellfun(@(x) floor(length(x)/Fs),y_sections,'Uni',0);
% if 0
  for k = 1:length(y_sections), y_sections{k}(end-rr{k}+1:end) = []; end
% end
y_sect_mat = cellfun(@(x,y) reshape(x,Fs,y),y_sections,rl,'Uni',0);
y_mean = cellfun(@(x) mean(x,2),y_sect_mat,'Uni',0);

% use this mean to find event locations
dy_mean = cellfun(@(x) moving(diff(moving(x,movavgwinsize)),movavgwinsize),y_mean,'Uni',0);
dy_mean2 = cellfun(@(x) x*(-1),dy_mean,'Uni',0); 
[~,mloc_on] = cellfun(@max,dy_mean2,'Uni',0);
[~,mloc_off] = cellfun(@max,dy_mean,'Uni',0);

end_sect = cellfun(@(x,y) x(y-(Fs/4):y),y_mean,mloc_off,'Uni',0);
end_sect2 = cellfun(@(x,y) x(y:end),y_mean,mloc_off,'Uni',0);
[~,IhPk2] = cellfun(@max,end_sect2,'Uni',0);
IhPk2 = cellfun(@(x,y) x+y,IhPk2,mloc_off,'Uni',0);
IhSta12 = cellfun(@(x) round(x-20),IhPk2,'Uni',0);
IhSta22 = cellfun(@(x) round(x+80),IhPk2,'Uni',0);
for k = 1:length(IhSta22)
    if IhSta22{k} > length(y_mean{k})
        IhSta22{k} = length(y_mean{k});
    end
end
[~,lend] = cellfun(@max,end_sect,'Uni',0);
lend = cellfun(@(x,y) x-(Fs/4)+y,mloc_off,lend,'Uni',0);
StepFin1 = cellfun(@(x) round(x-200),lend,'Uni',0);
StepFin2 = cellfun(@(x) round(x-100),lend,'Uni',0);

start_sect = cellfun(@(x,y,z) x(y:round(((z-y)/2)+y)),y_mean,mloc_on,StepFin1,'Uni',0);
[~,IhPk] = cellfun(@(x) max(x*(-1)),start_sect,'Uni',0);
IhPk = cellfun(@(x,y) x+y,IhPk,mloc_on,'Uni',0);
IhSta1 = cellfun(@(x) round(x-20),IhPk,'Uni',0);
for k = 1:length(IhSta1)
    if IhSta1{k} < 1
        IhSta1{k} = 1;
    end
end
IhSta2 = cellfun(@(x) round(x+80),IhPk,'Uni',0);
O.IhStd = cell2num(cellfun(@(x,y,z) std(x(y:z)),y_mean,IhSta1,IhSta2,'Uni',0));

bl = round((baseline_start_sec-x(1))*Fs:(baseline_start_sec+baseline_length_sec-x(1))*Fs);
O.Baseline_mV = mean(y(bl));

O.Ih_Peak_mV = cellfun(@(x,y,z) mean(x(y:z)),y_mean,IhSta1,IhSta2,'Uni',0);
O.Ih_Peak2_mV = cellfun(@(x,y,z) mean(x(y:z)),y_mean,IhSta12,IhSta22,'Uni',0);
Ih_FallOffTrace = cellfun(@(x,y,z) x(y:z),y_mean,IhPk,StepFin2,'Uni',0);
O.Ih_End_mV = cellfun(@(x,y,z) mean(x(y:z)),y_mean,StepFin1,StepFin2,'Uni',0);
O.Ih_Diff1 = cell2num(cellfun(@(x,y) -(x-y),O.Ih_Peak_mV,O.Ih_End_mV,'Uni',0));
O.Ih_Diff2 = cell2num(cellfun(@(x) (x-O.Baseline_mV),O.Ih_Peak2_mV,'Uni',0));
for k = 1:length(O.Ih_Diff2)
    if O.Ih_Diff2(k) > 15
        O.Ih_Diff2(k) = NaN;
    end
end
O.Ih_mV = O.Ih_Diff1 + O.Ih_Diff2;
for k = 1:length(O.Ih_Diff1)
    if O.Ih_Diff1(k) <= O.IhStd(k)*2.5
        O.IsThereIhYNeach(k) = 0;
    else O.IsThereIhYNeach(k) = 1;
    end
end

VDiff = cellfun(@(x) (x-O.Baseline_mV),O.Ih_End_mV,'Uni',0);
I_ForDiff = num2cell(sections_label_num);
O.Resistance_Mohms = cellfun(@(x,y) x/y,VDiff,I_ForDiff,'Uni',0);

O.IsThereIhYNany = any(O.IsThereIhYNeach);
O.IsThereIhYNlast = O.IsThereIhYNeach(end);

if O.IsThereIhYNlast==0
    cmap = [0.7,0,0;0.7,0,0;0.5,0.5,0.5;0.5,0.5,0.5;0,0,0.7;0,0,0.7;0.5,0.5,0.5;0.5,0.5,0.5];
else cmap = [0.7,0,0;0.7,0,0;0,0.7,0;0,0.7,0;0,0,0.7;0,0,0.7;0,0.7,0;0,0.7,0];
end
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
  legend(cellfun(@(x) num2str(x),(num2cell(sections_label_num)),'Uni',0))
  gx = get(gca,'XLim'); gy = get(gca,'YLim');
  GX = gx(1)+(range(gx)/20); GY = gy(1)+(range(gy)/20);
  if isnan(offset_voltage), text(GX,GY,'(offset voltage unknown)'), end
  for k = 1:length(IhSta1)
      scatter([bl(1)/Fs bl(end)/Fs IhSta1{k}/Fs IhSta2{k}/Fs StepFin1{k}/Fs StepFin2{k}/Fs IhSta12{k}/Fs IhSta22{k}/Fs ],...
          [O.Baseline_mV O.Baseline_mV O.Ih_Peak_mV{k} O.Ih_Peak_mV{k} O.Ih_End_mV{k} O.Ih_End_mV{k} O.Ih_Peak2_mV{k} O.Ih_Peak2_mV{k}],[],cmap,'filled');
  end
  axis tight, set(gca,'Box','on'), xlabel('Time (s)'), ylabel('Membrane Potential (mV)'), title('Plot to show accuracy in getting event points')
end
O.Fs = Fs;
O.step_sections_x = [1/Fs:1/Fs:length(y_mean{1})/Fs];
for k = 1:length(y_mean)
    O.step_sections(:,k) = y_mean{k};
    O.Ih_FallOffTrace(1:length(Ih_FallOffTrace{k}),k) = Ih_FallOffTrace{k};
end
O.Ih_FallOffTrace(O.Ih_FallOffTrace==0) = NaN;
try
  FOsg = sgolayfilt(O.Ih_FallOffTrace,1,201);
catch
  FOsg = smooth(O.Ih_FallOffTrace);
end
try
  for k = 1:size(FOsg,2)
      mi = nanmax(FOsg(:,k));
      ma = nanmin(FOsg(:,k));
      mh = mi+((ma-mi)/3);
      mctemp = crossing(FOsg(:,k),1:length(FOsg(:,k)),mh);
      if ~isempty(mctemp)
        O.IhDecayToThird_ms(k) = mctemp(1)/Fs*1000;
      end
  end
  O.DecayToThirdlast = O.IhDecayToThird_ms(end);
catch
  O.DecayToThirdlast = [];
end
O.Ih_End_mV = cell2num(O.Ih_End_mV);
O.Ih_Peak_mV = cell2num(O.Ih_Peak_mV);
O.Ih_Peak2_mV = cell2num(O.Ih_Peak2_mV);
O.Resistance_Mohms = cell2num(O.Resistance_Mohms)*(-1);

warning('on')

end
