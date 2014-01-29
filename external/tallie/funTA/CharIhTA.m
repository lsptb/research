function O = CharIhTA(x,y,offset_voltage,tonic_injected_current,sections_label_num,sections_start_sec,sections_length_sec,baseline_start_sec,baseline_length_sec)
% This function is part of a set of scripts for characterizing cells:
% CharHyperpolStepTA.m
% CharDepolStepTA.m
% CharDepolTonicSpikeTA.m
% CharIhTA.m

Fs = fix(1/(x(2)-x(1))); %Fs = round(1/(data{1}(2,1)-data{1}(1,1)));
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
if ~isnan(offset_voltage), y = (y*1e3) - offset_voltage; else y = (y*1e3); end

if length(sections_length_sec) < length(sections_label_num)
    if isnan(sections_start_sec), sections_start_sec = x(1); end
    sections_length_sec = ones(size(sections_label_num))*sections_length_sec;
end
if length(sections_start_sec) < length(sections_label_num)
    for k = 1:length(sections_length_sec)-1
        sections_start_sec(end+1) = sections_start_sec(k) + sections_length_sec(k);
    end
end
for k = 1:length(sections_start_sec)
    sl{k} = round((sections_start_sec(k)-x(1))*Fs:(sections_start_sec(k)+sections_length_sec(k)-x(1))*Fs);
end

% each {k} below is a section
y_sections = cellfun(@(x) y(x),sl,'Uni',0);
rr = cellfun(@(x) rem(length(x),Fs),y_sections,'Uni',0);
rl = cellfun(@(x) floor(length(x)/Fs),y_sections,'Uni',0);
for k = 1:length(y_sections), y_sections{k}(end-rr{k}+1:end) = []; end
y_sect_mat = cellfun(@(x,y) reshape(x,Fs,y),y_sections,rl,'Uni',0);
y_mean = cellfun(@(x) mean(x,2),y_sect_mat,'Uni',0);

% use this mean to find event locations
dy_mean = cellfun(@(x) moving(diff(moving(x,100)),100),y_mean,'Uni',0);
dy_mean2 = cellfun(@(x) x*(-1),dy_mean,'Uni',0); 
[~,mloc_on] = cellfun(@max,dy_mean2,'Uni',0);
[~,mloc_off] = cellfun(@max,dy_mean,'Uni',0);

end_sect = cellfun(@(x,y) x(y-(Fs/4):y),y_mean,mloc_off,'Uni',0);
[~,lend] = cellfun(@max,end_sect,'Uni',0);
lend = cellfun(@(x,y) x-(Fs/4)+y,mloc_off,lend,'Uni',0);
StepFin1 = cellfun(@(x) round(x-150),lend,'Uni',0);
StepFin2 = cellfun(@(x) round(x-50),lend,'Uni',0);

start_sect = cellfun(@(x,y,z) x(y:round(z/2)),y_mean,mloc_on,StepFin1,'Uni',0);
[~,IhPk] = cellfun(@(x) max(x*(-1)),start_sect,'Uni',0);
IhPk = cellfun(@(x,y) x+y,IhPk,mloc_on,'Uni',0);
IhSta1 = cellfun(@(x) round(x-20),IhPk,'Uni',0);
IhSta2 = cellfun(@(x) round(x+80),IhPk,'Uni',0);

bl = round((baseline_start_sec-x(1))*Fs:(baseline_start_sec+baseline_length_sec-x(1))*Fs);
O.baseline_mV = mean(y(bl));

O.Ih_Peak_mV = cellfun(@(x,y,z) mean(x(y:z)),y_mean,IhSta1,IhSta2,'Uni',0);
O.Ih_FallOffTrace = cellfun(@(x,y,z) x(y:z),y_mean,IhPk,StepFin2,'Uni',0);
O.Ih_End_mV = cellfun(@(x,y,z) mean(x(y:z)),y_mean,StepFin1,StepFin2,'Uni',0);
O.Ih_Diff = cellfun(@(x,y) (x-y)*(-1),O.Ih_Peak_mV,O.Ih_End_mV,'Uni',0);

cmap = [0.7,0,0;0.7,0,0;0,0.7,0;0,0.7,0;0,0,0.7;0,0,0.7];
if I<0, cmap2 = autumn(length(y_sections)+3); elseif I>0, cmap2 = winter(length(y_sections)+3); elseif I==0, cmap2 = bone(length(y_sections)+3); end
cmap2 = flipud(cmap2); cmap2(1:3,:) = []; cmap2 = num2cell(cmap2,2);
figure, hold on, cellfun(@(x,y) plot(1/Fs:1/Fs:length(x)/Fs,x,'Color',y),y_mean,cmap2','Uni',0);
legend([num2str(sections_label_num) ' nA'])
gx = get(gca,'XLim'); gy = get(gca,'YLim');
GX = gx(1)+(range(gx)/20); GY = gy(1)+(range(gy)/20);
if isnan(offset_voltage), text(GX,GY,'(offset voltage unknown)'), end
for k = 1:length(IhSta1)
    scatter([bl(1)/Fs bl(end)/Fs IhSta1{k}/Fs IhSta2{k}/Fs StepFin1{k}/Fs StepFin2{k}/Fs],...
        [O.baseline_mV O.baseline_mV O.Ih_Peak_mV{k} O.Ih_Peak_mV{k} O.Ih_End_mV{k} O.Ih_End_mV{k}],[],cmap,'filled');
end
axis tight, set(gca,'Box','on'), xlabel('Time (s)'), ylabel('Membrane Potential (mV)'), title('Plot to show accuracy in getting Ih data')

O.Fs = Fs;

    
end
