%% New visualization tool for continuous data with meaningful layouts
% 19-Apr-2012
function events = vismap(data,layout)
% data = TimeSurfer structure
% layout = layout file (*.lay, ...)
% events = optional structure containing event information

cfg.layout = layout;
lay = prepare_layout(cfg);
[seldat, sellay] = match_str({data.sensor_info.label}, lay.label);
labels = lay.label(sellay);
xpos = lay.pos(sellay, 1); 
ypos = lay.pos(sellay, 2); 
w = lay.width(sellay);  
h = lay.height(sellay); 
scale = 1/max([(xpos(:)+w(:))' (ypos(:)+h(:))']);
xpos = (xpos)*scale; % X (normalized)
ypos = (ypos)*scale; % Y (normalized)
xsiz = w*scale; % width (normalized)
ysiz = h*scale; % height (normalized)

I = 1:1000;
x = data.epochs.time(I);
y = data.epochs.data(seldat,I);

pos = get(0,'screensize');
handles.fig = figure('color','w','position',[0 0 pos(3) pos(4)],'units','normalized');%'WindowScrollWheelFcn',@ZoomFunction,'WindowKeyPressFcn',@arrowpress);
handles.btn = uicontrol('style','pushbutton','string','hi','callback',@update);

for k = 1:size(y,1)
  handles.subplot(k) = subplot('position',[xpos(k) ypos(k) xsiz(k) ysiz(k)],'units','normalized','parent',handles.fig);
  handles.plot(k) = plot(x,y(k,:),'b');
  axis tight; set(gca,'xtick',[],'ytick',[]); % maxmin limits
  handles.properties.ylim{k} = ylim;
%   title(labels{k},'fontsize',8);
  text(min(xlim)+.05*diff(xlim),max(ylim)-.15*diff(ylim),labels{k},'fontsize',10,'color','k','fontweight','bold')
end
set(handles.fig,'userdata',handles);

function update(src,evnt)
handles = get(gcf,'userdata');
k = 1;
h = handles.subplot(k);
p = handles.plot(k);
y = get(p,'YData');
set(p,'YData',y.*rand(size(y)));
set(h,'ylim',handles.properties.ylim{k});
pause(.1);
