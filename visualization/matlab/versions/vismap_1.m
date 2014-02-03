%% New visualization tool for continuous data with meaningful layouts
% 19-Apr-2012
function events = vismap(data,layout)
% data = TimeSurfer structure
% layout = layout file (*.lay, ...)
% events = optional structure containing event information

cfg.layout = layout;
lay = prepare_layout(cfg);
[seldat, sellay] = match_str({data.sensor_info.label}, lay.label);
labels = lay(p).label(sellay)';
xpos = lay.pos(sellay, 1)'; 
ypos = lay.pos(sellay, 2)'; 
w = lay.width(sellay)';  
h = lay.height(sellay)'; 
scale = 1/max([xpos(:)+w(:) ypos(:)+h(:)]);
xpos = xpos*scale; % X (normalized)
ypos = ypos*scale; % Y (normalized)
xsiz = w*scale; % width (normalized)
ysiz = h*scale; % height (normalized)

I = 1:1000;
x = data.epochs.time(I);
y = data.epochs.data(seldat,I);

f = figure('color','w','position',[400 200 900 600]);
b = uicontrol('style','pushbutton','string','hi','callback',@update);

for k = 1:size(y,1)
  

h = subplot('position',[0 0 1 1],'units','normalized','parent',f);
p = plot(x,data);

userdata.handles.f = f;
userdata.handles.b = b;
userdata.handles.h = h;
userdata.handles.p = p;
userdata.properties.ylim = ylim;
set(f,'userdata',userdata);

function update(src,evnt)
userdata = get(gcf,'userdata');
h = userdata.handles.h;
p = userdata.handles.p;

y = get(p,'YData');
set(p,'YData',y.*rand(size(y)));

set(h,'ylim',userdata.properties.ylim);
pause(.1);
