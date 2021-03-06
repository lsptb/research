%% New visualization tool for continuous data with meaningful layouts
% 19-Apr-2012
function events = vismap(alldata,layout)
% alldata = TimeSurfer structure
% layout = layout file (*.lay, ...)
% events = optional structure containing event information
clear global alldata
global alldata

cfg.layout = layout;
lay = prepare_layout(cfg);
[seldat, sellay] = match_str({alldata.sensor_info.label}, lay.label);
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
handles.layout.xpos = xpos;
handles.layout.ypos = ypos;
handles.layout.xsiz = xsiz;
handles.layout.ysiz = ysiz;
handles.layout.file = layout;
handles.properties.curr = 1; % index to first sample in current view
handles.properties.L = min(round(5*alldata.sfreq),length(alldata.epochs.time)); % number of samples to display
I = handles.properties.curr:(handles.properties.L-1); % samples to view
x = alldata.epochs.time(I);
y = alldata.epochs.data(seldat,I);
xmin = min(x(:)); xmax = max(x(:));
ymin = min(y(:)); ymax = max(y(:));
handles.properties.sfreq = alldata.sfreq;
% handles.properties.time = alldata.epochs.time;
% handles.properties.seltime = x;
% handles.properties.selchan = 1:alldata.num_sensors; % indices to channels to view
% handles.properties.numchan = alldata.num_sensors;

pos = get(0,'screensize');
handles.fig = figure('color','w','position',[0 0 pos(3) pos(4)],'units','normalized','nextplot','add','WindowScrollWheelFcn',@ZoomFunction);%,'WindowKeyPressFcn',@arrowpress);

BorderWidth = .2;
BorderType  = 'line';%none,etchedin,etchedout,beveledin,beveledout,line
handles.pview  = uipanel('Title','View','FontSize',10,'parent',handles.fig,'BackgroundColor','white','Position',[0 .07 1 .9]   ,'tag','plots','BorderWidth',BorderWidth,'BorderType',BorderType); % plots
handles.pcontrol  = uipanel('Title','Control','FontSize',10,'parent',handles.fig,'BackgroundColor','white','Position',[0 0 1 .07],'BorderWidth',BorderWidth,'BorderType',BorderType  ); % time
handles.leftbtn = uicontrol('style','pushbutton','string','<=','fontsize',12,'parent',handles.pcontrol,'units','normalized','position',[.05 .05 .05 .45],'callback',{@update,'moveleft','wave'});
handles.rightbtn = uicontrol('style','pushbutton','string','=>','fontsize',12,'parent',handles.pcontrol,'units','normalized','position',[.11 .05 .05 .45],'callback',{@update,'moveright','wave'});
handles.leftedit = uicontrol('style','edit','string',num2str(x(1)),'fontsize',12,'parent',handles.pcontrol,'units','normalized','position',[.05 .55 .05 .45],'callback',{@update,'setleft','wave'});
handles.rightedit = uicontrol('style','edit','string',num2str(x(end)),'fontsize',12,'parent',handles.pcontrol,'units','normalized','position',[.11 .55 .05 .45],'callback',{@update,'setright','wave'});
handles.topochk = uicontrol('style','checkbox','string','topo','value',0,'parent',handles.pcontrol,'units','normalized','position',[.26 .55 .05 .45],'callback',{@update,'topo','wave'});
handles.ncoledit = uicontrol('style','edit','string','1','fontsize',12,'parent',handles.pcontrol,'units','normalized','position',[.32 .05 .05 .45]);%,'callback',{@update,'topo','wave'});
handles.nrowedit = uicontrol('style','edit','string','1','fontsize',12,'parent',handles.pcontrol,'units','normalized','position',[.32 .55 .05 .45]);%,'callback',{@update,'topo','wave'});

handles.properties.nplot = size(y,1);
for k = 1:size(y,1)
  handles.subplot(k) = subplot('position',[xpos(k) ypos(k) xsiz(k) ysiz(k)],'units','normalized','parent',handles.pview,'nextplot','add');
  handles.plot(k) = plot(x,y(k,:),'b');
  axis([xmin xmax ymin ymax]); axis off % set(gca,'xtick',[],'ytick',[]); % maxmin limits
  handles.properties.ylim{k} = ylim;
%   handles.text(k) = text(min(xlim)+.05*diff(xlim),max(ylim)-.1*diff(ylim),labels{k},'fontsize',9,'color','k','fontweight','normal');
  handles.text(k) = text(.05,.9,labels{k},'fontsize',9,'color','k','fontweight','normal','units','normalized');
end
handles.properties.clim = ylim;
handles.labels = labels;
set(handles.plot,'ButtonDownFcn',@onclick);
set(handles.fig,'userdata',handles);

function update(src,evnt,action,type)
global alldata
handles = get(gcf,'userdata');
N = length(alldata.epochs.time);
% L = handles.properties.L;
c = handles.properties.curr;
I1 = nearest(alldata.epochs.time,str2num(get(handles.leftedit,'string')));
I2 = nearest(alldata.epochs.time,str2num(get(handles.rightedit,'string')));
L = I2 - I1 + 1;
% determine times to show
switch action
  case 'moveleft'
    I1 = max(1,I1-L+1);
    I2 = min(N,I1+L-1);
  case 'moveright'
    I2 = min(N,I2+L-1);
    I1 = max(1,I2-L+1);
  case {'setleft','setright'}
    I1 = max(1,I1);
    I2 = min(N,I2);    
    L = I2-I1+1;
  otherwise
end
% I1 = max(1,I1);
% I2 = min(N,I2);
I = I1:I2;
c = I1;
L = I2 - I1 + 1;
t1 = alldata.epochs.time(I1);
t2 = alldata.epochs.time(I2);
handles.properties.L = L;
handles.properties.curr = c;
set(handles.leftedit,'string',num2str(t1));
set(handles.rightedit,'string',num2str(t2));
% delete(handles.text);
if strcmp(type,'wave'); %1 % plot waveform
  try delete(handles.topoaxes); end
  if get(handles.topochk,'value')==1
    try
      set(handles.subplot,'visible','off');
      set(handles.plot,'visible','off');
      set(handles.text,'visible','off');
    end
    handles.topoaxes = axes('Position',[0 0 1 1],'units','normalized','visible','on','parent',handles.pview,'nextplot','add');
    axes(handles.topoaxes);
    tmpdat = alldata;
    tmpdat.epochs.data = alldata.epochs.data(:,I);
    tmpdat.epochs.time = alldata.epochs.time(I);
    nrows = str2num(get(handles.nrowedit,'string'));
    ncols = str2num(get(handles.ncoledit,'string'));
    iEEG_flag = double(strcmp(tmpdat.sensor_info(1).typestring,'eeg'));
    showlabels = double((nrows<=1 && ncols<=1));
%     yoffset=0; xoffset=0;
    if nrows*ncols == 1
      yoffset = 0;
      xoffset = 0;
    else
      yoffset = 0;%2                                    % correction for axes extending outside of frame vertically      
      xoffset = 0;%1                                    % correction for showing y-axis labels on in first column
    end    
    xstepsize = 1 / (ncols+xoffset);
    ystepsize = 1 / (nrows+yoffset);
    xstep = mod((1:ncols+xoffset)-1,ncols);
    ystep = mod((1:nrows+yoffset)-1,nrows)+1;
    xpos  = xstep*xstepsize;
    ypos  = 1-ystep*ystepsize;
    xpos  = xpos(1:end-xoffset) + (xoffset/2)*xstepsize;  % correction; does nothing if offset=0
    ypos  = ypos(1:end-yoffset) - (yoffset/2)*ystepsize;  % correction; does nothing if offset=0      
    % set(gcf,'Name',sprintf('%s',cfg.chantype));
    cnt = 0;
    for r = 1:nrows           % row index
      for c = 1:ncols         % column index
        cnt = cnt + 1;
        nn = floor(length(I)/(nrows*ncols));
        samp = (1:nn) + (cnt-1)*nn;
        tmptmp = tmpdat;
        tmptmp.epochs.time = tmptmp.epochs.time(samp);
        tmptmp.epochs.data = tmptmp.epochs.data(:,samp);
        handles.topoaxes(cnt) = subplot('position',[xpos(r) ypos(c) .9*xstepsize .9*ystepsize],'parent',handles.pview,'units','normalized','nextplot','add');
        handles.topoplot(cnt) = topo(tmptmp,'avgovertime','no','zlim',handles.properties.clim,'nrows',1,'ncols',1,...
            'layout',handles.layout.file,'iEEG_flag',iEEG_flag,'showlabels',showlabels,'fig',0,'parent',0,'electrodes',1,'interp','v4');     

    %     hh = topo(tmpdat,'avgovertime','yes','zlim',handles.properties.clim,'nrows',nrows,'ncols',ncols,...
    %         'layout',handles.layout.file,'iEEG_flag',iEEG_flag,'showlabels',showlabels,'fig',handles.fig,'parent',handles.pview,'electrodes',1); 
    %     handles.topoaxes = hh.subplot;
        pause(.1);
      end
    end
    colorbar
    set(handles.topoplot,'ButtonDownFcn',@onclick);
  else
%     if ~ishandle(handles.plot(1))
    try
      set(handles.plot,'visible','on'); 
      set(handles.text,'visible','on');
      set(handles.subplot,'visible','off'); 
      for k = 1:handles.properties.nplot
        h = handles.subplot(k);
        p = handles.plot(k);
        y = alldata.epochs.data(k,I); %get(p,'YData');
        if any(strcmp('setleft',{'setleft','setright'}))
          set(p,'XData',alldata.epochs.time(I),'YData',y);
        else
          set(p,'YData',y);
        end
        set(h,'ylim',handles.properties.ylim{k},'xlim',[t1 t2]);
        axis off
      end
    catch
      xpos = handles.layout.xpos;
      ypos = handles.layout.ypos;      
      x = alldata.epochs.time(I);
      for k = 1:handles.properties.nplot
        y = alldata.epochs.data(k,I); %get(p,'YData');
        handles.subplot(k) = subplot('position',[xpos(k) ypos(k) handles.layout.xsiz(k) handles.layout.ysiz(k)],'units','normalized','parent',handles.pview,'nextplot','add');
        handles.plot(k) = plot(x,y,'b');
        axis([t1 t2 handles.properties.ylim{k}]); axis off
        handles.text(k) = text(.05,.9,handles.labels{k},'fontsize',9,'color','k','fontweight','normal','units','normalized');
      end
      set(handles.plot,'ButtonDownFcn',@onclick);
    end
%     end
    pause(.1);
  end
end
if strcmp(type,'spectrogram') % stft % plot spectrogram
  tscale = 100; % ms
  Nfft = round(.75*handles.properties.L);
  Noverlap = Nfft - 10;
  fs = handles.properties.sfreq;
  t = (-Nfft/2+.5):(Nfft/2-.5); % time in samples for defining window
  f = 0:.5:30;%fs/2;
  sigma = (tscale/1000)*fs;
  W = exp(-(t/sigma).^2); % window function for gabor transform  
  for k = 1:handles.properties.nplot
    h = handles.subplot(k);
    p = handles.plot(k);
    y = get(p,'YData');
%     [Y,F,T,P] = spectrogram(double(y),double(W),[],f,fs);
    [Y,F,T,P] = spectrogram(double(y),double(W),Noverlap,f,fs);
    Z = 10*log10(abs(P));
    axes(h);
    surf(T,F,10*log10(abs(P)),'EdgeColor','none');   
    axis xy; axis tight; colormap(jet); view(0,90); set(gca,'xtick',[]);
%     set(p,'YData',F,'XData',T,'ZData',Z);
%     axis xy; axis tight;
  end
end
set(handles.fig,'userdata',handles);


function redraw_axes
handles = get(gcf,'userdata');
set(handles.fig,'userdata',handles);
    y = alldata.epochs.data(k,I); %get(p,'YData');
      if any(strcmp('setleft',{'setleft','setright'}))
        set(p,'XData',alldata.epochs.time(I),'YData',y);
      else
        set(p,'YData',y);
      end
      set(h,'ylim',handles.properties.ylim{k},'xlim',[t1 t2]);
      
function onclick(src,evnt)
handles = get(gcf,'userdata');
P  = get(gca,'Currentpoint');
if isfield(handles,'vline') && ~isempty(handles.vline)
  try
    delete(handles.vline); handles.vline=[];
    delete(handles.hline); handles.hline=[];
  end
end
if get(handles.topochk,'value')==1
  for k = 1:length(handles.topoaxes)
    handles.vline(k) = line([P(1,1) P(1,1)],ylim,'color','k','parent',handles.topoaxes(k)); % vline
    handles.hline(k) = line(xlim,[P(1,2) P(1,2)],'color','k','parent',handles.topoaxes(k)); % hline
  end  
  if length(handles.vline)>length(handles.topoaxes)
    handles.vline(k+1:end) = [];
  end
  if length(handles.hline)>length(handles.topoaxes)
    handles.hline(k+1:end) = [];
  end
else
  for k = 1:length(handles.subplot)
    handles.vline(k) = line([P(1,1) P(1,1)],ylim,'color','k','parent',handles.subplot(k)); % vline
    handles.hline(k) = line(xlim,[P(1,2) P(1,2)],'color','g','parent',handles.subplot(k)); % hline
  end
end
set(gcf,'userdata',handles);

function ZoomFunction(src,evnt)
handles = get(gcf,'userdata');
if get(handles.topochk,'value')==1
  CLIM = get(gca,'clim');
  if evnt.VerticalScrollCount < 0           % zoom in
    set(handles.topoaxes,'clim',CLIM/1.5);
  else                                      % zoom out
    set(handles.topoaxes,'clim',CLIM*1.5);
  end
  handles.properties.clim = get(gca,'clim');
else
  YLIM = get(handles.subplot(1),'ylim');
  if evnt.VerticalScrollCount < 0           % zoom in
    set(handles.subplot,'ylim',YLIM/1.5);
  else                                      % zoom out
    set(handles.subplot,'ylim',YLIM*1.5);
  end
  [handles.properties.ylim{:}] = deal(get(handles.subplot(1),'ylim'));
end
set(gcf,'userdata',handles);

function [i] = nearest(array, val)
array = array(:);
if isnan(val), error('incorrect value'); end
if val>max(array) % return the last occurence of the nearest number  
  [dum, i] = max(flipud(array));
  i = length(array) + 1 - i;
else % return the first occurence of the nearest number  
  [mindist, i] = min(abs(array(:) - val));
end

function varargout = vline(x, varargin)
abc = axis;
x = [x x];
y = abc([3 4]);
if length(varargin)==1
  varargin = {'color', varargin{1}};
end
h = line(x, y);
set(h, varargin{:});
if nargout > 0
  varargout{1} = h;
end

% function draw_topoplot
% topo(tmptmp,'avgovertime','yes','zlim','maxmin','nrows',1,'ncols',1,'layout',handles.layout.file,'iEEG_flag',1,'showlabels','yes','fig',0); 

% handles = get(gcf,'userdata');
% X = handles.layout.xpos;
% Y = handles.layout.ypos;
%   if any(X(:)<0), X = X - min(X(:)); end
%   if any(Y(:)<0), Y = Y - min(Y(:)); end
%   
% cfg.normal = setdiff(1:length(X),cfg.highlight);
% % Scale the data to a circle with x-axis and y-axis: -0.45 to 0.45
% y = 0.9*((X-min(X))/(max(X)-min(X))-0.5); % NOTE: x becomes y and y becomes x, in griddata is also reversed which makes x x and y y again
% x = 0.9*((Y-min(Y))/(max(Y)-min(Y))-0.5);
% interplimits = cfg.interplimits;
%   % 'electrodes' to furthest electrode
%   % 'head' to edge of head
% % Find limits for interpolation:
% if strcmp(interplimits,'head') || strcmp(interplimits,'grid')
%   xmin = min(-.5,min(x)); xmax = max(0.5,max(x));
%   ymin = min(-.5,min(y)); ymax = max(0.5,max(y));
% else
%   xmin = max(-.5,min(x)); xmax = min(0.5,max(x));
%   ymin = max(-.5,min(y)); ymax = min(0.5,max(y));
% end
% if ischar(cfg.zlim)
%   if strcmp(cfg.zlim,'absmax')
%     zmax     = max(abs(data.(datafield).data(:)));
%     zmin     = -zmax;
%   elseif strcmp(cfg.zlim,'maxmin')
%     zmin     = min(data.(datafield).data(:));
%     zmax     = max(data.(datafield).data(:));
%   end
% else
%   zmin     = min(data.(datafield).data(:));
%   zmax     = max(data.(datafield).data(:));  
% end
% gridscale  = cfg.gridscale;                   % resolution
% interp     = cfg.interp;                      % 'linear','cubic','nearest','v4'
% xi         = linspace(xmin,xmax,gridscale);   % x-axis description (row vector)
% yi         = linspace(ymin,ymax,gridscale);   % y-axis description (row vector)
% delta      = xi(2)-xi(1);   
% cfg.labels = {data.sensor_info.label};

  