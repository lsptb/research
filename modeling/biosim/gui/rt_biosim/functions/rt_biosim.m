function rt_biosim(varargin)
% specpath='/space/mdeh3/9/halgdev/projects/jsherfey/model/rt_biosim/specs';
% rt_biosim(specpath);

% ----------------------------------------------------------
% get specification
if nargin>0 && isstruct(varargin{1}) % biosim(spec,...)
  spec = varargin{1};
  if nargin>1, varargin = varargin(2:end); end
elseif nargin>0
  if isstr(varargin{1}) && exist(varargin{1},'file')
    spec = loadspec(varargin{1});
    if nargin>1, varargin = varargin(2:end); end
  elseif isstr(varargin{1}) && exist(varargin{1},'dir')
    specpath=varargin{1};
    if nargin>1 && isstr(varargin{2}) % biosim(fpath,prefix,...)
      spec = loadspec(varargin{1},varargin{2});
      if nargin>2, varargin = varargin(3:end); end
    else
      spec = loadspec(varargin{1});
      if nargin>1, varargin = varargin(2:end); end
    end    
  else
    spec = loadspec(varargin{:});
  end
else
  error('You must supply at least one input.');
end
if isfield(spec,'cells') && ~isfield(spec,'entities')
  spec.entities = spec.cells;
  spec = rmfield(spec,'cells');
end
if ~isfield(spec,'files') || isempty(spec.files)
  DBPATH = '/space/mdeh3/9/halgdev/projects/jsherfey/code/modeler/database';
  d=dir(DBPATH);
  spec.files = {d(cellfun(@(x)any(regexp(x,'.txt$')),{d.name})).name};
  spec.files = cellfun(@(x)fullfile(DBPATH,x),spec.files,'unif',0);
end

% ----------------------------------------------------------
% prepare parameters
if nargin <= 1, varargin = {}; end
parms = mmil_args2parms( varargin, ...
                           {  'Iext',@(t) 0,[],...
                              'timelimits',[0 40],[],...
                              'dsfact',1,[],...
                              'logfid',1,[],...
                              'SOLVER','euler',[],...
                              'dt',.01,[],...
                              'output_list',[],[],...
                              'override',[],[],...
                              'IC',[],[],...
                           }, false);

timelimits = parms.timelimits;
fileID = parms.logfid;
dt = parms.dt;

global pauseflag quitflag tlast t points handles F buffer
%global a b c d I v u VU L % params

pauseflag=-1;
quitflag=0;
tlast=-inf; 
points=[];
handles.img=[];
handles.imgax=[];
buffer = 500;

% ----------------------------------------------------------
% get model
%[model,ic,functions,auxvars,spec,readable,StateIndex] = buildmodel(spec,'logfid',parms.logfid,'override',parms.override);
[model,IC,functions,auxvars,spec,readable,StateIndex] = buildmodel2(spec,'logfid',parms.logfid,'override',parms.override,'dt',parms.dt);
if ~isempty(parms.IC) && numel(parms.IC)==numel(IC)
  IC = parms.IC;
end

% evaluate auxiliary variables (ie., adjacency matrices)
for k = 1:size(auxvars,1)
  %try  % added to catch mask=mask-diag(diag(mask)) when mask is not square
    eval(sprintf('%s = %s;',auxvars{k,1},auxvars{k,2}) );
  %end
end
% evaluate anonymous functions
for k = 1:size(functions,1)
  eval(sprintf('%s = %s;',functions{k,1},functions{k,2}) );
end

F = eval(model);

Elabel = spec.entities(1).label;
Npop=spec.entities(1).multiplicity; 
Npre=sqrt(Npop); Npost=sqrt(Npop); nr=sqrt(Npop);
clims=[-80 50];%[0 1]; %[min(F(dt,IC)) max(F(dt,IC))];
% ----------------------------------------------------------
init_gui(nr,clims,Elabel,model,dt,buffer);

% Simulate
X=IC; t=0; dt=parms.dt; cnt=0; record=zeros(length(IC),buffer);
while quitflag==0
  cnt=cnt+1;
  if get(findobj('tag','speed'),'value')~=0
    pause(0.1*get(findobj('tag','speed'),'value')^1);
  end
 
  p=findobj('tag','pause');
  while pauseflag>0
    pause(0);drawnow; %needed to overcome MATLAB7 bug (found by Gerardo Lafferriere) 
    set(p,'string','resume');
  end;
  set(p,'string','pause');
  
  F = eval(model);
  X = X+dt*F(t,X);
  t = t + dt;
  if cnt<=buffer
    record(:,cnt)=X;
  else
    record = [record(:,2:end) X];
  end
  
  % get output to plot
  e=1;
  var = spec.entities(e).ode_labels{1};
  var = find(cellfun(@(x)isequal(x,var),spec.entities(e).var_list));
  Z = X(var);
  switch Elabel
    case {'LeakyIntegrator','LI'}
      Y = reshape(Z(:),sqrt([Npop Npop]));
    case {'SpikingNeuron','SN'}
      Y = reshape(Z(:),sqrt([Npop Npop]));
      % TODO: add firing rate filter
      % ...
    otherwise
      Y = zeros(sqrt([Npop Npop]));
  end
  set(handles.img,'cdata',Y); xlabel(sprintf('time = %g',t));  
  set(handles.avgtrace,'ydata',mean(record(var,:),1));
  if numel(var)>20
    sel=randperm(numel(var));
    var=var(sel(1:20));
  end
  for k=1:length(var)
    set(handles.alltrace(k),'ydata',record(var(k),:));
  end
  drawnow;
end

%{
in rt_biosim():
  ...
  tlast=-inf; target=[];
  % then numeric integration...
    integrate(step)
    if stateevent, tlast(end+1)=tlastevent (ex. spike in target node or LFP)
    if onclick,    tlast(end+1)=tlastclick; [u,v] = getpts(fig); points=[u,v];
      alt: use ginput(1)
    if callback:optionselect
      arglabel=[Elabel '_' opttype '_args']; % ex) opttype='stim'
      args = eval(arglabel);
      [UPDATE ARGS W/ NEW OPTIONS; user modified val of opt (string)]
      args{find(cellfun(@(x)isequal(x,opt)))+1} = val;
      eval([arglabel '=args;']); % ex) for 'stim': opt='type','dur','amp'
note: opttype's in GUI should match mech file names and each given its own
subpanel in the set-param area of the GUI figure; opt's for a type are keys
with defaults set in the mech file var 'args' passed to auxfuncs (which
set the absolute defaults); user opt values are set w/ edit and radio controls.
%}


% SUBFUNCTIONS

function changeparms(src,evnt,parm,varargin)
global handles
  obj = findobj('tag',parm);
  nv = str2num(get(obj,'string'));
  %imgax=findobj('tag','imgax');
  switch parm
    case 'clims'
      if ~isempty(varargin) && isequal(varargin{1},'autoscale')
        findobj('tag','img');
        cdata=get(findobj('tag','img'),'cdata'); 
        nv = [min(cdata(:)) max(cdata(:))];
        set(obj,'string',num2str(nv));
      end
      caxis(nv);%set(imgax,'zlim',nv);
      set(handles.h1,'ylim',nv);
      set(handles.h2,'ylim',nv);
  end


function init_gui(nr,clims,Elabel,model,dt,buffer)
  % Prepare GUI
  global handles
  T=(0:buffer-1)*dt;
  figNumber = figure(1);
  clf;
  set(figNumber,'NumberTitle','off','doublebuffer','on',...
          'Name','Interactive BioSim',...
          'Units','normalized','toolbar','figure',...
          'Position',[0.05 0.1 0.8 0.8]);
  handles.h1=subplot(4,2,1);
  set(handles.h1,'Position',[0.05 0.75 0.27 0.2])
  handles.avgtrace=line('color','k','LineStyle','-','erase','background','xdata',T,'ydata',zeros(1,buffer),'zdata',[]);
  axis([T(1) T(end) -100 30])
  title('avg activity')
  xlabel('time');

  handles.h2=subplot(4,2,3);
  set(handles.h2,'Position',[0.05 0.5 0.27 0.15])
  for k=1:nr^2
    handles.alltrace(k)=line('color','k','LineStyle','-','erase','background','xdata',T,'ydata',zeros(1,buffer),'zdata',[]);
  end
  axis([T(1) T(end) -40 60])
  title('activity overlay')
  xlabel('time')

  imgax=subplot(2,2,2);
  set(imgax,'Position',[0.4 0.5 0.45 0.45],'tag','imgax')
  img = imagesc(zeros(nr,nr)); axis xy; axis square
  set(img,'tag','img');
  handles.img=img;
  handles.imgax=imgax;
  head = line('color','r','Marker','.','markersize',20,'erase','xor','xdata',[],'ydata',[],'zdata',[]);
  tail=line('color','k','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
  vnull=line('color','b','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
  unull=line('color','b','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
  cnull=line('color','g','LineStyle',':','erase','xor','xdata',[],'ydata',[],'zdata',[]);
  thresh=line('color','r','LineStyle','-','erase','xor','xdata',[],'ydata',[],'zdata',[]);
  axis([1 nr 1 nr]); caxis(clims); colorbar;
  title([Elabel ' activity'])
  xlabel('Nx');ylabel('Ny');

  % parameter controls
  a=nan;b=nan;c=nan;d=nan;I=nan;
  uicontrol('Style','frame', 'Units','normalized', ...
            'Position',[0.27  0.05 0.12 0.37]);
  uicontrol('Style','text', 'Units','normalized', 'HorizontalAlignment','left',...
            'tag','parameters','Position',[0.29  0.38 0.07 0.03],'string','parameters:');
%     uicontrol('Style','text', 'Units','normalized', ...
%               'Position',[0.28  0.2 0.03 0.03],'string','d=');
%     uicontrol('Style','edit', 'Units','normalized', ...
%               'Position',[0.31  0.2 0.05 0.03],...
  uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.28  0.2 0.03 0.03],'string','clims=');
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.31  0.2 0.05 0.03],...
            'string',num2str(clims),'tag','clims','Callback',{@changeparms,'clims'});
  uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.28  0.3 0.03 0.03],'string','Npop=');
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.31  0.3 0.05 0.03],...
            'string',num2str(nr^2),'tag','Npop','Callback',['rt_biosim(specpath,''override'',{''' Elabel ''',''multiplicity'',str2num(get(gcbo,''string''))});']);
  uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.28  0.25 0.03 0.03],'string','c=');
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.31  0.25 0.05 0.03],...
            'string',num2str(c),'tag','c','Callback','izhikevich(''changepars'')');
  uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.28  0.35 0.03 0.03],'string','d=');
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.31  0.35 0.05 0.03],...
            'string',num2str(d),'tag','d','Callback','izhikevich(''changepars'')');
  uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.28  0.15 0.03 0.03],'string','I=');
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.31  0.15 0.05 0.03],...
            'string',num2str(I),'tag','I','Callback','izhikevich(''changepars'')');
  uicontrol('Style','text', 'Units','normalized', ...
            'Position',[0.28  0.1  0.05 0.03],...
            'String','current');
  uicontrol('Style','radio', 'Units','normalized', ...
            'Position',[0.33  0.1  0.03 0.03],...
            'tag','current','value',0,'callback','global drawnull;drawnull=1;');

  % model description
  uicontrol('Style','frame', 'Units','normalized', ...
            'Position',[0.4  0.25 0.16 0.17]);
  uicontrol('Style','text', 'Units','normalized', 'HorizontalAlignment','left',...
            'Position',[0.41  0.26 0.14 0.15],'string',['model: ' model]);

  % % actions
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.4  0.19 0.1 0.05],...
            'String','autoscale','Callback',{@changeparms,'clims','autoscale'});
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.4  0.12 0.1 0.05],...
            'String','E-pulse onclick','Callback', 'global t tlast points F; tlast(end+1)=t; [u,v]=getpts; points=[v,u];');
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.4  0.05 0.1 0.05],...
            'String','E-pulse onrect','Callback','global t tlast points F; tlast(end+1)=t; points=getrect;');
  % uicontrol('Style','pushbutton', 'Units','normalized', ...
  %           'Position',[0.57  0.37 0.1 0.05],...
  %           'tag','threshold','String','show threshold','Callback','izhikevich(''showthreshold'')');

  % % slider control
  % uicontrol('Style','frame', 'Units','normalized', ...
  %           'Position',[0.52  0.12 0.31 0.05]);
  % uicontrol('Style','text', 'Units','normalized',...
  %           'Position',[0.525  0.145 0.3 0.02],'string','input current');
  % uicontrol('Style','text', 'Units','normalized','HorizontalAlignment','left',...
  %           'Position',[0.525  0.145 0.05 0.02],'string','-100');
  % uicontrol('Style','text', 'Units','normalized','HorizontalAlignment','right',...
  %           'Position',[0.775  0.145 0.05 0.02],'string','+100');
  % uicontrol('Style','slider', 'Units','normalized', ...
  %           'Position',[0.525  0.125 0.3 0.015],...
  %           'min',-100,'max',100,...
  %           'value',I, 'tag','inputcurrent','callback','izhikevich(''changepars(1)'')');

  % % slider control
  uicontrol('Style','frame', 'Units','normalized', ...
            'Position',[0.52  0.05 0.31 0.05]);
  uicontrol('Style','text', 'Units','normalized',...
            'Position',[0.525  0.075 0.3 0.02],'string','visualization speed');
  uicontrol('Style','slider', 'Units','normalized', ...
            'Position',[0.525  0.055 0.3 0.015],...
            'value',0.0, 'tag','speed');

  % pause/quit
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.85  0.12 0.1 0.05],...
            'String','pause','tag','pause','Callback','global pauseflag;pauseflag=-pauseflag;');
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.85  0.05 0.1 0.05],...
            'String','quit','Callback','global quitflag;quitflag=1;');
