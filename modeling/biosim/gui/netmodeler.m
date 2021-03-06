function netmodeler(net)
global cfg H CURRSPEC LASTSPEC
if nargin<1, net=[]; end
CURRSPEC = net;
LASTSPEC = net;
%cfg = [];
cfg.focusconn = 1;
cfg.pauseflag = -1;
cfg.quitflag = -1;
cfg.tlast=-inf; 
cfg.buffer = 20000;%10000;
cfg.dt = .01;

if ~isfield(net,'cells')
  if isfield(net,'entities')
    net.cells = net.entities;
  else
    net.cells = net;
  end
end
if isfield(net.cells,'files'), net.files = net.cells(1).files; end
if ~isfield(net.cells,'parent'), net.cells(1).parent=[]; end

% get list of all known mechs (stored in DB)
% TODO: read list fom MySQL DB (see http://introdeebee.wordpress.com/2013/02/22/connecting-matlab-to-mysql-database-using-odbc-open-database-connector-for-windows-7/)
DBPATH = '/space/mdeh3/9/halgdev/projects/jsherfey/code/modeler/database';
if ~exist(DBPATH,'dir')
  DBPATH = 'C:\Users\jsherfey\Desktop\My World\Code\modelers\database';
end
[allmechlist,allmechfiles]=get_mechlist(DBPATH);
% use stored mechs if user did not provide list of mech files
if ~isfield(net,'files') || isempty(net.files)
  selmechlist=allmechlist;
  selmechfiles=allmechfiles;
elseif ischar(net.files) && exist(net.files,'dir')
  % user provided a directory that contains the user mech files
  d=dir(net.files);
  selmechlist = {d(cellfun(@(x)any(regexp(x,'.txt$')),{d.name})).name};
  allmechlist = {selmechlist{:} allmechlist{:}}; % prepend user mechs
  selmechfiles = cellfun(@(x)fullfile(DBPATH,x),net.files,'unif',0);
  allmechfiles = {selmechfiles{:} allmechfiles{:}};
elseif iscellstr(net.files)
  % todo: get mechs from spec
%   selmechlist = {d(cellfun(@(x)any(regexp(x,'.txt$')),{d.name})).name};
%   allmechlist = {selmechlist{:} allmechlist{:}}; % prepend user mechs
%   selmechfiles = cellfun(@(x)fullfile(DBPATH,x),net.files,'unif',0);
%   allmechfiles = {selmechfiles{:} allmechfiles{:}};  
end

% only keep mech selection matching those included in the cell
if isfield(net.cells,'mechanisms')
  cellmechs={net.cells.mechanisms}; cellmechs=[cellmechs{:}];
  connmechs={net.connections.mechanisms}; 
  connmechs=connmechs(~cellfun(@isempty,connmechs));
  cellmechs={cellmechs{:} connmechs{:}};
  selmechlist = cellfun(@(x)strrep(x,'.txt',''),selmechlist,'unif',0);
  sel=cellfun(@(x)any(strmatch(x,cellmechs,'exact')),selmechlist);
  selmechlist = selmechlist(sel);
  selmechfiles = selmechfiles(sel);
  net.files = selmechfiles;
  updatemodel(net); % prepare model
end
cfg.selmechlist = selmechlist;
cfg.allmechlist = allmechlist;
cfg.allmechfiles = allmechfiles;
cfg.focuscolor = [.7 .7 .7];
% load all mech data
global allmechs
for i=1:length(allmechfiles)
  this = parse_mech_spec(allmechfiles{i},[]);
  [fpath,fname,fext]=fileparts(allmechfiles{i});
  this.label = fname;
  this.file = allmechfiles{i};
  if i==1, allmechs = this;
  else allmechs(i) = this;
  end
end

%% Create main figure
maxcomp = 5; c=1.5;
%H = [];
sz = get(0,'ScreenSize'); 
H.f_net = figure('position',[.005*sz(3) .01*sz(4) .94*sz(3) .88*sz(4)],...%'position',[125 100 1400 800],...
  'WindowScrollWheelFcn',@ZoomFunction,'CloseRequestFcn','delete(gcf); clear global H');
set(H.f_net,'MenuBar','none');
file_m = uimenu(H.f_net,'Label','File');
uimenu(file_m,'Label','Load sim_data','Callback',@Load_Data);
uimenu(file_m,'Label','Load spec','Callback',@Load_Spec);
uimenu(file_m,'Label','Save spec','Callback',@Save_Spec);
uimenu(file_m,'Label','Grab spec','Callback','global CURRSPEC; assignin(''base'',''spec'',CURRSPEC);');
uimenu(file_m,'Label','Update spec from base','Callback',{@refresh,1});%'if ismember(''spec'',evalin(''base'',''who'')), updatemodel(evalin(''base'',''spec'')); SelectCells; end');
uimenu(file_m,'Label','Exit','Callback','global CURRSPEC H cfg; close(H.f_net); clear CURRSPEC H cfg;');
plot_m = uimenu(H.f_net,'Label','Plot');
condition='global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), ';
uimenu(plot_m,'Label','plotv','Callback',[condition 'plotv(evalin(''base'',''sim_data''),CURRSPEC); else disp(''load data to plot''); end']);
uimenu(plot_m,'Label','plotpow','Callback',[condition 'plotpow(evalin(''base'',''sim_data''),CURRSPEC); else disp(''load data to plot''); end']);
uimenu(plot_m,'Label','plotspk','Callback',[condition 'plotspk(evalin(''base'',''sim_data''),CURRSPEC,''window_size'',30/1000,''dW'',5/1000); else disp(''load data to plot''); end']);

%jobj=findjobj(H.f_net); 
%set(jobj,'MouseEnteredCallback','global H; figure(H.f_net)');
% Panels
H.p_net_select = uipanel('parent',H.f_net,'Position',[.02 .67 .35 .3],'BackgroundColor','white','BorderWidth',.2,'BorderType','line'); % cell morphology
H.p_net_connect = uipanel('parent',H.f_net,'Position',[.02 .44 .35 .34/c],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','Population connections'); % cell specification
H.p_net_kernel = uipanel('parent',H.f_net,'Position',[.02 .01 .35 .42],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','Cell connections'); % cell specification
H.p_state_plot = uipanel('parent',H.f_net,'Position',[.4 .01 .59 .97],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','plots'); % cell specification
% basic controls (load cell model; undo; ...)
%H.btn_loadcell = uicontrol('parent',H.f_net,'units','normalized',...
%  'style','pushbutton','fontsize',10,'string','load','callback',[],...
%  'position',[.02+.5*.2-.01-.06-.03 .97 .05 .03]);%,'BackgroundColor','white');
H.btn_undo = uicontrol('parent',H.f_net,'units','normalized','position',[.02+.5*.2+.05-.06 .97 .05 .03],'style','pushbutton','string','undo','fontsize',10,'callback',@undo);
H.btn_print = uicontrol('parent',H.f_net,'units','normalized','position',[.02+.5*.2+.05 .97 .05 .03],'style','pushbutton','string','print','fontsize',10,'callback',@printmodel);
H.btn_update = uicontrol('parent',H.f_net,'units','normalized','position',[.02+.5*.2+.05+.06 .97 .05 .03],'style','pushbutton','string','refresh','fontsize',10,'callback',{@refresh,0});
% Selection panel
H.txt_lstlabel = uicontrol('parent',H.p_net_select,'units','normalized',...
  'style','text','position',[.02 .87 .25 .06],'string','Compartments','ListboxTop',0);%,'HorizontalAlignment','left','backgroundcolor','w');
if isfield(net.cells,'label')
  l={net.cells.label}; i=1:length(l);
else
  l={}; i=[];
end
H.lst_comps = uicontrol('parent',H.p_net_select,'units','normalized',...
  'style','listbox','position',[.02 .05 .25 .8],'value',i,'string',l,...
  'backgroundcolor','w','Max',maxcomp,'Min',0,'Callback',@SelectCells);

if ~isfield(net.cells,'mechanisms'), return; end

% Selection and Connection panels
SelectCells;
DrawAuxView;
DrawSimPlots;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_Data(src,evnt)
[filename,pathname] = uigetfile({'*.mat'},'Pick a file.','MultiSelect','off');
if isequal(filename,0) || isequal(pathname,0), return; end
if iscell(filename)
  datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
  filename = filename{1};
else
  datafile = [pathname filename];
end
if exist(datafile,'file')
  fprintf('Loading data: %s\n',datafile);
  o=load(datafile);
  if isfield(o,'sim_data') && isfield(o,'spec')
    assignin('base','sim_data',o.sim_data);
    assignin('base','spec',o.spec);
    global CURRSPEC
    CURRSPEC = o.spec;
    refresh(src,evnt);
%     global H
%     if isfield(H,'f_cell') && ishandle(H.f_cell), close(H.f_cell); end
%     if isfield(H,'f_net') && ishandle(H.f_net), close(H.f_net); end
%     netmodeler(o.spec);
  else
    fprintf('select file does not contain sim_data and spec. no data loaded\n');
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_Spec(src,evnt)
[filename,pathname] = uigetfile({'*.mat'},'Pick a file.','MultiSelect','off');
if isequal(filename,0) || isequal(pathname,0), return; end
if iscell(filename)
  datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
  filename = filename{1};
else
  datafile = [pathname filename];
end
if exist(datafile,'file')
  fprintf('Loading model specification: %s\n',datafile);
  o=load(datafile);
  if isfield(o,'spec')
    global CURRSPEC
    CURRSPEC = o.spec;
    refresh(src,evnt);
  else
    fprintf('select file does not contain spec structure. no specification loaded\n');
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Save_Spec(src,evnt)
[filename,pathname] = uiputfile({'*.mat;'},'Save as','model_specification.mat');
if isequal(filename,0) || isequal(pathname,0)
  return;
end
outfile = fullfile(pathname,filename);
[fpath,fname,fext] = fileparts(outfile);
global CURRSPEC
spec=CURRSPEC;
fprintf('Saving model specification: %s\n',outfile);
save(outfile,'spec');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SelectCells(src,evnt)
global H CURRSPEC
v=get(H.lst_comps,'value'); 
l=get(H.lst_comps,'string');
set(H.lst_comps,'string',{CURRSPEC.cells.label},'value',v(v<=length(CURRSPEC.cells)));
DrawCellInfo(CURRSPEC); % Selection panel
DrawNetGrid(CURRSPEC);  % Connection panel

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawCellInfo(net)
global H
c=1.5; dy=-.07*c; 
sel = get(H.lst_comps,'value');
l={net.cells(sel).label}; 
N=[net.cells(sel).multiplicity];
mechs={net.cells(sel).mechanisms};
for i=1:length(sel)
  m=mechs{i};
  str=m{1}; for j=2:length(m), str=[str ', ' m{j}]; end
  if ~isfield(H,'edit_comp_label') || length(H.edit_comp_label)<length(sel) || ~ishandle(H.edit_comp_label(i))
    H.edit_comp_label(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.33 .8+dy*(i-1) .1 .06],'backgroundcolor','w','string',l{i},...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'label'});
    H.edit_comp_N(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.44 .8+dy*(i-1) .06 .06],'backgroundcolor','w','string',N(i),...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'multiplicity'});
    H.edit_comp_mechs(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.51 .8+dy*(i-1) .42 .06],'backgroundcolor','w','string',str,...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'mechanisms'});
    H.btn_comp_delete(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','-','callback',{@DeleteCell,l{i}},...
      'position',[.295 .8+dy*(i-1) .03 .06]);%,'BackgroundColor','white');    
    H.btn_comp_copy(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','+','callback',{@CopyCell,l{i}},...
      'position',[.93 .8+dy*(i-1) .03 .06]);%,'BackgroundColor','white');    
    H.btn_comp_edit(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','o','callback',{@OpenCellModeler,l{i}},...
      'position',[.965 .8+dy*(i-1) .03 .06]);%,'BackgroundColor','white');        
  else
    % update properties
    set(H.edit_comp_label(i),'string',l{i},'visible','on','Callback',{@UpdateCells,l{i},'label'});
    set(H.edit_comp_N(i),'string',N(i),'visible','on','Callback',{@UpdateCells,l{i},'multiplicity'});
    set(H.edit_comp_mechs(i),'string',str,'visible','on','Callback',{@UpdateCells,l{i},'mechanisms'});
    set(H.btn_comp_copy(i),'callback',{@CopyCell,l{i}},'visible','on');
    set(H.btn_comp_delete(i),'callback',{@DeleteCell,l{i}},'visible','on');
    set(H.btn_comp_edit(i),'callback',{@OpenCellModeler,l{i}},'visible','on');
  end
  if length(H.edit_comp_label)>i
    set(H.edit_comp_label(i+1:end),'visible','off');
    set(H.edit_comp_N(i+1:end),'visible','off');
    set(H.edit_comp_mechs(i+1:end),'visible','off');
    set(H.btn_comp_copy(i+1:end),'visible','off');
    set(H.btn_comp_delete(i+1:end),'visible','off');
    set(H.btn_comp_edit(i+1:end),'visible','off');    
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawNetGrid(net)%,forceflag)
global H
dx=.15; x=.13; c=1.5; dy=-.07*c; 
sel = get(H.lst_comps,'value');
net.cells = net.cells(sel);
net.connections = net.connections(sel,:);
net.connections = net.connections(:,sel);
l={net.cells.label}; 
for i=1:length(sel)
  for j=1:length(sel)
    m = net.connections(i,j).mechanisms;
    pos = [x+dx*(i-1) .8+dy*(j-1) .9*dx .06*c];
    u.from=l{j}; u.to=l{i};
    connname = net.connections(i,j).label;
    if ~isempty(m)
      if ischar(m), m={m}; end
      str=m{1}; for k=2:length(m), str=[str ', ' m{k}]; end
    else
      str = '';
    end
    if ~isfield(H,'txt_to') || i>length(H.txt_from) || j>length(H.txt_to) || ~ishandle(H.edit_conn_mechs(i,j)) || H.edit_conn_mechs(i,j)==0
      if i==1
        H.txt_to(j) = uicontrol('parent',H.p_net_connect,'units','normalized',...
          'style','text','position',[x+dx*(j-1) .91 .11 .06*c],'string',l{j});%,'HorizontalAlignment','left','backgroundcolor','w');
      end
      if j==1
        H.txt_from(i) = uicontrol('parent',H.p_net_connect,'units','normalized',...
          'style','text','position',[.01 .8+dy*(i-1) .11 .06*c],'string',l{i});%,'HorizontalAlignment','left','backgroundcolor','w');      
      end
      H.edit_conn_mechs(i,j) = uicontrol('parent',H.p_net_connect,'units','normalized',...
        'style','edit','position',pos,'backgroundcolor','w',...
        'string',str,'HorizontalAlignment','left','UserData',u);
      set(H.edit_conn_mechs(i,j),'ButtonDownFcn',...
        {@Display_Mech_Info,pos,H.edit_conn_mechs(i,j),connname},...
        'Callback',{@UpdateNet,connname})
    else
      set(H.txt_to(i),'string',l{i},'visible','on');
      set(H.txt_from(i),'string',l{i},'visible','on');
      set(H.edit_conn_mechs(i,j),'string',str,'UserData',u,'Callback',{@UpdateNet,connname},'ButtonDownFcn',{@Display_Mech_Info,pos,H.edit_conn_mechs(i,j),connname},'visible','on');
    end
  end
end
if isfield(H,'txt_to') && length(H.txt_to)>i
  set(H.txt_to(i+1:end),'visible','off');
  set(H.txt_from(i+1:end),'visible','off');
  set(H.edit_conn_mechs(i+1:end,:),'visible','off');
  set(H.edit_conn_mechs(:,i+1:end),'visible','off');
end    

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawAuxView
global H CURRSPEC
% evaluate auxiliary variables
a = CURRSPEC.model.auxvars;
for i=1:size(a,1)
  key = a{i,1};
  try
    eval(sprintf('%s = %s;',a{i,1},a{i,2}));
    val = eval(a{i,2});
  catch
    val = nan;
  end
  CURRSPEC.model.eval.(key) = val;
end
% make list of connections
l={CURRSPEC.cells.label};
ll={CURRSPEC.connections.label};
sel=find(~cellfun(@isempty,ll));
if isempty(sel), return; end
lst={}; % = list of connections: from-to.mechanism
for i=1:length(sel)
  m=CURRSPEC.connections(sel(i)).mechanisms;
  if ~iscell(m), m={m}; end
  for j=1:length(m)
    lst{end+1} = [ll{sel(i)} '.' m{j}];
  end
  if i==1
    [from,to]=ind2sub(size(CURRSPEC.connections),sel(i));
    prefix = [l{from} '_' l{to}];
  end
end
% get list of auxvars for the first connection
sel=strmatch([prefix '_'],CURRSPEC.model.auxvars(:,1));
auxlist=CURRSPEC.model.auxvars(sel,1);
% get expression for the last auxvar
a=CURRSPEC.connections(from,to).mechs(1).auxvars(end,:);
auxeqn=[a{1} ' = ' a{2}];
% get list of params for the first connection
parmlist=fieldnames(CURRSPEC.connections(from,to).mechs(1).params);
% get value of first parameter
parmval=CURRSPEC.connections(from,to).mechs(1).params.(parmlist{1});
% get auxvar matrix 
auxmat=CURRSPEC.model.eval.(auxlist{end});
midpt=round(size(auxmat,2)/2); lims=[.9*min(auxmat(:)) 1.1*max(auxmat(:))];

% dropdown box to select the connection to examine
H.pop_sel_conn = uicontrol('style','popupmenu','string',lst,'parent',H.p_net_kernel,...
  'units','normalized','position',[.1 .94 .33 .06],'value',1);
% listbox to select which auxvar to plot
H.lst_auxvars = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','listbox','position',[.1 .1 .3 .3],'string',auxlist,'value',length(auxlist),...
  'Callback',@SelectAuxVar); % 'backgroundcolor','w',
% edit box to modify defining expressions
H.edit_auxvar_eqn = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','edit','position',[.1 .03 .3 .06],'backgroundcolor','w','string',auxeqn,...
  'HorizontalAlignment','left');%,'Callback',{@UpdateParams,'change','cells'});
% listbox to select a parameter to modify
H.lst_auxparams = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','listbox','position',[.44 .1 .3 .3],'string',parmlist,...
  'Callback',@SelectAuxVar); % 'backgroundcolor','w',
% edit box to modify parameter
H.edit_auxparams = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','edit','position',[.44 .03 .3 .06],'backgroundcolor','w','string',num2str(parmval),...
  'HorizontalAlignment','left');%,'Callback',{@UpdateParams,'change','cells'});

% plots
H.ax_conn_img = subplot('position',[.1 .5 .3 .4],'parent',H.p_net_kernel); 
H.img_connect = imagesc(auxmat); axis xy; vline(midpt,'k'); caxis(lims);

% callback should use: ginput() - select N points interactively...

H.ax_conn_line = subplot('position',[.44 .5 .3 .4],'parent',H.p_net_kernel); 
H.line_connect = line('ydata',1:size(auxmat,1),'xdata',auxmat(:,midpt),'color','k','LineStyle','-','erase','background');
xlim(lims);
% create button group for predefined adjacency matrices
H.rad_adj = uibuttongroup('visible','off','units','normalized','Position',[.77 .5 .2 .3],'parent',H.p_net_kernel);
H.rad_adj_1 = uicontrol('Style','radiobutton','String','one-to-one','units','normalized',...
    'pos',[.05 .75 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
H.rad_adj_2 = uicontrol('Style','radiobutton','String','all-to-all','units','normalized',...
    'pos',[.05 .45 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
H.rad_adj_3 = uicontrol('Style','radiobutton','String','random','units','normalized',...
    'pos',[.05 .15 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
set(H.rad_adj,'SelectionChangeFcn',@seladj);
set(H.rad_adj,'SelectedObject',[]);  % No selection
set(H.rad_adj,'Visible','on');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulate(src,evnt,action)
global CURRSPEC H cfg
DrawSimPlots;
functions = CURRSPEC.model.functions;
auxvars = CURRSPEC.model.auxvars;
ode = CURRSPEC.model.ode;
IC = CURRSPEC.model.IC;
cfg.T = (0:cfg.buffer-1)*cfg.dt; % ms
fh_pow = @(x,y) PowerSpecTA(x,y,[10 80],min(8000,cfg.buffer/2),'Normalized',[]);

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

F = eval(ode);

% Simulate (simple forward Euler integration)
X=IC; t=0; cnt=0; cfg.record=zeros(length(IC),cfg.buffer);
while cfg.quitflag<0
  cnt=cnt+1;
  % speed?
  if get(findobj('tag','speed'),'value')~=0
    pause(0.1*get(findobj('tag','speed'),'value')^1);
  end
  % pause?
  p=findobj('tag','pause');
  while cfg.pauseflag>0
    pause(0);drawnow; %needed to overcome MATLAB7 bug (found by Gerardo Lafferriere) 
    set(p,'string','resume');
  end;
  set(p,'string','pause');
  % integrate
  X = X+cfg.dt*F(t,X);
  t = t + cfg.dt;
  if cnt<=cfg.buffer
    cfg.record(:,cnt)=X;
  else
    cfg.record = [cfg.record(:,2:end) X];
    %if mod(cnt,100)==0
      %set(H.ax_state_plot,'xtick',get(H.ax_state_plot(1),'xtick')+100*cfg.dt);
    %end
  end
  % get output to plot
  list=get(H.lst_comps,'string');
  sel=get(H.lst_comps,'value');
  if length(sel)>3, sel=sel(1:3); end
  show=list(sel);
  for k=1:length(show)
    this=sel(k);
    numcell=CURRSPEC.cells(this).multiplicity;
    if numcell > cfg.ncellshow
      tmp=randperm(numcell);
      inds = sort(tmp(1:cfg.ncellshow));
    else
      inds = 1:numcell;
    end
    var = CURRSPEC.cells(this).ode_labels{1};
    var = find(cellfun(@(x)isequal(x,var),CURRSPEC.variables.labels));%CURRSPEC.cells(this).var_list));
    for j=1:length(inds)
      set(H.simdat_alltrace(k,j),'ydata',cfg.record(var(inds(j)),:));
    end    
    if 0 
      lfp=mean(cfg.record(var,:),1);
      set(H.simdat_LFP(k),'ydata',lfp,'linewidth',4,'color','k','linestyle','-');
    end
    if 0 && mod(cnt,cfg.buffer)==0 && cnt>1
      % update power spectrum
      try
        res = feval(fh_pow,cfg.T,lfp);
        set(H.simdat_LFP_power(k),'xdata',res.f,'ydata',log10(res.Pxx));
      catch
        fprintf('power spectrum calc failed.\n');
      end
    end
  end
  drawnow
end
cfg.quitflag=-1;
p=findobj('tag','pause');
set(p,'string','pause'); 
set(findobj('tag','start'),'string','restart');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawSimPlots
global H cfg CURRSPEC
cfg.T = (0:cfg.buffer-1)*cfg.dt; % ms
cfg.f = 1:100;
cfg.ncellshow=10;
cfg.colors  = 'kbrgmy';
cfg.lntype  = {'-',':','-.','--'};
dy=-.25;
list=get(H.lst_comps,'string');
sel=get(H.lst_comps,'value');
show=list(sel);
if length(show)>3, show=show(1:3); end
% data plots
if isfield(H,'ax_state_plot')
  delete(H.simdat_alltrace(ishandle(H.simdat_alltrace)));
  delete(H.simdat_LFP(ishandle(H.simdat_LFP)));
  delete(H.ax_state_plot(ishandle(H.ax_state_plot)));
  delete(H.simdat_LFP_power(ishandle(H.simdat_LFP_power)));
  delete(H.ax_state_power(ishandle(H.ax_state_power)));
  H=rmfield(H,{'simdat_alltrace','simdat_LFP','ax_state_plot','simdat_LFP_power','ax_state_power'});
end
for i=1:length(show) % i=1:ncomp
  H.ax_state_plot(i) = subplot('position',[.05 .7+(i-1)*dy .6 -.8*dy],...
    'parent',H.p_state_plot,'linewidth',3,'color','w'); 
  for k=1:cfg.ncellshow
    H.simdat_alltrace(i,k)=line('color',cfg.colors(max(1,mod(k,length(cfg.colors)))),'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',cfg.T,'ydata',zeros(1,cfg.buffer),'zdata',[]);
  end
  H.simdat_LFP(i)=line('color',cfg.colors(max(1,mod(k,length(cfg.colors)))),'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',cfg.T,'ydata',zeros(1,cfg.buffer),'zdata',[],'linewidth',2);
  axis([cfg.T(1) cfg.T(end) -100 30]); xlabel('time (ms)');
  title([show{i} '.V']);  %ylabel([show{i} '.V']);  
  H.ax_state_power(i) = subplot('position',[.7 .7+(i-1)*dy .25 -.8*dy],'parent',H.p_state_plot); 
  H.simdat_LFP_power(i)=line('color','k','LineStyle','-','erase','background','xdata',cfg.f,'ydata',zeros(size(cfg.f)),'zdata',[],'linewidth',2);
  %H.ax_state_img(i) = subplot('position',[.7 .7+(i-1)*dy .25 -.8*dy],'parent',H.p_state_plot); 
  %H.img_state = imagesc(zeros(10,40)); axis xy; axis square
end
% % slider control
if isempty(findobj('tag','speed'))
  dx=-.1;
  uicontrol('Style','frame', 'Units','normalized', ...
            'Position',[0.52+dx  0.05 0.31 0.05]);
  uicontrol('Style','text', 'Units','normalized',...
            'Position',[0.525+dx  0.075 0.3 0.02],'string','visualization speed');
  uicontrol('Style','slider', 'Units','normalized', ...
            'Position',[0.525+dx  0.055 0.3 0.015],...
            'value',0.0, 'tag','speed');
end
if isempty(findobj('tag','start'))
  % btn: start <=> reset        
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.75  0.05 0.04 0.05],...
            'String','start','tag','start','Callback',{@simulate,'restart'}); % start <=> pause
end
if isempty(findobj('tag','pause'))
  % btn: pause <=> resume
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.80  0.05 0.04 0.05],...
            'String','pause','tag','pause','Callback','global cfg;cfg.pauseflag;cfg.pauseflag=-cfg.pauseflag;');
end
if isempty(findobj('tag','stop'))
  % btn: stop
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.85  0.05 0.04 0.05],...
            'String','stop','tag','stop','Callback','global cfg;cfg.quitflag=1;');
end     
if isempty(findobj('tag','dt'))
  % btn: start <=> reset        
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.75  0.11 0.04 0.05],...
            'String',num2str(cfg.dt),'tag','dt','Callback','global cfg; cfg.dt=str2num(get(gcbo,''string''));'); % start <=> pause
end        
if isempty(findobj('tag','buffer'))
  % btn: start <=> reset        
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.8  0.11 0.04 0.05],...
            'String',num2str(cfg.buffer),'tag','buffer','Callback','global cfg; cfg.buffer=str2num(get(gcbo,''string''));'); % start <=> pause
end     

% autoscale       
uicontrol('Style','pushbutton', 'Units','normalized', ...
          'Position',[0.9 0.05 0.075 0.05],...
          'String','autoscale','Callback',{@setlimits,'autoscale'});
% simstudy       
uicontrol('Style','pushbutton', 'Units','normalized', ...
          'Position',[0.9 0.11 0.075 0.05],...
          'String','sim study','Callback','StudyDriverUI;');
        
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits(src,evnt,action)
global cfg H
if ~isfield(cfg,'record'), return; end
switch action
  case 'autoscale'
    ymin = min(cfg.record(:));
    ymax = max(cfg.record(:));
    if ymin~=ymax
      set(H.ax_state_plot,'ylim',[ymin ymax]);
    end
end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OpenCellModeler(src,evnt,compname)
global CURRSPEC
ind = strmatch(compname,{CURRSPEC.cells.label},'exact');
if isfield(CURRSPEC.cells,'parent') && ischar(CURRSPEC.cells(ind).parent)
  ind = find(cellfun(@(x)isequal(CURRSPEC.cells(ind).parent,x),{CURRSPEC.cells.parent}));
end
cellmodeler(CURRSPEC.cells(ind));
%keyboard

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CopyCell(src,evnt,compname)
global CURRSPEC H
n=length(CURRSPEC.cells); 
ind = strmatch(compname,{CURRSPEC.cells.label},'exact');
lold=CURRSPEC.cells(ind).label;
lnew=sprintf('%s%g',lold,n+1); cnt=1;
while any(strmatch(lnew,{CURRSPEC.cells.label},'exact'))
  lnew=sprintf('%s%g',lold,n+1+cnt); cnt=cnt+1;
end
newspec = CURRSPEC;
newspec.cells(n+1) = CURRSPEC.cells(ind);
newspec.cells(n+1).label = lnew;
newspec.connections(n+1,1:n) = newspec.connections(ind,1:n);
newspec.connections(1:n,n+1) = newspec.connections(1:n,ind);
for i=1:n+1
  l=newspec.connections(n+1,i).label;
  if ischar(l), l=strrep(l,lold,lnew); end
  newspec.connections(n+1,i).label=l;
  l=newspec.connections(i,n+1).label;
  if ischar(l), l=strrep(l,lold,lnew); end  
  newspec.connections(i,n+1).label=l;
end
set(H.lst_comps,'string',{newspec.cells.label},'value',[get(H.lst_comps,'value') n+1]);
updatemodel(newspec);
SelectCells;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteCell(src,evnt,compname)
global CURRSPEC H
%forceflag = 1;
ind = strmatch(compname,{CURRSPEC.cells.label},'exact');
newspec = CURRSPEC;
newspec.cells(ind) = [];
newspec.connections(ind,:) = [];
newspec.connections(:,ind) = [];
cfg.focusconn = 1;
%cfg.focuscomp = 1;
%cfg.focusmech = 1;
updatemodel(newspec);
SelectCells(src,evnt);%,forceflag);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Display_Mech_Info(src,evnt,mechpos,hmech,connname)
% purpose: display the mech model in a readable form
% try
global H CURRSPEC cfg
uitype = 'edit';
cfg.focusconn = find(cellfun(@(x)isequal(connname,x),{CURRSPEC.connections.label}));
spec = CURRSPEC.connections(cfg.focusconn);
mechstr = get(hmech,'string'); % only show the first mech in list if n>1
if isempty(mechstr), return; end
mechstr = splitstr(mechstr,',');
mechstr = mechstr{1};
focusmech = strmatch(mechstr,spec.mechanisms,'exact');
if iscell(spec.mechanisms)
  mechname = spec.mechanisms{focusmech};
else
  mechname = spec.mechanisms;
end
mech = spec.mechs(focusmech);
showname = [connname '.' mechname];
%Update(src,evnt,'controls',connname,mechname);
fpos=get(H.f_net,'position');
panpos=get(H.p_net_connect,'position'); lo=.5;
pos(1)=fpos(1)+fpos(3)*panpos(1);
pos(2)=fpos(2)+fpos(4)*(panpos(2)-lo);
pos(3)=fpos(3)*panpos(3);
pos(4)=fpos(4)*lo+.86*fpos(4)*mechpos(2)*panpos(4);
if isfield(H,'f_net_hover') && ishandle(H.f_net_hover)
  figure(H.f_net_hover); 
  set(gcf,'name',showname,'position',pos);
else
  H.f_net_hover = figure('MenuBar','none','name',showname,'NumberTitle','off','position',pos);  
end
clf
fontname = 'courier'; % fixedwidth, courier, helvetica (default)
H.mech_txt = uicontrol('style',uitype,'string','','FontName',fontname,'FontSize',10,...
  'BackgroundColor',[.9 .9 .9],'Max',50,'HorizontalAlignment','Left','unit','normalized','Position',[0 0 1 1]);
try
  jobj=findjobj(H.mech_txt);
  set(jobj,'MouseClickedCallback',{@Display_Mech_Info,mechpos,hmech,connname});
end
if strcmp(uitype,'edit')
  set(H.f_net_hover,'MenuBar','none'); 
  mech_apply = uimenu(H.f_net_hover,'Label','Apply','Callback',{@UpdateMech,H.mech_txt,connname,mechname});
  file_write = uimenu(H.f_net_hover,'Label','Write','Callback',{@SaveMech,H.mech_txt,connname,mechname});
  set(H.mech_txt,'BackgroundColor','w');
end
txt = mech_spec2str(mech);
set(H.mech_txt,'string',txt,'position',[0 0 1 1]);
% catch
  %disp('ERROR!');
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMech(src,evnt,htxt,connname,mechname)
global CURRSPEC
% purpose: write mech to new file
outfile = '';
txt = get(htxt,'string');
%thiscell = strmatch(compname,{CURRSPEC.cells.label},'exact');
%comp = CURRSPEC.cells(thiscell);
%thismech = strmatch(mechname,comp.mechanisms,'exact');

% - popup dialog box
% - use fprintf to write each line of txt{:}
% - store new mech name in spec and update the model

% CURRSPEC.files{end+1} = outfile;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateMech(src,evnt,htxt,connname,mechname)
% purpose: apply user changes to the mech model
global CURRSPEC
spec=CURRSPEC;
txt = get(htxt,'string');
newmech = parse_mech_spec(txt);
newmech.label = mechname;
thisconn = find(cellfun(@(x)isequal(connname,x),{spec.connections.label}));
conn = spec.connections(thisconn);
if ~iscell(conn.mechanisms), conn.mechanisms={conn.mechanisms}; end
thismech = strmatch(mechname,conn.mechanisms,'exact');
if ~isempty(thismech)
  conn.mechs(thismech) = newmech;
else
  conn.mechs(end+1) = newmech;
  conn.mechanisms{end+1} = mechname;
end
spec.connections(thisconn) = conn;
updatemodel(spec);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatemodel(newspec) % maybe use same specification as for "override"
% purpose: update the integratable model after each change to its specification
global CURRSPEC LASTSPEC currspec
LASTSPEC = CURRSPEC;
CURRSPEC = newspec;
% check for updated cell model from CellModeler
if isempty(currspec)
  clear global currspec
else
  l={currspec.cells.label}; L={CURRSPEC.cells.label};
  % update compartments in both models
  [s1,s2]=match_str(L,l);
  if ~isempty(s1) && (...
      ~isequal({CURRSPEC.cells(s1).mechs},{currspec.cells(s2).mechs}) || ...
      (isfield(CURRSPEC.cells,'parent') && ~isequal({CURRSPEC.cells(s1).parent},{currspec.cells(s2).parent})) || ...
      ~isequal({CURRSPEC.cells(s1).parameters},{currspec.cells(s2).parameters})...
                     )
    CURRSPEC.cells(s1).mechanisms = currspec.cells(s2).mechanisms;
    CURRSPEC.cells(s1).parameters = currspec.cells(s2).parameters;
    CURRSPEC.cells(s1).dynamics = currspec.cells(s2).dynamics;
    CURRSPEC.cells(s1).mechs = currspec.cells(s2).mechs;
    if isfield(currspec.cells,'parent')
      CURRSPEC.cells(s1).parent = currspec.cells(s2).parent;
    end
    warning('Incorporating external changes to cell model.');
  end
  % add new compartments to the network model
  newcomps = setdiff(l,L);
  if ~isempty(newcomps)
    [s1,s2] = match_str(l,newcomps);
    CURRSPEC.cells = cat(2,CURRSPEC.cells,currspec.cells(s1));
    L={CURRSPEC.cells.label};
    % add connections associated with the new compartments
    CURRSPEC.connections(length(L),length(L)).label=[];
    [I,J]=ind2sub(size(currspec.connections),find(~cellfun(@isempty,{currspec.connections.label})))
    kp = ismember(I,s1) | ismember(J,s1);
    I=I(kp); J=J(kp); % <== these are sub-indices to the new connections in the cell model
    for k=1:length(I)
      ii=strmatch(currspec.cells(I(k)).label,L); % <== these are indices in the network model
      jj=strmatch(currspec.cells(J(k)).label,L);
      CURRSPEC.connections(ii,jj) = currspec.connections(I(k),J(k));
    end    
  end
%   if exist('H','var') && isfield(H,'f_cell') && ishandle(H.f_cell)
%     close(H.f_cell);
%   end
end
if 1
  % remove all connections for which gCOM=0 before building the model
  % ...
  [model,IC,functions,auxvars,CURRSPEC] = buildmodel2(CURRSPEC,'verbose',0);
  try CURRSPEC.cells=CURRSPEC.entities; CURRSPEC=rmfield(CURRSPEC,'entities'); end
  fprintf('model updated successfully\n');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undo(src,evnt)
% revert to the last working model
global LASTSPEC
updatemodel(LASTSPEC);
SelectCells;

function refresh(src,evnt,where)
if nargin<3 || where==0
  global CURRSPEC
  updatemodel(CURRSPEC);
elseif where==1 % get spec from base workspace
  if ismember('spec',evalin('base','who'))
    updatemodel(evalin('base','spec'));
  else
    return;
  end
end
SelectCells;
DrawSimPlots; % ??? should this be here???

function printmodel(src,evnt)
global CURRSPEC
buildmodel2(CURRSPEC,'verbose',1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = mech_spec2str(mech)
% Purpose: prepare text to display mech model parameters and equations
txt = {}; n=0;
% print parameters
keys=fieldnames(mech.params);
vals=struct2cell(mech.params);
for i=1:length(keys)
  if i==1, n=n+1; txt{n}=sprintf('%% Parameters:'); end
  n=n+1; txt{n}=sprintf('%s = %s',keys{i},dat2str(vals{i}));
  if i==length(keys), n=n+1; txt{n}=sprintf(' '); end
end
% print auxiliary functions
for i=1:size(mech.auxvars,1)
  if i==1, n=n+1; txt{n}=sprintf('%% Auxiliary variables:'); end
  n=n+1; txt{n}=sprintf('%s = %s',mech.auxvars{i,1},dat2str(mech.auxvars{i,2}));
  if i==size(mech.auxvars,1), n=n+1; txt{n}=sprintf(''); end
end
% print functions
for i=1:size(mech.functions,1)
  if i==1, n=n+1; txt{n}=sprintf('%% Functions:'); end
  % put in a form that parse_mech_spec() can process
  tmp = sprintf('%s = %s',mech.functions{i,:});
  lhs = regexp(tmp,'^\w+','match'); lhs=lhs{1};
  var = regexp(tmp,'@\([\w,]+\)','match'); var=var{1};
  rhs = regexp(tmp,'@\([\w,]+\).+$','match'); rhs=rhs{1};
  rhs = strtrim(strrep(rhs,var,''));
  var = var(2:end);
  tmp = sprintf('%s%s = %s',lhs,var,rhs);
  n=n+1; txt{n}=tmp;%sprintf('%s = %s',tmp);
  if i==size(mech.functions,1), n=n+1; txt{n}=sprintf(' '); end
end
% print odes
for i=1:size(mech.odes,1)
  if i==1, n=n+1; txt{n}=sprintf('%% ODEs:'); end
  n=n+1; txt{n}=sprintf('%s'' = %s',mech.statevars{i},mech.odes{i});
  n=n+1; txt{n}=sprintf('%s(0) = %s',mech.statevars{i},dat2str(mech.ic{i}));
  if i==size(mech.odes,1), n=n+1; txt{n}=sprintf(' '); end
end
% print interface statements
for i=1:size(mech.substitute,1)
  if i==1, n=n+1; txt{n}=sprintf('%% Substitution:'); end%Expose and/or insert into compartment dynamics:'); end
  n=n+1; txt{n}=sprintf('%s => %s',mech.substitute{i,:});
  if i==size(mech.substitute,1), n=n+1; txt{n}=sprintf(' '); end
end

function val=dat2str(val)
% purpose: convert various data classes into a character form for readable display
if isnumeric(val)
  val = num2str(val);
elseif ischar(val)
  % do nothing
else
  val = 'unrecognized';
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CURRSPEC = addmech(CURRSPEC,mechadded,type,index)
% purpose: add a mechanism to the compartment model
global allmechs cfg
if isempty(mechadded), return; end
if ~iscell(mechadded), mechadded={mechadded}; end
for i=1:length(mechadded)
  newmech = mechadded{i};
  mechind = find(strcmp({allmechs.label},newmech));
  newmech = allmechs(mechind);
  if isempty(newmech)
    warndlg([mechadded{i} ' not found. Check spelling and case.']);
    disp('known mechanisms include: '); disp(get_mechlist');
  else
    newmech = rmfield(newmech,'file');
    CURRSPEC.(type)(index).mechs(end+1)=newmech;
    CURRSPEC.(type)(index).mechanisms{end+1}=mechadded{i};
    CURRSPEC.files{end+1} = cfg.allmechfiles{mechind};
  end
end

function CURRSPEC = removemech(CURRSPEC,mechremoved,type,index)
% purpose: remove a mechanism from the compartment model
global cfg
if isempty(mechremoved), return; end
if ~iscell(mechremoved), mechremoved={mechremoved}; end
for i=1:length(mechremoved)
  oldmech = mechremoved{i};
  oldmech = find(strcmp(CURRSPEC.(type)(index).mechanisms,oldmech));
  if ~isempty(oldmech)   
    CURRSPEC.(type)(index).mechs(oldmech)=[];
    CURRSPEC.(type)(index).mechanisms(oldmech)=[];
  end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateCells(src,evnt,compname,field)
% apply user changes to the compartment mechanisms
global cfg CURRSPEC H 
%forceflag = 1;
this = strmatch(compname,{CURRSPEC.cells.label},'exact');
if numel(this)>1, this = this(1); end
newspec = CURRSPEC;
switch field
  case 'label'
    newspec.cells(this).(field) = get(src,'string');
    updatemodel(newspec);
    set(H.lst_comps,'string',{newspec.cells.label});
    SelectCells;
  case 'multiplicity'
    newspec.cells(this).(field) = str2num(get(src,'string'));
    updatemodel(newspec);
    % DrawKernels;
  case 'mechanisms'
    %pause(.5);
    mechlist = newspec.cells(this).mechanisms;
    str = get(src,'string');
    if ~isempty(str)
      newmechlist = strtrim(splitstr(str,','));
      mechadded = setdiff(newmechlist,mechlist);
      mechremoved = setdiff(mechlist,newmechlist);
      if ~isempty(mechadded) || ~isempty(mechremoved)
        newspec = removemech(newspec,mechremoved,'cells',this);
        newspec = addmech(newspec,mechadded,'cells',this);
        updatemodel(newspec);
      end
    end    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateNet(src,evnt,connname)
% apply user changes to the connection mechanisms
global cfg CURRSPEC
newspec = CURRSPEC; %pause(.5);
if isempty(connname) % first connection b/w these compartments...
  u=get(src,'UserData');
  connname = [u.from '-' u.to];
  from=find(strcmp({newspec.cells.label},u.from));
  to=find(strcmp({newspec.cells.label},u.to));
  newspec.connections(from,to).label = connname;
end
cfg.focusconn = find(cellfun(@(x)isequal(connname,x),{newspec.connections.label}));
mechlist = newspec.connections(cfg.focusconn).mechanisms;
str = get(src,'string');
newmechlist = strtrim(splitstr(str,','));
mechadded = setdiff(newmechlist,mechlist);
mechremoved = setdiff(mechlist,newmechlist);
if ~isempty(mechadded) || ~isempty(mechremoved)
  newspec = removemech(newspec,mechremoved,'connections',cfg.focusconn);
  newspec = addmech(newspec,mechadded,'connections',cfg.focusconn);
  if isempty(newspec.connections(cfg.focusconn).mechanisms)
    newspec.connections(cfg.focusconn).label = [];
  end
  updatemodel(newspec);
end

function seladj(src,evnt)

function selkernel(src,evnt)

function ZoomFunction(src,evnt)
