function modeler(varargin)
clear global cfg H CURRSPEC
global cfg H CURRSPEC LASTSPEC BIOSIMROOT %currspec lastspec
prepare_spec;
updatemodel(CURRSPEC);

%% set up GUI

% main figure
sz = get(0,'ScreenSize'); 
fig = figure('position',[.005*sz(3) .07*sz(4) .938*sz(3) .83*sz(4)],'color','w','name','Biosim Modeler','NumberTitle','off','WindowScrollWheelFcn',@ZoomFunction,'CloseRequestFcn','delete(gcf); clear global H'); % [320 240 920 560]

% global controls (i.e., always present in main figure in all views)
  uicontrol('parent',fig,'style','text','string','biosim','fontsize',24,'units','normalized','position',[.1 .9 .2 .07],'backgroundcolor','w');
% tabs:
  bnet =uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[0 .85 .1 .04],'string','net','backgroundcolor',[.7 .7 .7],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pnet''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bcell=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[.2 .85 .1 .04],'string','cell','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pcell''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bmech=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[.1 .85 .1 .04],'string','mech','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pmech''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bview=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[.3 .85 .1 .04],'string','model','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pview''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
% model controls:
  bsave=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[0 .98 .05 .025],'string','save','backgroundcolor',[.8 .8 .8],'callback',@Save_Spec);
  bload=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.0525 .98 .02 .025],'string','o','backgroundcolor',[.8 .8 .8],'callback',@Load_Spec);
  bclear=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.098 .98 .02 .025],'string','x','backgroundcolor',[.8 .8 .8],'callback',{@SelectCells,1});
  %bclear=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.0+.075 .98 .02 .025],'string','x','backgroundcolor',[.8 .8 .8],'callback',[]);
  %bprint=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.0+0 .955 .05 .02],'string','print','backgroundcolor',[.8 .8 .8],'callback',@printmodel);
% change controls:
  %bapply=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.2+.1 .98 .05 .025],'string','refresh','backgroundcolor',[.8 .8 .8],'callback',{@refresh,0});
  bapply=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.34 .98 .05 .025],'string','refresh','backgroundcolor',[.8 .8 .8],'callback',{@refresh,0});
  bundo=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.075 .98 .02 .025],'string','<-','backgroundcolor',[.8 .8 .8],'callback',@undo);
  %bundo=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.2+.1525 .98 .02 .025],'string','<-','backgroundcolor',[.8 .8 .8],'callback',@undo);
  %brevert=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.2+.175 .98 .02 .025],'string','<=','backgroundcolor',[.8 .8 .8],'callback',{@refresh,2});

% left panels for cell, network, and mechanism controls
pnet =uipanel('parent',fig,'title','network builder','visible','on','tag','ptoggle','userdata','pnet','units','normalized','position',[0 0 .4 .85]);
pcell=uipanel('parent',fig,'title','cell builder','visible','off','tag','ptoggle','userdata','pcell','units','normalized','position',[0 0 .4 .85]);
pmech=uipanel('parent',fig,'title','mechanism editor','visible','off','tag','ptoggle','userdata','pmech','units','normalized','position',[0 0 .4 .85]);
pview=uipanel('parent',fig,'title','','visible','off','tag','ptoggle','userdata','pview','units','normalized','position',[0 0 .4 .85]);

% left panel: model view
txt_model = uicontrol('parent',pview,'style','edit','units','normalized','tag','modeltext',...
  'position',[0 0 1 1],'string',cfg.modeltext,'FontName','Monospaced','FontSize',9,'HorizontalAlignment','Left','Max',100,'BackgroundColor',[.9 .9 .9]);
  % enable horizontal scrolling
  jEdit = findjobj(txt_model);
  jEditbox = jEdit.getViewport().getComponent(0);
  jEditbox.setWrapping(false);                % turn off word-wrapping
  jEditbox.setEditable(false);                % non-editable
  set(jEdit,'HorizontalScrollBarPolicy',30);  % HORIZONTAL_SCROLLBAR_AS_NEEDED
  % maintain horizontal scrollbar policy which reverts back on component resize 
  hjEdit = handle(jEdit,'CallbackProperties');
  set(hjEdit, 'ComponentResizedCallback','set(gcbo,''HorizontalScrollBarPolicy'',30)')

% left panel: network builder %GUI_netpanel;
p_net_select  = uipanel('parent',pnet,'Position',[.02 .67 .96 .3],'BorderWidth',.2,'BorderType','line'); % cell morphology
p_net_connect = uipanel('parent',pnet,'Position',[.02 .44 .96 .34/1.5],'BorderWidth',.2,'BorderType','line','title','Population connections         [targets]'); % cell specification
p_net_kernel  = uipanel('parent',pnet,'Position',[.02 .01 .96 .42],'BorderWidth',.2,'BorderType','line','title','Cell connections'); % cell specification
% compartment controls
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[0 .88 .25 .06],'string','Compartments','ListboxTop',0,'HorizontalAlignment','left');
if ~isempty(net.cells) && ischar(net.cells(1).label)
  l1={net.cells.label}; 
  l2={net.cells.parent};
  l=cellfun(@(x,y)[x '.' y],l2,l1,'uni',0);
  i=1:length(l);
else
  l={}; i=[];
end
lst_comps = uicontrol('parent',p_net_select,'units','normalized',...
  'style','listbox','position',[0 .07 .19 .8],'value',i,'string',l,...
  'backgroundcolor','w','Max',5,'Min',0,'Callback',@SelectCells);
% headers for cell info
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.24 .88 .09 .06],'string','cell','ListboxTop',0,'HorizontalAlignment','left');
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.33 .88 .1 .06],'string','comp','ListboxTop',0,'HorizontalAlignment','left');
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.44 .88 .06 .06],'string','n','ListboxTop',0,'HorizontalAlignment','left');
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.51 .88 .15 .06],'string','mechanisms','ListboxTop',0,'HorizontalAlignment','left');

% left panel: cell builder %GUI_cellpanel;
p_cell_morph = uipanel('parent',pcell,'Position',[0 .6 .35 .35],'BorderWidth',.2,'BorderType','line','visible','off'); % cell morphology
p_cell_parms = uipanel('parent',pcell,'Position',[0 .4 1 .6],'BorderWidth',.2,'BorderType','line','title','Parameters');
p_cell_spec = uipanel('parent',pcell,'Position',[0 0 1 .4],'BorderWidth',.2,'BorderType','line','title','','fontangle','italic'); % cell specification
edit_comp_dynamics = uicontrol('parent',pcell,'style','edit','string','',...
    'units','normalized','position',[.1 .1 .8 .05],'BackgroundColor','w','HorizontalAlignment','left','tooltipstring','Mechanism functions will be substituted here.');
  
% left panel: mechanism editor %GUI_mechpanel;
% compartment label
if ~isempty(CURRSPEC.(cfg.focustype))
  tmp=CURRSPEC.(cfg.focustype)(cfg.focus);
  str1=tmp.mechanisms;
  str2=mech_spec2str(tmp.mechs(cfg.focusmech));
  cl=tmp.label;
  u.focustype=cfg.focustype; 
  u.focus=cfg.focus; 
  u.mechlabel=tmp.mechanisms{cfg.focusmech};  
else
  str1='';
  str2='';
  cl='';
  u=[];
end
txt_comp = uicontrol('style','text','string',cl,...
  'units','normalized','position',[0 .95 .1 .05],'parent',pmech);
% dropdown list of mechanisms for this compartment
lst_mechs = uicontrol('style','popupmenu','value',min(1,length(str1)),'string',str1,...
  'units','normalized','position',[.11 .95 .2 .05],'parent',pmech,'callback',@Display_Mech_Info);
% button to apply changes to mech text
uicontrol('parent',pmech,'style','pushbutton','units','normalized','position',[.35 .97 .1 .03],...
  'string','apply','callback',@UpdateMech);
% button to display list of mechs in DB
uicontrol('parent',pmech,'style','pushbutton','units','normalized','position',[.95 .97 .05 .03],'string','DB','callback','global allmechs; msgbox({allmechs.label},''available'');');%msgbox(get_mechlist,''available'')');%'get_mechlist');
% edit box with mech info
txt_mech = uicontrol('parent',pmech,'style','edit','units','normalized','BackgroundColor','w',... % [.9 .9 .9]
  'position',[0 .6 1 .35],'string',str2,'userdata',u,'FontName','courier','FontSize',10,'HorizontalAlignment','Left','Max',100);
% mech plots associated w/ this compartment
p_static_plots = uipanel('parent',pmech,'Position',[0 0 1 .6],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','');
ax_static_plot = subplot('position',[.04 .45 .9 .5],'parent',p_static_plots,'linewidth',3,'color','w','fontsize',6); box on; 
lst_static_funcs = uicontrol('units','normalized','position',[.04 .02 .9 .35],'parent',p_static_plots,...
  'style','listbox','value',1:5,'string',{},'Max',50,'Callback',@DrawAuxFunctions);
uicontrol('Style','edit', 'Units','normalized','Position',[0.8 0.45 0.2 0.05],...
          'String',sprintf('[%g,%g]',min(cfg.V),max(cfg.V)),'Callback',[],'parent',p_static_plots);
if ~isempty(CURRSPEC.cells)
  maxlhs=20; maxlen=150; % limit how much is shown in the listbox
  funcs = CURRSPEC.cells(cfg.focuscomp).functions;
  len = min(maxlhs,max(cellfun(@length,funcs(:,1))));
  str = {};
  for i=1:size(funcs,1)
    str{i} = sprintf(['%-' num2str(len) 's  = %s'],funcs{i,1},strrep(funcs{i,2},' ',''));
    if length(str{i})>maxlen, str{i}=str{i}(1:maxlen); end
  end
  val=get(lst_static_funcs,'value');
  set(lst_static_funcs,'string',str);
  set(lst_static_funcs,'value',val(val<=length(str)));
end

% right panel: simulation plots and controls %GUI_simpanel;
psims=uipanel('parent',fig,'title','','visible','on','units','normalized','position',[.4 0 .6 1]);

% menu %GUI_menu;
set(fig,'MenuBar','none');
file_m = uimenu(fig,'Label','File');
uimenu(file_m,'Label','Load sim_data','Callback',@Load_Data);
uimenu(file_m,'Label','Load spec','Callback',@Load_Spec);
uimenu(file_m,'Label','Save spec','Callback',@Save_Spec);
uimenu(file_m,'Label','Grab spec','Callback','global CURRSPEC; assignin(''base'',''spec'',CURRSPEC);');
uimenu(file_m,'Label','Update spec from base','Callback',{@refresh,1});
uimenu(file_m,'Label','Exit','Callback','global CURRSPEC H cfg; close(H.fig); clear CURRSPEC H cfg;');
plot_m = uimenu(fig,'Label','Plot');
uimenu(plot_m,'Label','plotv','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), plotv(evalin(''base'',''sim_data''),CURRSPEC); else disp(''load data to plot''); end');
uimenu(plot_m,'Label','plotpow','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), plotpow(evalin(''base'',''sim_data''),CURRSPEC,''spectrogram_flag'',0); else disp(''load data to plot''); end');
uimenu(plot_m,'Label','plotspk','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), plotspk(evalin(''base'',''sim_data''),CURRSPEC,''window_size'',30/1000,''dW'',5/1000); else disp(''load data to plot''); end');
uimenu(plot_m,'Label','visualizer','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), visualizer(evalin(''base'',''sim_data'')); else disp(''load data to plot''); end');

% collect object handles
% figures
H.fig   = fig;
% panels
H.pcell = pcell;
H.pnet  = pnet;
H.pmech = pmech;
H.pview = pview;
H.psims = psims;
H.p_net_select  = p_net_select;
H.p_net_connect = p_net_connect;
H.p_net_kernel  = p_net_kernel;
H.p_cell_morph  = p_cell_morph;
H.p_cell_spec   = p_cell_spec;
H.p_cell_parms  = p_cell_parms;
H.p_static_plots= p_static_plots;
% buttons
H.bcell = bcell; 
H.bnet  = bnet;
H.bmech = bmech;
H.bview = bview; % @callback (name callback functions here...)
H.bsave = bsave; % @callback
H.bload = bload; % @callback
%H.bclear= bclear;
H.bapply= bapply;
H.bundo = bundo;
%H.brevert=brevert;
% list boxes
H.lst_comps = lst_comps; % compartment listbox in cell view
H.lst_mechs = lst_mechs; % mechanism listbox in mech view
H.txt_mech = txt_mech; % edit box with mech info
H.txt_comp = txt_comp; % text label in mech view
H.txt_model = txt_model; % text control for readable model text
% menu items
H.file_m = file_m;
H.plot_m = plot_m;

H.ax_static_plot = ax_static_plot;
H.lst_static_funcs = lst_static_funcs;
H.edit_comp_dynamics = edit_comp_dynamics; % cell dynamics for focuscomp

% populate controls
if ~isempty(CURRSPEC.cells)
  SelectCells;
  DrawAuxView;
  DrawSimPlots;
  DrawUserParams([],[],[],0);
end

%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    if isfield(o.spec,'entities') 
      if ~isfield(o.spec,'cells')
        o.spec.cells=o.spec.entities;
      end
      o.spec=rmfield(o.spec,'entities');
      if ~isfield(o.spec.cells,'parent')
        for i=1:length(o.spec.cells)
          o.spec.cells(i).parent=o.spec.cells(i).label;
        end
      end
    end
    global CURRSPEC
    if isfield(CURRSPEC,'cells') && ~isempty(CURRSPEC.cells)
      n=length(o.spec.cells);
      [addflds,I]=setdiff(fieldnames(CURRSPEC.cells),fieldnames(o.spec.cells));
      [jnk,I]=sort(I);
      addflds=addflds(I);
      for i=1:length(addflds)
        o.spec.cells(1).(addflds{i})=[];
      end
      CURRSPEC.cells(end+1:end+n) = o.spec.cells;
      for i=1:n
        CURRSPEC.connections(end+i,end+1:end+n) = o.spec.connections(i,:);
      end
      if isfield(o.spec,'files')
        CURRSPEC.files = unique({CURRSPEC.files{:},o.spec.files{:}});
      end
    else
      CURRSPEC = o.spec;
    end
    refresh(src,evnt);
  else
    fprintf('select file does not contain spec structure. no specification loaded\n');
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Save_Spec(src,evnt)
[filename,pathname] = uiputfile({'*.mat;'},'Save as','model-specification.mat');
if isequal(filename,0) || isequal(pathname,0)
  return;
end
outfile = fullfile(pathname,filename);
[fpath,fname,fext] = fileparts(outfile);
global CURRSPEC
spec=CURRSPEC;
fprintf('Saving model specification: %s\n',outfile);
save(outfile,'spec');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SelectCells(src,evnt,null)
global H CURRSPEC cfg
if nargin<3, null=0; end
if null
  CURRSPEC=[];
  set(H.lst_comps,'string',{},'value',[]);
  try
    set(H.edit_comp_parent,'visible','off');
    set(H.edit_comp_label,'visible','off');
    set(H.edit_comp_N,'visible','off');
    set(H.edit_comp_mechs,'visible','off');
    set(H.btn_comp_copy,'visible','off');
    set(H.btn_comp_delete,'visible','off');
    set(H.btn_comp_edit,'visible','off');    
    set(H.p_comp_mechs,'visible','off');    
  end
  return; 
end
v=get(H.lst_comps,'value'); 
l=get(H.lst_comps,'string');
if ischar(CURRSPEC.cells(1).label)
  if isfield(CURRSPEC.cells,'parent') && ~isempty(CURRSPEC.cells(1).parent)
    l1={CURRSPEC.cells.label}; 
    l2={CURRSPEC.cells.parent};
    l=cellfun(@(x,y)[x '.' y],l2,l1,'uni',0);
  else
    l={CURRSPEC.cells.label}; 
  end
else
  l={};
end 
set(H.lst_comps,'string',l,'value',v(v<=length(CURRSPEC.cells)));
if isfield(H,'p_comp_mechs')
  set(H.p_comp_mechs,'visible','off');
end
DrawCellInfo(CURRSPEC); % Selection panel
DrawNetGrid(CURRSPEC);  % Connection panel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawCellInfo(net)
global H
c=1.5; dy=-.07*c; 
sel = get(H.lst_comps,'value');
if isfield(net.cells,'parent')
  p={net.cells(sel).parent};
else
  [p{1:length(sel)}]=deal('');
end
l={net.cells(sel).label}; 
N=[net.cells(sel).multiplicity];
mechs={net.cells(sel).mechanisms};
for i=1:length(sel)
  m=mechs{i};
  str=m{1}; for j=2:length(m), str=[str ', ' m{j}]; end
  if ~isfield(H,'edit_comp_label') || length(H.edit_comp_label)<length(sel) || ~ishandle(H.edit_comp_label(i))
    H.edit_comp_parent(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.24 .8+dy*(i-1) .09 .08],'backgroundcolor','w','string',p{i},...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'parent'},'ButtonDownFcn',{@DrawUserParams,sel(i)});    
    H.edit_comp_label(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.33 .8+dy*(i-1) .1 .08],'backgroundcolor','w','string',l{i},...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'label'},'ButtonDownFcn',{@DrawUserParams,sel(i)});
    H.edit_comp_N(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.44 .8+dy*(i-1) .06 .08],'backgroundcolor','w','string',N(i),...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'multiplicity'});
    H.btn_comp_delete(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','-','callback',{@DeleteCell,l{i}},...
      'position',[.205 .8+dy*(i-1) .03 .08]);%,'BackgroundColor','white');    
    H.btn_comp_copy(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','+','callback',{@CopyCell,l{i}},...
      'position',[.93 .8+dy*(i-1) .03 .08]);%,'BackgroundColor','white');    
    H.btn_comp_edit(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','...','callback',{@ShowClickMechList,i,'cells'},...%{@OpenCellModeler,l{i}},...
      'position',[.965 .8+dy*(i-1) .03 .08]);%,'BackgroundColor','white');            
    H.edit_comp_mechs(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.51 .8+dy*(i-1) .42 .08],'backgroundcolor','w','string',str,...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'mechanisms'},...
      'ButtonDownFcn',{@Display_Mech_Info,l{i},{},'cells'});
    H.p_comp_mechs(i) = uipanel('parent',H.p_net_select,'units','normalized',...
      'position',[.51 .8+dy*(i-1) .42 .08],'visible','off');    
  else
    % update properties
    set(H.edit_comp_parent(i),'string',p{i},'visible','on','Callback',{@UpdateCells,l{i},'parent'});
    set(H.edit_comp_label(i),'string',l{i},'visible','on','Callback',{@UpdateCells,l{i},'label'});
    set(H.edit_comp_N(i),'string',N(i),'visible','on','Callback',{@UpdateCells,l{i},'multiplicity'});
    set(H.edit_comp_mechs(i),'string',str,'visible','on','Callback',{@UpdateCells,l{i},'mechanisms'},'ButtonDownFcn',{@Display_Mech_Info,l{i},{},'cells'});
    set(H.btn_comp_copy(i),'callback',{@CopyCell,l{i}},'visible','on');
    set(H.btn_comp_delete(i),'callback',{@DeleteCell,l{i}},'visible','on');
    set(H.btn_comp_edit(i),'callback',{@ShowClickMechList,i,'cells'},'visible','on');
    set(H.p_comp_mechs(i),'visible','off');    
  end
  if length(H.edit_comp_label)>i
    set(H.edit_comp_parent(i+1:end),'visible','off');
    set(H.edit_comp_label(i+1:end),'visible','off');
    set(H.edit_comp_N(i+1:end),'visible','off');
    set(H.edit_comp_mechs(i+1:end),'visible','off');
    set(H.btn_comp_copy(i+1:end),'visible','off');
    set(H.btn_comp_delete(i+1:end),'visible','off');
    set(H.btn_comp_edit(i+1:end),'visible','off');    
    set(H.p_comp_mechs(i),'visible','off');    
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawNetGrid(net)
global H
dx=.15; x=.13; c=1.5; dy=-.07*c; 
sel = get(H.lst_comps,'value');
net.cells = net.cells(sel);
net.connections = net.connections(sel,:);
net.connections = net.connections(:,sel);
l={net.cells.label}; 
for i=1:length(sel)
  for j=1:length(sel)
    m = net.connections(j,i).mechanisms;
    pos = [x+dx*(i-1) .8+dy*(j-1) .9*dx .08*c];
    u.from=l{j}; u.to=l{i};
    connname = net.connections(j,i).label;
    if ~isempty(m)
      if ischar(m), m={m}; end
      str=m{1}; for k=2:length(m), str=[str ', ' m{k}]; end
      ml=m{1};
    else
      str = '';
      ml={};
    end
    if ~isfield(H,'txt_to') || i>length(H.txt_from) || j>length(H.txt_to) || ~ishandle(H.edit_conn_mechs(i,j)) || H.edit_conn_mechs(i,j)==0
      if i==1 % to
        this=zeros(max(sel),1);
        this(sel)=j;
        H.txt_to(j) = uicontrol('parent',H.p_net_connect,'units','normalized',...
          'style','text','position',[x+dx*(j-1) .91 .11 .08*c],'string',l{j},...
          'callback',{@ShowClickMechList,this,'connections'});
      end
      if j==1 % from
        this=ones(1,max(sel));
        this(sel)=i;
        H.txt_from(i) = uicontrol('parent',H.p_net_connect,'units','normalized',...
          'style','text','position',[.01 .8+dy*(i-1) .11 .08*c],'string',l{i},...
          'callback',{@ShowClickMechList,this,'connections'});
      end
      H.edit_conn_mechs(i,j) = uicontrol('parent',H.p_net_connect,'units','normalized',...
        'style','edit','position',pos,'backgroundcolor','w',...
        'string',str,'HorizontalAlignment','left','UserData',u);
      set(H.edit_conn_mechs(i,j),'Callback',{@UpdateNet,connname},...
        'ButtonDownFcn',{@Display_Mech_Info,connname,ml,'connections'});
    H.p_conn_mechs(i,j) = uipanel('parent',H.p_net_connect,'units','normalized',...
      'position',pos,'visible','off');
    else
      set(H.txt_to(i),'string',l{i},'visible','on');
      set(H.txt_from(i),'string',l{i},'visible','on');
      set(H.p_conn_mechs(i,j),'visible','off');
      set(H.edit_conn_mechs(i,j),'string',str,'UserData',u,'Callback',{@UpdateNet,connname},...
        'ButtonDownFcn',{@Display_Mech_Info,connname,ml,'connections'},'visible','on');
    end
  end
end
if isfield(H,'txt_to') && length(H.txt_to)>length(sel)
  set(H.txt_to(i+1:end),'visible','off');
  set(H.txt_from(i+1:end),'visible','off');
  set(H.edit_conn_mechs(i+1:end,:),'visible','off');
  set(H.edit_conn_mechs(:,i+1:end),'visible','off');
  set(H.p_conn_mechs(i+1:end,:),'visible','off');
  set(H.p_conn_mechs(:,i+1:end),'visible','off');
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ShowClickMechList(src,evnt,this,type)
global H CURRSPEC cfg
sel=get(H.lst_comps,'value');
if strcmp(type,'cells')
  % cell mech selected in comp mechlist
  h=H.edit_comp_mechs(this);
  p=H.p_comp_mechs(this);
  cfg.focuscomp=this;
else
  % connection mech selected in NetGrid
  if numel(this)>1
    % clicked to or from label; toggle panel...
    % ...
    return; % not implemented yet
  else
    % clicked connections
    % ...
    cfg.focusconn=this;
  end
end
cfg.focustype=type;
cfg.focus=sel(this);
if strcmp(get(h,'visible'),'off')
  set(p,'visible','off');
  set(h,'visible','on');
  return;
end
set(p,'visible','on');
set(h,'visible','off');
children=get(p,'Children');
if ~isempty(children)
  delete(children);
end
l=CURRSPEC.(type)(cfg.focus).label;
m=CURRSPEC.(type)(cfg.focus).mechanisms;
n=length(m);
if length(n)<1, return; end
w=1/n;
x=(0:n-1)*w;
for i=1:n
  uicontrol('parent',p,'style','text','string',m{i},'units','normalized','position',[x(i) 0 w 1],...
    'tooltip',m{i},'ButtonDownFcn',{@Display_Mech_Info,l,m{i}});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawAuxView
global H CURRSPEC
if isempty(CURRSPEC.cells), return; end
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
mech_u=[];
for i=1:length(sel)
  m=CURRSPEC.connections(sel(i)).mechanisms;
  if ~iscell(m), m={m}; end
  [from,to]=ind2sub(size(CURRSPEC.connections),sel(i));
  for j=1:length(m)
    mech_u(end+1).conn_i=sel(i);
    mech_u(end).mech_i=j;
    mech_u(end).from=from;
    mech_u(end).to=to;
    mech_u(end).label=[l{from} '_' l{to}];
    lst{end+1} = [ll{sel(i)} '.' m{j}];
  end
  if i==1
    prefix = [l{from} '_' l{to}];
    % get expression for the last auxvar
    a=CURRSPEC.connections(from,to).mechs(1).auxvars(end,:);
    auxeqn=[a{1} ' = ' a{2}];
    % get list of params for the first connection
    parmlist=fieldnames(CURRSPEC.connections(from,to).mechs(1).params);
    % get value of first parameter
    parmval=CURRSPEC.connections(from,to).mechs(1).params.(parmlist{1});    
  end
end
% get list of auxvars for the first connection
sel=strmatch([prefix '_'],CURRSPEC.model.auxvars(:,1));
auxlist=CURRSPEC.model.auxvars(sel,1);
if isempty(auxlist), return; end
% get auxvar matrix 
auxmat=CURRSPEC.model.eval.(auxlist{end});
midpt=round(size(auxmat,2)/2); 
lims=[min(auxmat(:)) max(auxmat(:))];%[.9*min(auxmat(:)) 1.1*max(auxmat(:))];

% dropdown box to select the connection to examine
H.pop_sel_conn = uicontrol('style','popupmenu','string',lst,'parent',H.p_net_kernel,...
  'units','normalized','position',[.1 .92 .3 .08],'value',1,'userdata',mech_u,...
  'Callback',@UpdateAuxView);
% listbox to select which auxvar to plot
H.lst_auxvars = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','listbox','position',[.1 .1 .3 .3],'string',auxlist,'value',length(auxlist),...
  'Callback',@UpdateAuxView); % 'backgroundcolor','w',
% edit box to modify defining expressions
H.edit_auxvar_eqn = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','edit','position',[.1 .03 .3 .08],'backgroundcolor','w','string',auxeqn,...
  'HorizontalAlignment','left','fontsize',8);%,'Callback',{@UpdateParams,'change','cells'});
% listbox to select a parameter to modify
H.lst_auxparams = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','listbox','position',[.44 .1 .3 .3],'string',parmlist,'value',1,...
  'Callback',@UpdateAuxView); % 'backgroundcolor','w',
% edit box to modify parameter
H.edit_auxparams = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','text','position',[.44 .03 .3 .08],'backgroundcolor',[.9 .9 .9],'string',num2str(parmval),...
  'HorizontalAlignment','left');%,'Callback',{@UpdateParams,'change','cells'});

% plots
H.ax_conn_img = subplot('position',[.1 .5 .3 .4],'parent',H.p_net_kernel); 
H.img_connect = imagesc(auxmat); axis xy; vline(midpt,'k'); 
if lims(2)>lims(1), caxis(lims); end
colorbar

H.ax_conn_line = subplot('position',[.44 .5 .3 .4],'parent',H.p_net_kernel); 
H.line_connect = line('ydata',1:size(auxmat,1),'xdata',auxmat(:,midpt),'color','k','LineStyle','-','erase','background');
if lims(2)>lims(1), xlim(lims); end
% create button group for predefined adjacency matrices
H.rad_adj = uibuttongroup('visible','off','units','normalized','Position',[.77 .5 .2 .3],'parent',H.p_net_kernel);
H.rad_adj_1 = uicontrol('Style','radiobutton','String','1-to-1','units','normalized',...
    'pos',[.05 .75 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
H.rad_adj_2 = uicontrol('Style','radiobutton','String','all-to-all','units','normalized',...
    'pos',[.05 .45 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
H.rad_adj_3 = uicontrol('Style','radiobutton','String','random','units','normalized',...
    'pos',[.05 .15 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
set(H.rad_adj,'SelectionChangeFcn',@seladj);
set(H.rad_adj,'SelectedObject',[]);  % No selection
set(H.rad_adj,'Visible','on');

function UpdateAuxView(src,evnt)
global H CURRSPEC
str=get(H.pop_sel_conn,'string');
val=get(H.pop_sel_conn,'value');
mech_u=get(H.pop_sel_conn,'userdata');
mech_u=mech_u(val);
from=mech_u.from;
to=mech_u.to;
label=mech_u.label;
mech_i=mech_u.mech_i;
% update aux list
if strcmp(get(src,'style'),'popupmenu')
  sel=strmatch([strrep(label,'-','_') '_'],CURRSPEC.model.auxvars(:,1));
  auxlist=CURRSPEC.model.auxvars(sel,1);
  auxind=length(auxlist);
  set(H.lst_auxvars,'string',auxlist,'value',auxind);
else
  auxlist=get(H.lst_auxvars,'string');
  auxind=get(H.lst_auxvars,'value');
end
if isempty(auxlist), return; end
% get expression for aux var
a=CURRSPEC.connections(from,to).mechs(mech_i).auxvars(auxind,:);
auxeqn=[a{1} ' = ' a{2}];
% get list of params
%parmlist=fieldnames(CURRSPEC.connections(from,to).mechs(mech_i).params);
parmlist=get(H.lst_auxparams,'string');
parmsel=get(H.lst_auxparams,'value');
% get value of first parameter
try
  parmval=CURRSPEC.connections(from,to).mechs(mech_i).params.(parmlist{parmsel});
catch
  parmval=[];
end
% get auxvar matrix 
auxmat=CURRSPEC.model.eval.(auxlist{auxind});
midpt=round(size(auxmat,2)/2); 
lims=[min(auxmat(:)) max(auxmat(:))];
% update uicontrols
set(H.edit_auxvar_eqn,'string',auxeqn);
set(H.lst_auxparams,'string',parmlist);
set(H.edit_auxparams,'string',num2str(parmval));
set(H.img_connect,'cdata',auxmat);
set(H.line_connect,'ydata',1:size(auxmat,1),'xdata',auxmat(:,midpt));
if lims(2)>lims(1)
  set(H.ax_conn_img,'clim',lims);
  set(H.ax_conn_line,'xlim',lims);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulate(src,evnt,action)
global CURRSPEC H cfg
DrawSimPlots;
functions = CURRSPEC.model.functions;
auxvars = CURRSPEC.model.auxvars;
ode = CURRSPEC.model.ode;
IC = CURRSPEC.model.IC;
allvars=CURRSPEC.variables.labels;
allinds=CURRSPEC.variables.entity;
cfg.T = (0:cfg.buffer-1)*cfg.dt; % ms
fh_pow = @(x,y) PowerSpecTA(x,y,[10 80],min(8000,cfg.buffer/2),'Normalized',[]);

% get list of indices of vars to plot
list=get(H.lst_comps,'string');
sel=get(H.lst_comps,'value');
if length(sel)>3, sel=sel(1:3); end
plotvars=cell(size(sel));
show=list(sel);
for k=1:length(show)
  this=sel(k);
  list=get(H.lst_vars(k),'string');
  vind=get(H.lst_vars(k),'value');
  var = list{vind};% CURRSPEC.cells(this).ode_labels{1};
  plotvars{k}=find(cellfun(@(x)isequal(x,var),allvars));
end
plotfuncs=cell(size(sel));

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
X=IC; 
t=0; 
cnt=0; 
var_flag=ones(size(sel));
tvec=zeros(1,cfg.buffer);
cfg.record=zeros(length(IC),cfg.buffer);
while cfg.quitflag<0 && (length(IC)==length(CURRSPEC.model.IC))
  if ~isequal(ode,CURRSPEC.model.ode)
    ode=CURRSPEC.model.ode;
    allvars=CURRSPEC.variables.labels;
    allinds=CURRSPEC.variables.entity;    
    F=eval(ode);
    functions = CURRSPEC.model.functions;
    auxvars = CURRSPEC.model.auxvars;
    % update auxiliary variables (ie., adjacency matrices)
    for k = 1:size(auxvars,1)
      %try  % added to catch mask=mask-diag(diag(mask)) when mask is not square
        eval(sprintf('%s = %s;',auxvars{k,1},auxvars{k,2}) );
      %end
    end
    % update anonymous functions
    for k = 1:size(functions,1)
      eval(sprintf('%s = %s;',functions{k,1},functions{k,2}) );
    end      
  end
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
    tvec(cnt)=t;
  else
    cfg.record = [cfg.record(:,2:end) X];
    tvec(2:end)=t;
    %if mod(cnt,100)==0
      %set(H.ax_state_plot,'xtick',get(H.ax_state_plot(1),'xtick')+100*cfg.dt);
    %end
  end
  if mod(cnt,round(cfg.buffer/50))==0
    % get output to plot
    list=get(H.lst_comps,'string');
    sel=get(H.lst_comps,'value');
    if length(sel)>3, sel=sel(1:3); end
    show=list(sel);
    if cfg.changeflag>0
      % get vars to plot from listboxes
      for k=1:length(show)
        this=sel(k);
        numcell=CURRSPEC.cells(this).multiplicity;
        list=get(H.lst_vars(k),'string');
        vind=get(H.lst_vars(k),'value');
        var = list{vind};% CURRSPEC.cells(this).ode_labels{1};
        if ismember(var,allvars) % plot state var
          var_flag(k)=1;
          plotvars{k}=find(cellfun(@(x)isequal(x,var),allvars));
        elseif ismember(var,functions(:,1)) % plot aux function
          var_flag(k)=0;
        end
        if var_flag(k)==1
          set(get(H.ax_state_plot(k),'title'),'string',sprintf('%s (n=%g/%g)',strrep(list{vind},'_','\_'),min(numcell,cfg.ncellshow),numcell));
        else
          set(get(H.ax_state_plot(k),'title'),'string',sprintf('%s (n=%g/%g): %s',strrep(list{vind},'_','\_'),min(numcell,cfg.ncellshow),numcell,strrep(functions{find(strcmp(var,functions(:,1))),2},'_','\_')));
        end
      end    
      cfg.changeflag=-1;
    end
    for k=1:length(show)
      this=sel(k);
      numcell=CURRSPEC.cells(this).multiplicity;
      if numcell > cfg.ncellshow
        tmp=randperm(numcell);
        inds = sort(tmp(1:cfg.ncellshow));
      else
        inds = 1:numcell;
      end
      %for j=1:length(inds)
        if var_flag(k)==1 % plot state var
          for j=1:length(inds)
            set(H.simdat_alltrace(k,j),'ydata',cfg.record(plotvars{k}(inds(j)),:));
          end
        elseif var_flag(k)==0 % plot aux function
          % -----------------------------------------------------------
          list=get(H.lst_vars(k),'string');
          vind=get(H.lst_vars(k),'value');
          var = list{vind};% CURRSPEC.cells(this).ode_labels{1};          
          thisfunc=find(strcmp(var,functions(:,1)));
          basename=functions{thisfunc,3};
          funceqn=functions{thisfunc,2};
          args=regexp(funceqn,'^@\([\w,]*\)','match');
          args=splitstr(args{1}(3:end-1),',');
          plotfargs=cell(1,length(args));
          % loop over args
          plotfunc_flag=1;
          for a=1:length(args)
            plotfargs{a}=zeros(numcell,cfg.buffer);
            arg=args{a}; 
            if isequal(arg,'t') % current time vector
              for b=1:numcell
                plotfargs{a}(b,:)=tvec;%t+(0:cfg.buffer-1)*cfg.dt;
              end
            elseif ismember(arg,allvars) % known state var
              tmp=find(cellfun(@(x)isequal(x,arg),allvars));
              plotfargs{a}(:,:)=cfg.record(tmp,:);
            else % unknown var. need to map var labels to get state var indices
              try
                spc=CURRSPEC.cells(this);
                ii=find(~cellfun(@isempty,{spc.mechs.functions}));
                jj=arrayfun(@(x)any(ismember(basename,x.functions(:,1))),spc.mechs(ii));
                if any(jj)
                  m=ii(jj);
                  args2=regexp(spc.mechs(m).substitute(:,2),[basename '\([\w,]*\)'],'match');
                else
                  ii=find(~cellfun(@isempty,{spc.connection_mechs.functions}));
                  jj=arrayfun(@(x)any(ismember(basename,x.functions(:,1))),spc.connection_mechs(ii));
                  m=ii(jj);
                  args2=regexp(spc.connection_mechs(m).substitute(:,2),[basename '\([\w,]*\)'],'match');
                end
                args2=args2{1};
                if iscell(args2), args2=args2{1}; end
                args2=splitstr(args2((length(basename)+2):end-1),',');
                for b=a:length(args2)
                  origvar=args2{a};
                  thisarg=spc.var_list{strmatch(origvar,spc.orig_var_list)};
                  if ismember(thisarg,allvars)
                    tmp=find(cellfun(@(x)isequal(x,thisarg),allvars));
                    plotfargs{a}(:,:)=cfg.record(tmp,:);
                  end
                end
                plotfunc_flag=1;
              catch
                plotfunc_flag=0;
              end
            end
          end
          for j=1:length(inds)
            if plotfunc_flag
              plotfuncs=eval(sprintf('%s(plotfargs{:})',var)); %feval(var,plotfargs{:});      
              set(H.simdat_alltrace(k,j),'ydata',plotfuncs(inds(j),:));
            else
              set(H.simdat_alltrace(k,j),'ydata',zeros(1,cfg.buffer));
            end
          end
          % -----------------------------------------------------------
        end
      %end
  %     if 0 
  %       lfp=mean(cfg.record(plotvars{k},:),1);
  %       set(H.simdat_LFP(k),'ydata',lfp,'linewidth',4,'color','k','linestyle','-');
  %     end
  %     if 0 && mod(cnt,cfg.buffer)==0 && cnt>1
  %       % update power spectrum
  %       try
  %         res = feval(fh_pow,cfg.T,lfp);
  %         set(H.simdat_LFP_power(k),'xdata',res.f,'ydata',log10(res.Pxx));
  %       catch
  %         fprintf('power spectrum calc failed.\n');
  %       end
  %     end
    end
    drawnow
  end
end
cfg.quitflag=-1;
p=findobj('tag','pause');
set(p,'string','pause'); 
set(findobj('tag','start'),'string','restart');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
allvars=CURRSPEC.variables.labels;
allinds=CURRSPEC.variables.entity;
numcells=[CURRSPEC.cells.multiplicity];
% data plots
if isfield(H,'ax_state_plot')
  delete(H.simdat_alltrace(ishandle(H.simdat_alltrace)));
  delete(H.simdat_LFP(ishandle(H.simdat_LFP)));
  delete(H.ax_state_plot(ishandle(H.ax_state_plot)));
  delete(H.simdat_LFP_power(ishandle(H.simdat_LFP_power)));
  delete(H.lst_vars(ishandle(H.lst_vars)));
  %delete(H.ax_state_power(ishandle(H.ax_state_power)));
  H=rmfield(H,{'simdat_alltrace','simdat_LFP','ax_state_plot','simdat_LFP_power'});%,'ax_state_power'});
end
for i=1:length(show) % i=1:ncomp
  % get var labels
  vars=unique(allvars(allinds==sel(i)));
  % add function labels
  vars=cat(2,vars,CURRSPEC.cells(sel(i)).functions(:,1)');
  vind=find(~cellfun(@isempty,regexp(vars,'_V$','once')));
  if isempty(vind), vind=1; end
  H.lst_vars(i) = uicontrol('parent',H.psims,'units','normalized','backgroundcolor','w',...
    'style','listbox','position',[.8 .7+(i-1)*dy .2 -.8*dy],'value',vind,'string',vars,'Callback','global cfg;cfg.changeflag=1;');  
  H.ax_state_plot(i) = subplot('position',[.03 .7+(i-1)*dy .76 -.8*dy],...
    'parent',H.psims,'linewidth',3,'color','w'); 
  for k=1:cfg.ncellshow
    H.simdat_alltrace(i,k)=line('color',cfg.colors(max(1,mod(k,length(cfg.colors)))),'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',cfg.T,'ydata',zeros(1,cfg.buffer),'zdata',[]);
  end
  H.simdat_LFP(i)=line('color',cfg.colors(max(1,mod(k,length(cfg.colors)))),'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',cfg.T,'ydata',zeros(1,cfg.buffer),'zdata',[],'linewidth',2);
  axis([cfg.T(1) cfg.T(end) -100 30]); 
  if i==length(show), xlabel('time (ms)'); end
  titlestr=sprintf('%s (n=%g/%g)',strrep(vars{vind},'_','\_'),min(numcells(sel(i)),cfg.ncellshow),numcells(sel(i)));
  title(titlestr);
  %title(strrep(vars{vind},'_','\_'));%[show{i} '.V']);  %ylabel([show{i} '.V']);  
  %H.ax_state_power(i) = subplot('position',[.7 .7+(i-1)*dy .25 -.8*dy],'parent',H.psims); 
  H.simdat_LFP_power(i)=line('color','k','LineStyle','-','erase','background','xdata',cfg.f,'ydata',zeros(size(cfg.f)),'zdata',[],'linewidth',2);
  %H.ax_state_img(i) = subplot('position',[.7 .7+(i-1)*dy .25 -.8*dy],'parent',H.psims); 
  %H.img_state = imagesc(zeros(10,40)); axis xy; axis square
end
% % slider control
if isempty(findobj('tag','speed'))
  uicontrol('Style','frame', 'Units','normalized', ...
            'Position',[0.1  0.05 0.41 0.05],'parent',H.psims);
  uicontrol('Style','text', 'Units','normalized',...
            'Position',[0.15  0.075 0.3 0.02],'parent',H.psims,'string','visualization speed');
  uicontrol('Style','slider', 'Units','normalized', ...
            'Position',[0.15  0.055 0.3 0.015],'parent',H.psims,...
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits(src,evnt,action)
global cfg H CURRSPEC
if ~isfield(cfg,'record'), return; end
allvars=CURRSPEC.variables.labels;
switch action
  case 'autoscale'
    for i=1:length(H.ax_state_plot)
%       list=get(H.lst_vars(i),'string');
%       ind=get(H.lst_vars(i),'value');
%       var = list{ind};
%       ind = find(cellfun(@(x)isequal(x,var),allvars));
%       if ~isempty(ind)
%         rec = cfg.record(ind,:);
%       else
        rec = get(H.simdat_alltrace(i,:),'ydata');
        rec = [rec{:}];
%       end
      ymin = min(rec(:));
      ymax = max(rec(:));
      if ymin~=ymax
        set(H.ax_state_plot(i),'ylim',[ymin ymax]);
      end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  if ischar(l) && ~isempty(regexp(l,['^' lold '-'],'once'))
    l=strrep(l,[lold '-'],[lnew '-']);
  end
  newspec.connections(n+1,i).label=l;
  l=newspec.connections(i,n+1).label;
  if ischar(l) && ~isempty(regexp(l,['-' lold '$'],'once'))
    l=strrep(l,['-' lold],['-' lnew]); 
  end  
  newspec.connections(i,n+1).label=l;
end
set(H.lst_comps,'string',{newspec.cells.label},'value',[get(H.lst_comps,'value') n+1]);
updatemodel(newspec);
SelectCells;
DrawAuxView;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DeleteCell(src,evnt,compname)
global CURRSPEC cfg
ind = strmatch(compname,{CURRSPEC.cells.label},'exact');
newspec = CURRSPEC;
newspec.cells(ind) = [];
newspec.connections(ind,:) = [];
newspec.connections(:,ind) = [];
cfg.focusconn = 1;
if cfg.focuscomp>=ind
  cfg.focuscomp = max(1,cfg.focuscomp-1);
end
%cfg.focusmech = 1;
updatemodel(newspec);
SelectCells(src,evnt);
DrawAuxView;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Display_Mech_Info(src,evnt,complabel,mechlabel,type)%mechpos,hmech,connname)
% purpose: display the mech model in a readable form
global H CURRSPEC cfg
if nargin<3 || isempty(complabel)
  complabel = get(H.txt_comp,'string');
end
if nargin<4 || isempty(mechlabel)
  s=get(src,'string');
  if isempty(s)
    return;
  elseif ischar(s)
    mechlabel=regexp(s,'[\w\d]+[^,]','once','match');
  elseif iscell(s)
    v=get(src,'value'); 
    mechlabel=s{v};
  end
end
if nargin<5 || isempty(type)
  type=cfg.focustype;
else
  cfg.focustype=type;
end
if isempty(mechlabel)
  return;
end
% switch to mechview:
set(findobj('tag','ptoggle'),'visible','off'); 
set(findobj('tag','tab'),'backgroundcolor',[1 1 1]); 
set(findobj('userdata','pmech'),'visible','on'); 
set(H.bmech,'backgroundcolor',[.7 .7 .7]);
% get mechanism info
if strcmp(cfg.focustype,'cells')
  comp_i = find(strcmp(complabel,{CURRSPEC.cells.label}),1,'first');
  cfg.focuscomp=comp_i;
  DrawAuxFunctions;
elseif strcmp(cfg.focustype,'connections')
  comp_i = find(cellfun(@(x)isequal(x,complabel),{CURRSPEC.connections.label}),1,'first');
  cfg.focusconn=comp_i;
end
mech_i = find(strcmp(mechlabel,CURRSPEC.(cfg.focustype)(comp_i).mechanisms));
cfg.focusmech=mech_i;
cfg.focus=comp_i;
mechs = CURRSPEC.(cfg.focustype)(cfg.focus).mechanisms;
mech = CURRSPEC.(cfg.focustype)(cfg.focus).mechs(mech_i);

% update compartment label
set(H.txt_comp,'string',CURRSPEC.(cfg.focustype)(cfg.focus).label);

% update dropdown list with mechanisms for this cell or cell-pair
set(H.lst_mechs,'string',mechs,'value',mech_i);

% show text of mech selected in dropdown list
u.focustype=cfg.focustype;
u.focus=cfg.focus;
u.mechlabel=mechlabel;
set(H.txt_mech,'string',mech_spec2str(mech),'userdata',u);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMech(src,evnt,htxt,connname,mechname)
global CURRSPEC
% purpose: write mech to new file
outfile = '';
txt = get(htxt,'string');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateMech(src,evnt)%,htxt,connname,mechname)
% purpose: apply user changes to the mech model
global CURRSPEC H
u=get(H.txt_mech,'userdata');
txt = get(H.txt_mech,'string');
newmech = parse_mech_spec(txt);
newmech.label = u.mechlabel;
spec=CURRSPEC;
this = spec.(u.focustype)(u.focus);
if ~iscell(this.mechanisms), this.mechanisms={this.mechanisms}; end
mech_i = find(strcmp(u.mechlabel,this.mechanisms));
if ~isempty(mech_i)
  this.mechs(mech_i) = newmech;
else
  this.mechs(end+1)=newmech;
  this.mechanisms{end+1}=u.mechlabel;
end
spec.(u.focustype)(u.focus) = this;
updatemodel(spec);
DrawAuxFunctions;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateCells(src,evnt,compname,field)
% apply user changes to the compartment mechanisms
global cfg CURRSPEC H 
this = strmatch(compname,{CURRSPEC.cells.label},'exact');
if numel(this)>1, this = this(1); end
newspec = CURRSPEC;
switch field
  case {'parent','label'}
    newspec.cells(this).(field) = get(src,'string');
    updatemodel(newspec);
    l1={newspec.cells.label}; 
    l2={newspec.cells.parent};
    l=cellfun(@(x,y)[x '.' y],l2,l1,'uni',0);    
    set(H.lst_comps,'string',l);
    SelectCells;    
  case 'multiplicity'
    newspec.cells(this).(field) = str2num(get(src,'string'));
    updatemodel(newspec);
  case 'mechanisms'
    mechlist = newspec.cells(this).mechanisms;
    str = get(src,'string');
    if isempty(str)
      newmechlist = {}
    else
      newmechlist = strtrim(splitstr(str,','));
    end
    mechadded = setdiff(newmechlist,mechlist);
    mechremoved = setdiff(mechlist,newmechlist);
    if ~isempty(mechadded) || ~isempty(mechremoved)
      newspec = removemech(newspec,mechremoved,'cells',this);
      newspec = addmech(newspec,mechadded,'cells',this);
      updatemodel(newspec);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateNet(src,evnt,connname)
% apply user changes to the connection mechanisms
global cfg CURRSPEC
newspec = CURRSPEC;
if isempty(connname) % first connection b/w these compartments...
  u=get(src,'UserData');
  connname = [u.from '-' u.to];
  source=find(strcmp({newspec.cells.label},u.from));
  target=find(strcmp({newspec.cells.label},u.to));
  newspec.connections(source,target).label = connname;
  set(src,'ButtonDownFcn',{@Display_Mech_Info,connname,{},'connections'});
end
cfg.focusconn = find(cellfun(@(x)isequal(connname,x),{newspec.connections.label}));
if ~isempty(cfg.focusconn)
  mechlist = newspec.connections(cfg.focusconn).mechanisms;
else
  cfg.focusconn = 1;
  mechlist = {};
end
str = get(src,'string');
if isempty(str)
  newmechlist = {};
else
  newmechlist = strtrim(splitstr(str,','));
end
mechadded = setdiff(newmechlist,mechlist);
mechremoved = setdiff(mechlist,newmechlist);
if ~isempty(mechadded) || ~isempty(mechremoved)
  newspec = removemech(newspec,mechremoved,'connections',cfg.focusconn);
  newspec = addmech(newspec,mechadded,'connections',cfg.focusconn);
  if isempty(newspec.connections(cfg.focusconn).mechanisms)
    newspec.connections(cfg.focusconn).label = [];
  end
  updatemodel(newspec);
  DrawAuxView;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CURRSPEC = addmech(CURRSPEC,mechadded,type,index)
% purpose: add a mechanism to the compartment model
global allmechs cfg
if isempty(mechadded), return; end
if ~iscell(mechadded), mechadded={mechadded}; end
for i=1:length(mechadded)
  newmech = mechadded{i};
  mechind = find(strcmp({allmechs.label},newmech));
  if isempty(mechind) && exist(fullfile(pwd,[newmech '.txt.']),'file')
    file=fullfile(pwd,[newmech '.txt.']);
    this = parse_mech_spec(file,[]);
    cfg.allmechfiles{end+1}=file;
    this.label = newmech;
    this.file = file;
    allmechs(end+1)=this;
    newmech=this;
    mechind=length(allmechs);
  else
    newmech = allmechs(mechind);
  end
  if isempty(newmech)
    warndlg([mechadded{i} ' not found. Check spelling and case.']);
    disp('known mechanisms include: '); disp(get_mechlist');
  else
    newmech = rmfield(newmech,'file');
    if ~isempty(CURRSPEC.(type)(index)) && isfield(CURRSPEC.(type)(index),'mechs') && isstruct(CURRSPEC.(type)(index).mechs)
      CURRSPEC.(type)(index).mechs(end+1)=newmech;
    else
      CURRSPEC.(type)(index).mechs = newmech;
    end
    CURRSPEC.(type)(index).mechanisms{end+1}=mechadded{i};
    CURRSPEC.files{end+1} = cfg.allmechfiles{mechind};
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatemodel(newspec) % maybe use same specification as for "override"
% purpose: update the integratable model after each change to its specification
global CURRSPEC LASTSPEC cfg
LASTSPEC = CURRSPEC;
CURRSPEC = newspec;
if 1
  % remove all connections for which gCOM=0 before building the model
  % ...
  [model,IC,functions,auxvars,CURRSPEC,sodes,svars,txt] = buildmodel2(CURRSPEC,'verbose',0);
  if isfield(CURRSPEC,'entities') && ~isfield(CURRSPEC,'cells')
    CURRSPEC.cells=CURRSPEC.entities; CURRSPEC=rmfield(CURRSPEC,'entities'); 
  end
  cfg.modeltext = txt;
  h=findobj('tag','modeltext');
  if ~isempty(h)
    set(h,'string',cfg.modeltext);
  end
  %fprintf('model updated successfully\n');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undo(src,evnt)
% revert to the last working model
global LASTSPEC
updatemodel(LASTSPEC);
SelectCells;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refresh(src,evnt,where)
if 1
  % load all mech data (if db has changed)
  global allmechs cfg
  [allmechlist,allmechfiles]=get_mechlist(cfg.DBPATH);
  if ~isequal(cfg.allmechfiles_db,allmechfiles)
    cfg.allmechfiles_db=allmechfiles;
    cfg.allmechfiles=unique({cfg.allmechfiles{:},allmechfiles{:}});
    for i=1:length(cfg.allmechfiles)
      this = parse_mech_spec(cfg.allmechfiles{i},[]);
      [fpath,fname,fext]=fileparts(cfg.allmechfiles{i});
      this.label = fname;
      this.file = cfg.allmechfiles{i};
      if i==1, allmechs = this;
      else allmechs(i) = this;
      end
    end
  end
end
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
DrawAuxView;
DrawUserParams([],[],[],0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printmodel(src,evnt)
global CURRSPEC
buildmodel2(CURRSPEC,'verbose',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=get_neighb(n,xi,yi)
% purpose: get neighborhood of a given grid element
% N = linear identifiers for neighbors of compartment at (x,y) in grid
% n = # of elements in grid
% x,y = row,col indices (increasing away from top-left corner)
id=flipud(reshape(1:n,sqrt([n n])));
N=[xi-1,yi-1; xi-1,yi; xi-1,yi+1; xi,yi+1; xi+1,yi+1; xi+1,yi; xi+1,yi-1; xi,yi-1];
N(N<1)=1; N(N>sqrt(n))=sqrt(n);
N=cellfun(@(x,y)id(x,y),num2cell(N(:,1)),num2cell(N(:,2)));
N(N==id(xi,yi))=[];
[jnk,I]=unique(N,'first');
N=N(sort(I));
function conninds=get_connected(spec,compind)
% purpose: get compartments connected to a target compartment
% conninds = linear indices in spec.connections(conninds) of cell/comp types
%            connected to spec.cells(compind).
conninds=[];
if size(spec.connections,1)>=compind
  conninds=[conninds find(~cellfun(@isempty,{spec.connections(compind,:).label}))];
end
if size(spec.connections,2)>=compind
  conninds=[conninds find(~cellfun(@isempty,{spec.connections(:,compind).label}))];
end
conninds=unique(conninds);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function changename(src,evnt)
global currspec
name=get(src,'string');
[currspec.cells.parent]=deal(name);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function seladj(src,evnt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ZoomFunction(src,evnt)
global H
hplot=H.ax_static_plot;
YLIM = get(hplot,'ylim'); 
if isempty(YLIM) || YLIM(1)==YLIM(2), return; end
if evnt.VerticalScrollCount < 0           % zoom in
  set(hplot,'ylim',YLIM/1.5);
else                                      % zoom out
  set(hplot,'ylim',YLIM*1.5);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawAuxFunctions(src,evnt)
global H cfg CURRSPEC
% get list of functions in focuscomp
maxlhs = 20; % limit how much is shown in the listbox
maxlen = 150;
funcs = CURRSPEC.cells(cfg.focuscomp).functions;
len = min(maxlhs,max(cellfun(@length,funcs(:,1))));
str = {};
for i=1:size(funcs,1)
  str{i} = sprintf(['%-' num2str(len) 's  = %s'],funcs{i,1},strrep(funcs{i,2},' ',''));
  if length(str{i})>maxlen, str{i}=str{i}(1:maxlen); end
end
% get list of functions in aux listbox
sel = get(H.lst_static_funcs,'value');
list = get(H.lst_static_funcs,'string');
if ~isequal(list,str)
  sel=sel(sel<=length(str));
  set(H.lst_static_funcs,'string',str);
  set(H.lst_static_funcs,'value',sel);
  list=str;
end
functions = CURRSPEC.cells(cfg.focuscomp).functions(sel,:);
% manage curves
if isfield(H,'static_traces')
  axislimits=[get(H.ax_static_plot,'xlim') get(H.ax_static_plot,'ylim')];
  try delete(H.static_traces); end
  H=rmfield(H,'static_traces'); 
  cla(H.ax_static_plot);
else
  axislimits='tight';
end
% only consider functions of one variable
keep = zeros(size(sel));
for k=1:length(sel)
  var = regexp(functions{k,2},'@\([\w,]+\)','match'); var=var{1};
  if length(splitstr(var,','))==1
    keep(k)=1;
  end
end
functions = functions(keep==1,:);
X = cfg.V; 
for k=1:size(functions,1)
  f = str2func(functions{k,2});
  try
    Y = f(X);
    H.static_traces(k)=line('parent',H.ax_static_plot,'color',cfg.colors(max(1,mod(k,length(cfg.colors)))),...
      'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',X,'ydata',Y,'zdata',[]);
  end
end
funclabels=cellfun(@(x)strrep(x,'_','\_'),functions(:,1),'uni',0);
if isfield(H,'static_traces') && ~isempty(H.static_traces)
  h=legend(H.ax_static_plot,funclabels); set(h,'fontsize',6,'location','EastOutside');
  if strcmp(axislimits,'tight')
    axes(H.ax_static_plot); axis(axislimits);
  else
    set(H.ax_static_plot,'xlim',axislimits(1:2),'ylim',axislimits(3:4));
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawUserParams(src,evnt,focuscomp,movetabs)
global H cfg CURRSPEC
if nargin<4, movetabs=1; end
if nargin>=3 && ~isempty(focuscomp)
  cfg.focuscomp=focuscomp;
end
this=CURRSPEC.cells(cfg.focuscomp);
if isempty(this), return; end
if movetabs
  % switch to cellview:
  set(findobj('tag','ptoggle'),'visible','off'); 
  set(findobj('tag','tab'),'backgroundcolor',[1 1 1]); 
  set(findobj('userdata','pcell'),'visible','on'); 
  set(H.bcell,'backgroundcolor',[.7 .7 .7]);
end  
  
% set panel title to indicate the current compartment
if isfield(this,'parent') && ~isempty(this.parent), prefix=[this.parent '.']; else prefix=''; end
set(H.p_cell_parms,'Title',[prefix this.label ': override parameters'],'FontWeight','normal');

% set dynamics edit box
set(H.edit_comp_dynamics,'string',[this.dynamics{:}],...
  'Callback',{@UpdateParams,'dynamics','cells'});

% add parameter boxes for global intrinsic and connection parameters

% intrinsic parameters
p = this.parameters;
if ~isempty(p)
  intparms = p(1:2:end); valind=1; val=num2str(p{2});
else
  intparms = ''; valind=[]; val='';
end
H.ui_cells_paramlabel = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','text','position',[0 .9 .4 .05],'string','intrinsic','HorizontalAlignment','left');
H.ui_cells_paramlist = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','listbox','position',[0 .2 .4 .7],'value',valind,'string',intparms,...
  'backgroundcolor','w','Max',1,'Min',0,'Callback',{@UpdateParams,'show','cells'});
H.ui_cells_paramedit = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[0 .1 .4 .05],'backgroundcolor','w','string',val,...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'change','cells'});
H.ui_cells_paramadd = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[0 0 .4 .05],'backgroundcolor','w','string','key = value',...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'add','cells'});

% connection parameters
conninds=get_connected(CURRSPEC,cfg.focuscomp);
connparms = {}; val=''; valind=[];
for i=1:length(conninds)
  this = CURRSPEC.connections(conninds(i),cfg.focuscomp);
  if isempty(this.parameters), continue; end
  if i==1, valind=1; val=num2str(this.parameters{2}); end
  thisp = this.parameters(1:2:end);
  for j=1:length(thisp), thisp{j}=[this.label '.' thisp{j}]; end
  connparms = {connparms{:} thisp{:}};
end
H.ui_connections_paramlabel = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','text','position',[.45 .9 .5 .05],'string','connections','HorizontalAlignment','left');
H.ui_connections_paramlist = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','listbox','position',[.45 .2 .5 .7],'value',valind,'string',connparms,...
  'backgroundcolor','w','Max',1,'Min',0,'Callback',{@UpdateParams,'show','connections'});
H.ui_connections_paramedit = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[.45 .1 .5 .05],'backgroundcolor','w','string',val,...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'change','connections'});
H.ui_connections_paramadd = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[.45 0 .5 .05],'backgroundcolor','w','string','key = value',...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'add','connections'});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateParams(src,evnt,action,type) 
% note: type=connections|cells, and is used w/ srcnum to get params from
% cells or connections fields using the same statements.
global cfg H CURRSPEC
h = H.(['ui_' type '_paramlist']);
g = H.(['ui_' type '_paramedit']);
f = H.(['ui_' type '_paramadd']);
newspec = CURRSPEC;
list = get(h,'string');
if isempty(list)
  list = {};
end
if strcmp(type,'connections')
  if strcmp(action,'show') || strcmp(action,'change')
    label = splitstr(list{get(h,'value')},'.');
  elseif strcmp(action,'add')
    label = regexp(get(f,'string'),'^.+=','match');
    label = splitstr(strtrim(label{1}(1:end-1)),'.');
  end
  srcnum = find(cellfun(@(x)isequal(x,label{1}),{CURRSPEC.(type)(:,cfg.focuscomp).label}));
  if isempty(srcnum) && strcmp(action,'add')
    fprintf('Improper connection parameter. Use format: from-to.param=value\n');
    return;
  end
else
  srcnum = 1;
end
keys = newspec.(type)(srcnum,cfg.focuscomp).parameters(1:2:end);
switch action
  case 'show' % get value of selected param from currspec
    values = [newspec.(type)(srcnum,cfg.focuscomp).parameters{2:2:end}];
    if strcmp(type,'connections')
      %length(find(~cellfun(@isempty,regexp(list,['^' label{1} '\.']))))
      sel = strmatch(label{2},newspec.(type)(srcnum,cfg.focuscomp).parameters(1:2:end),'exact');
    else
      sel = get(h,'value');
    end
    set(g,'string',num2str(values(sel(1))));
  case 'change' % update currspec param value w/ value of selected param
    pause(.5); % pause to let gui catch up; otherwise get() will retrieve the previous value and not the new one
    key = list{get(h,'value')};
    if strcmp(type,'connections')
      key = splitstr(key,'.');
      key = key{2};
    end
    val = str2num(get(g,'string'));
    ind = 2*find(strcmp(key,keys));
    newspec.(type)(srcnum,cfg.focuscomp).parameters{ind} = val;
    updatemodel(newspec);
  case 'add' % parse paramadd string for param and value to add to currspec
    str = get(f,'string');
    key = regexp(str,'^.+=','match');
    key = strtrim(key{1}(1:end-1));
    if strcmp(type,'connections')
      key = splitstr(key,'.');
      src=key{1}; key=key{2};
    end
    val = regexp(str,'=.+$','match');
    val = strtrim(val{1}(2:end));
    val = str2num(val);
    if any(~ismember(key,keys)) % add if a new parameter
      newspec.(type)(srcnum,cfg.focuscomp).parameters{end+1} = key;
      newspec.(type)(srcnum,cfg.focuscomp).parameters{end+1} = val;
    else % change value if param already in currspec
      ind = 2*find(strcmp(key,keys));
      newspec.(type)(srcnum,cfg.focuscomp).parameters{ind} = val;
    end
    if strcmp(type,'connections')
      list = {list{:}, [src '.' key]};
    else
      list = {list{:}, key};
    end
    %list = newspec.(type)(srcnum,cfg.focuscomp).parameters(1:2:end);
    set(h,'string',list);
    set(h,'value',find(strcmp(key,list)));
    set(g,'string',num2str(val));
    updatemodel(newspec);
  case 'dynamics'
    newspec.(type)(cfg.focuscomp).dynamics=splitstr(get(src,'string'),',');
    updatemodel(newspec);
end

%% add component labels to boxes (i.e., "draw" the multi-compartment cell)
%  - i.e., determine the arrangement of compartments from an arbitrary spec
%          structure (which does not contain compartment coordinate info).
% TODO: rewrite this block (it's ugly and will probably break for complex morphologies).
%       ex) have get_neighb() return linear indices (not identifiers) to grid elements...


% jobj=findjobj(H.ui_cellcomp(xi,yi));
% set(jobj,'MouseEnteredCallback',{@Update,'controls',compnames{compcnt}});

% con=currspec.connections; 
% [from,to]=ind2sub(size(con),cellfun(@isempty,{con.mechanisms}));
% %from=find(from); to=find(to);
% root=currspec.cells(1).label;
% level = nan(size(currspec.cells));
% level(1)=0;
% bucket=[1];%{root};
% parents=[0]; parentsi=[1];
% while ~isempty(bucket)
%      parent = bucket(1);%{1};
%      bucket(1)=[];
%      children = find((from(parent)&to)|(to(parent)&from))
%      children = {x | parent=>x or x=>parent} = nodes connected to parent
%      foreach child: level['child']=level['parent']+1
%      bucket=cat(2,bucket,children)
%      parents=cat(2,parents,repmat(parent,[1 numel(children)]))
% end
% for i=1:max(level)
%      nodes= {nodes w/ level==i}
%      for j=1:length(nodes)
%                goto parent of nodes(j)
%                get neighborhood of parent (set of surrounding grid-cells)
%                opt=find neighborhood cells surrounded by empty grid-cells
%                if ~isempty(opt)
%                     put nodes(j) in opt(1)
%                else
%                     put nodes(j) in first empty neighborhood grid-cell
%                end
%      end
% end

% HOW TO FILL CELL MORPHOLOGY GRID WITH COMPARTMENT LABELS:
% root=cells(1).label
% level = nan(#nodes)
% level[root]=0
% bucket={root}
% parents=[0]
% while ~isempty(bucket)
%      parent = bucket{1}
%      bucket(1)=[]
%      children = {x | parent=>x or x=>parent} = nodes connected to parent
%      foreach child: level['child']=level['parent']+1
%      bucket=cat(2,bucket,children)
%      parents=cat(2,parents,repmat(parent,[1 numel(children)]))
% end
% for i=1:max(level)
%      nodes= {nodes w/ level==i}
%      for j=1:length(nodes)
%                goto parent of nodes(j)
%                get neighborhood of parent (set of surrounding grid-cells)
%                opt=find neighborhood cells surrounded by empty grid-cells
%                if ~isempty(opt)
%                     put nodes(j) in opt(1)
%                else
%                     put nodes(j) in first empty neighborhood grid-cell
%                end
%      end
% end
% 
