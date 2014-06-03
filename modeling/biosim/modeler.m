function modeler(varargin)
clear global cfg H CURRSPEC
global cfg H CURRSPEC LASTSPEC BIOSIMROOT %currspec lastspec
warning off
prepare_spec;
updatemodel(CURRSPEC);

%% set up GUI

% main figure
sz = get(0,'ScreenSize'); sz0=sz;
sz = [.005*sz(3) .07*sz(4) .938*sz(3) .83*sz(4)];
% sz = [10 300 1350 700];
fig = figure('position',sz,'color','w','name','contact: jason@infinitebrain.org','NumberTitle','off','WindowScrollWheelFcn',@ZoomFunction,'CloseRequestFcn','delete(gcf); clear global H'); % [320 240 920 560]

% global controls (i.e., always present in main figure in all views)
titlestring = 'Dynamic Neural Simulator'; % DNSim
username = 'nak-crc';
  uicontrol('parent',fig,'style','text','string',titlestring,'fontsize',19,'units','normalized','position',[.08 .895 .25 .07],'backgroundcolor','w');
  uicontrol('parent',fig,'style','text','string',['user: ' username],'fontsize',12,'units','normalized','position',[.08 .86 .25 .07],'backgroundcolor','w');
% tabs:
  bbuild=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[0 .85 .1 .04],'string','build','backgroundcolor',[.7 .7 .7],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pbuild''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bmodel=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[.1 .85 .1 .04],'string','model','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pmodel''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bsimstudy=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[.2 .85 .1 .04],'string','batch','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''psimstudy''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bhistory=uicontrol('parent',fig,'style','pushbutton','tag','tab','units','normalized','position',[.3 .85 .1 .04],'string','history','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle''),''visible'',''off''); set(findobj(''tag'',''tab''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''phistory''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
% model controls:
  bsave=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[0 .98 .035 .025],'string','save','backgroundcolor',[.8 .8 .8],'callback',@Save_Spec);
  bload=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.035 .98 .035 .025],'string','load','backgroundcolor',[.8 .8 .8],'callback',@Load_File);
  bappend=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.07 .98 .04 .025],'string','append','backgroundcolor',[.8 .8 .8],'callback',{@Load_File,[],1});
  %bclear=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.098 .98 .02 .025],'string','x','backgroundcolor',[.8 .8 .8],'callback',{@SelectCells,1});
  bundo=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.11 .98 .035 .025],'string','undo','backgroundcolor',[.8 .8 .8],'callback',@undo);
  bapply=uicontrol('parent',fig,'style','pushbutton','units','normalized','position',[.35 .98 .05 .025],'string','refresh','backgroundcolor',[.8 .8 .8],'callback',{@refresh,0},'visible','off');

% left panels for cell, network, and mechanism controls
pbuild=uipanel('parent',fig,'title','model builder','visible','on','tag','ptoggle','userdata','pbuild','units','normalized','position',[0 0 .4 .85]);
  bmech=uicontrol('parent',pbuild,'style','pushbutton','tag','tab2','units','normalized','position',[.2 .65 .2 .04],'string','mechanisms','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle2''),''visible'',''off''); set(findobj(''tag'',''tab2''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pmech''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bnet=uicontrol('parent',pbuild,'style','pushbutton','tag','tab2','units','normalized','position',[.4 .65 .2 .04],'string','connections','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle2''),''visible'',''off''); set(findobj(''tag'',''tab2''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pnet''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  bcell=uicontrol('parent',pbuild,'style','pushbutton','tag','tab2','units','normalized','position',[.6 .65 .2 .04],'string','parameters','backgroundcolor',[1 1 1],'callback','set(findobj(''tag'',''ptoggle2''),''visible'',''off''); set(findobj(''tag'',''tab2''),''backgroundcolor'',[1 1 1]); set(findobj(''userdata'',''pcell''),''visible'',''on''); set(gcbo,''backgroundcolor'',[.7 .7 .7]);');
  pmech=uipanel('parent',pbuild,'title','mechanism editor','visible','off','tag','ptoggle2','userdata','pmech','units','normalized','position',[0 0 1 .65]);
  pnet=uipanel('parent',pbuild,'title','network editor','visible','off','tag','ptoggle2','userdata','pnet','units','normalized','position',[0 0 1 .65]);
  pcell=uipanel('parent',pbuild,'title','compartment/entity/connection parameters','visible','off','tag','ptoggle2','userdata','pcell','units','normalized','position',[0 0 1 .65]);
pmodel=uipanel('parent',fig,'title','review model','visible','off','tag','ptoggle','userdata','pmodel','units','normalized','position',[0 0 .4 .85]);
psimstudy=uipanel('parent',fig,'title','batch simulation','visible','off','tag','ptoggle','userdata','psimstudy','units','normalized','position',[0 0 .4 .85]);
  pbatchcontrols=uipanel('parent',psimstudy,'title','','units','normalized','position',[0 .8 1 .2]);
  pbatchspace=uipanel('parent',psimstudy,'title','search space','units','normalized','position',[0 .3 1 .5]);
  pbatchoutputs=uipanel('parent',psimstudy,'title','outputs','units','normalized','position',[0 0 1 .3]);
phistory=uipanel('parent',fig,'title','','visible','off','tag','ptoggle','userdata','phistory','units','normalized','position',[0 0 .4 .85]);
  pnotes=uipanel('parent',phistory,'title','project notes','units','normalized','position',[.22 .3 .78 .7]);%,'backgroundcolor',[.5 .5 .5]);
  pcomparison=uipanel('parent',phistory,'title','model comparison','units','normalized','position',[.22 0 .78 .3]);

% make notes scrollable:
%   pnotes=uipanel('parent',phistory,'title','project notes','units','normalized','position',[0 0 1 1]);%,'backgroundcolor',[.5 .5 .5]);
%   jEdit = findjobj(pnotes);
%   j=javax.swing.JScrollPane(jEdit);
%   [jj hh]=javacomponent(j,[1 1 .4 .85].*[.22 1 .78 .3].*[sz(3) sz(4) sz(3) sz(4)],fig);
%   set(pnotes,'position',[.22 .3 .78 .7]);
  
% left panel: model view
txt_model = uicontrol('parent',pmodel,'style','edit','units','normalized','tag','modeltext',...
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
p_net_select  = uipanel('parent',pbuild,'Position',[0 .7 1 .29],'BorderWidth',.2,'BorderType','line'); % cell morphology
p_net_connect = uipanel('parent',pnet,'Position',[0 .5 1 .5],'BorderWidth',.2,'BorderType','line','title','Population connections         [targets]'); % cell specification
p_net_kernel  = uipanel('parent',pnet,'Position',[0 0 1 .5],'BorderWidth',.2,'BorderType','line','title','Cell connections'); % cell specification
% compartment controls
if ~isempty(net.cells) && ischar(net.cells(1).label)
  l1={net.cells.label}; 
  %l2={net.cells.parent};
  %l=cellfun(@(x,y)[x '.' y],l2,l1,'uni',0);
  l=l1;
  i=1:length(l);
else
  l={}; i=[];
end
lst_comps = uicontrol('parent',p_net_select,'units','normalized','style','listbox','position',[0 0 .2 .9],'value',i,'string',l,'BackgroundColor',[.9 .9 .9],'Max',5,'Min',0,'Callback',@SelectCells,...
  'ButtonDownFcn',@RenameComponent,'TooltipString','Right-click to edit compartment name');
% headers for cell info
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[0 .91 .25 .09],'string','Compartments','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);
%uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.25 .91 .1 .09],'string','name','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.25 .91 .06 .09],'string','n','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.59 .91 .15 .09],'string','mechanisms','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);
uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.31 .91 .13 .09],'string','dynamics','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);
% uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.41 .91 .15 .09],'string','mechanisms','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);
% uicontrol('parent',p_net_select,'units','normalized','style','text','position',[.79 .91 .13 .09],'string','dynamics','ListboxTop',0,'HorizontalAlignment','left','fontsize',10);

% left panel: mechanism editor %GUI_mechpanel;
% compartment label
if ~isempty(CURRSPEC.(cfg.focustype))
  tmp=CURRSPEC.(cfg.focustype)(cfg.focus);
  str1=tmp.mechanisms;
  str2=mech_spec2str(tmp.mechs(cfg.focusmech));
  cl=tmp.label;
  u.focustype=cfg.focustype; 
  u.focus=cfg.focus; 
  u.mechlabel=[tmp.mechanisms{cfg.focusmech}]; 
else
  str1='';
  str2='';
  cl='';
  u=[];
end
txt_comp = uicontrol('style','text','string',cl,'units','normalized','position',[.05 .95 .1 .05],'parent',pmech,'FontWeight','bold','visible','off');
% dropdown list of mechanisms for this compartment
% lst_mechs = uicontrol('style','popupmenu','value',min(1,length(str1)),'string',str1,...
%   'units','normalized','position',[.11 .95 .2 .05],'parent',pmech,'callback',@Display_Mech_Info);
% button to apply changes to mech text
uicontrol('parent',pmech,'style','pushbutton','units','normalized','position',[.75 .97 .1 .04],'string','apply','callback',@UpdateMech,'visible','off');
% button to create a new mechanism
uicontrol('parent',pmech,'style','pushbutton','units','normalized','position',[.85 .97 .1 .04],'string','write','callback',@SaveMech);
% button to display list of mechs in DB
uicontrol('parent',pmech,'style','pushbutton','units','normalized','position',[.95 .97 .05 .04],'string','DB','callback','global allmechs; msgbox({allmechs.label},''available'');');%msgbox(get_mechlist,''available'')');%'get_mechlist');
% edit box with mech info
lst_mechs = uicontrol('units','normalized','position',[0 .42 .2 .55],'parent',pmech,'BackgroundColor',[.9 .9 .9],...
  'style','listbox','value',1,'string',str1,'Max',1,'Callback',@Display_Mech_Info,'ButtonDownFcn',@RenameMech,'TooltipString','Right-click to edit mechanism name');
txt_mech = uicontrol('parent',pmech,'style','edit','units','normalized','BackgroundColor','w','callback',@UpdateMech,... % [.9 .9 .9]
  'position',[.2 .42 .8 .55],'string',str2,'userdata',u,'FontName','courier','FontSize',10,'HorizontalAlignment','Left','Max',100);
% mech plots associated w/ this compartment
p_static_plots = uipanel('parent',pmech,'Position',[0 0 1 .4],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','');
lst_static_funcs = uicontrol('units','normalized','position',[0 0 .2 .95],'parent',p_static_plots,'BackgroundColor',[.9 .9 .9],...
  'style','listbox','value',1:5,'string',{},'Max',50,'Callback',@DrawAuxFunctions);
ax_static_plot = subplot('position',[.2 0 .8 .95],'parent',p_static_plots,'linewidth',3,'color','w','fontsize',6); box on; 
% lst_static_funcs = uicontrol('units','normalized','position',[.04 .02 .9 .35],'parent',p_static_plots,...
%   'style','listbox','value',1:5,'string',{},'Max',50,'Callback',@DrawAuxFunctions);
edit_static_lims=uicontrol('Style','edit', 'Units','normalized','Position',[0.85 0.075 0.13 0.1],...
          'String',sprintf('[%g,%g]',min(cfg.V),max(cfg.V)),'Callback',{@DrawAuxFunctions,1},'parent',p_static_plots);
if ~isempty(CURRSPEC.cells)
  maxlhs=20; maxlen=150; % limit how much is shown in the listbox
  funcs = CURRSPEC.cells(cfg.focuscomp).functions;
  len = min(maxlhs,max(cellfun(@length,funcs(:,1))));
  str = {};
  for i=1:size(funcs,1)
    str{i} = sprintf(['%-' num2str(len) 's  = %s'],funcs{i,1},strrep(funcs{i,2},' ',''));
    if length(str{i})>maxlen, str{i}=str{i}(1:maxlen); end
  end
  val=1:4;
  set(lst_static_funcs,'string',str);
  set(lst_static_funcs,'value',val(val<=length(str)));
end

% left panel: cell builder %GUI_cellpanel;
% p_cell_morph = uipanel('parent',pcell,'Position',[0 .6 .35 .35],'BorderWidth',.2,'BorderType','line','visible','off'); % cell morphology
% p_cell_parms = uipanel('parent',pcell,'Position',[0 0 1 1],'BorderWidth',.2,'BorderType','line','title','Parameters');
% p_cell_spec = uipanel('parent',pcell,'Position',[0 0 1 .4],'BorderWidth',.2,'BorderType','line','title','as','fontangle','italic'); % cell specification
% edit_comp_dynamics = uicontrol('parent',pcell,'style','edit','string','',...
%     'units','normalized','position',[.1 .1 .8 .05],'BackgroundColor','w','HorizontalAlignment','left','tooltipstring','Mechanism functions will be substituted here.');
  
% right panel: simulation plots and controls %GUI_simpanel;
psims=uipanel('parent',fig,'title','','visible','on','units','normalized','position',[.4 0 .6 1]);

% menu %GUI_menu;
set(fig,'MenuBar','none');
file_m = uimenu(fig,'Label','File');
uimenu(file_m,'Label','Load model','Callback',@Load_File);
% load_m = uimenu(file_m,'Label','Load');
% uimenu(load_m,'Label','model','Callback',@Load_Spec);
% uimenu(load_m,'Label','sim_data','Callback',@Load_File);
uimenu(file_m,'Label','Save model','Callback',@Save_Spec);
ws_m = uimenu(file_m,'Label','Interact');
uimenu(ws_m,'Label','Pass model (''spec'') to command window','Callback','global CURRSPEC; assignin(''base'',''spec'',CURRSPEC);');
uimenu(ws_m,'Label','Update model (''spec'') from command window','Callback',{@refresh,1});
uimenu(ws_m,'Label','Pass ''sim_data'' (during interactive simulation) to command window','Callback','global cfg;cfg.publish=1;');
import_m = uimenu(file_m,'Label','Import');
uimenu(import_m,'Label','XPP (wip)','Callback','not implemented yet');
export_m = uimenu(file_m,'Label','Export');
uimenu(export_m,'Label','XPP (wip)','Callback','not implemented yet');
uimenu(export_m,'Label','NEURON (wip)','Callback','not implemented yet');
uimenu(export_m,'Label','CellML (wip)','Callback','not implemented yet');
uimenu(file_m,'Label','Exit','Callback','global CURRSPEC H cfg; close(H.fig); clear CURRSPEC H cfg; warning on');
plot_m = uimenu(fig,'Label','Plot');
uimenu(plot_m,'Label','plotv','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), plotv(evalin(''base'',''sim_data''),CURRSPEC); else disp(''load data to plot''); end');
uimenu(plot_m,'Label','plotpow','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), plotpow(evalin(''base'',''sim_data''),CURRSPEC,''spectrogram_flag'',0); else disp(''load data to plot''); end');
uimenu(plot_m,'Label','plotspk','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), plotspk(evalin(''base'',''sim_data''),CURRSPEC,''window_size'',30/1000,''dW'',5/1000); else disp(''load data to plot''); end');
uimenu(plot_m,'Label','visualizer','Callback','global CURRSPEC; if ismember(''sim_data'',evalin(''base'',''who'')), visualizer(evalin(''base'',''sim_data'')); else disp(''load data to plot''); end');

% collect object handles
% figures
H.fig   = fig;
% panels
H.pbuild = pbuild;
H.pcell  = pcell;
H.pnet  = pnet;
H.pmech = pmech;
H.pmodel = pmodel;
H.psimstudy = psimstudy;
H.phistory = phistory;
H.psims = psims;
H.p_net_select  = p_net_select;
H.p_net_connect = p_net_connect;
H.p_net_kernel  = p_net_kernel;
H.pbatchcontrols = pbatchcontrols;
H.pbatchspace = pbatchspace;
H.pbatchoutputs = pbatchoutputs;
H.pnotes = pnotes;
H.pcomparison = pcomparison;
% H.p_cell_morph  = p_cell_morph;
% H.p_cell_spec   = p_cell_spec;
% H.p_cell_parms  = p_cell_parms;
H.p_static_plots= p_static_plots;
H.edit_static_lims = edit_static_lims;
% buttons
H.pbuild = pbuild; 
H.bnet  = bnet;
H.bmech = bmech;
H.bmodel = bmodel; % @callback (name callback functions here...)
H.bsimstudy = bsimstudy;
H.bhistory = bhistory;
H.bsave = bsave; % @callback
H.bload = bload; % @callback
H.bapply= bapply;
H.bappend= bappend;
H.bundo = bundo;
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
% H.edit_comp_dynamics = edit_comp_dynamics; % cell dynamics for focuscomp

% populate controls
if ~isempty(CURRSPEC.cells)
  SelectCells;
  DrawAuxView;
  DrawAuxFunctions;
  DrawUserParams;%([],[],[],0);
  DrawSimPlots;
end
DrawStudyInfo;
UpdateHistory;
Display_Mech_Info;

%% FUNCTIONS

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Load_File(src,evnt,path,append_flag)
% Usage:
% @Load_File; % to load and override current data and model
% Load_File([],[],path); % to load from path (e.g., batch result saved in rootdir) and override current data and model
% Load_File([],[],[],1); % to load and concatenate models
% Load_File([],[],path,1); % to load from path and concatenate
if nargin<4, append_flag=0; end
if nargin>=3 && isdir(path) % cd if given a directory where to find files
  cwd=pwd;
  cd(path);
end
[filename,pathname] = uigetfile({'*.mat'},'Pick a model or sim_data file.','MultiSelect','off');
if nargin>=3 && isdir(path) % move back to original directory
  cd(cwd);
end
if isequal(filename,0) || isequal(pathname,0), return; end
if iscell(filename)
  datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
  filename = filename{1};
else
  datafile = [pathname filename];
end
if exist(datafile,'file')
  fprintf('Loading file: %s\n',datafile);
  try
    o=load(datafile); % load file
  catch
    fprintf('failed to load file. check that it is a valid matlab file: %s\n',datafile);
    return;
  end
  if isfield(o,'modelspec') % standardize spec name
    o.spec=o.modelspec;
    o=rmfield(o,'modelspec');
  end
  if ~isfield(o,'sim_data') && ~isfield(o,'spec')
    fprintf('select file does not contain sim_data or spec structure. no data loaded\n');
    return;
  end
  if isfield(o,'sim_data')
    % pass data to base workspace
    assignin('base','sim_data',o.sim_data);
    fprintf('sim_data assigned to base workspace.\n');
  end
  if isfield(o,'spec') % has model
    % standarize spec structure
    if isfield(o.spec,'entities') 
      if ~isfield(o.spec,'cells')
        o.spec.cells=o.spec.entities;
      end
      o.spec=rmfield(o.spec,'entities');
    end
    if ~isfield(o.spec,'history')
      o.spec.history=[];
    end
    if ~isfield(o.spec.cells,'parent')
      for i=1:length(o.spec.cells)
        o.spec.cells(i).parent=o.spec.cells(i).label;
      end      
    end    
    global CURRSPEC
    newspec=CURRSPEC;
    if append_flag
      if isempty(newspec) || isempty(newspec.cells)
        test = 0;
      else
        test = ismember({o.spec.cells.label},{newspec.cells.label});
      end
      if any(test)
        dup={o.spec.cells.label};
        dup=dup(test);
        str='';
        for i=1:length(dup)
          str=[str dup{i} ', '];
        end
        fprintf('failed to concatenate models. duplicate names found: %s. rename and try again.\n',str(1:end-2));
      else
        if isfield(newspec,'cells') && ~isempty(newspec.cells)
          n=length(o.spec.cells);
          [addflds,I]=setdiff(fieldnames(newspec.cells),fieldnames(o.spec.cells));
          [jnk,I]=sort(I);
          addflds=addflds(I);
          for i=1:length(addflds)
            o.spec.cells(1).(addflds{i})=[];
          end
          o.spec.cells=orderfields(o.spec.cells,newspec.cells);
          o.spec.connections=orderfields(o.spec.connections,newspec.connections);
          newspec.cells(end+1:end+n) = o.spec.cells;
          for i=1:n
            newspec.connections(end+i,end+1:end+n) = o.spec.connections(i,:);
          end
          if isfield(o.spec,'files')
            newspec.files = unique({newspec.files{:},o.spec.files{:}});
          end
        else
          newspec = o.spec;
        end
      end
    else
      newspec=o.spec;
    end
    % update model
    updatemodel(newspec);
    refresh;
    % pass model to base workspace
    assignin('base','spec',o.spec);
    fprintf('model specification loaded and assigned to base workspace as ''spec''.\n');    
  end
else
  fprintf('file does not exist: %s\n',datafile);
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
if null==1
  CURRSPEC=[];
  set(H.lst_comps,'string',{},'value',[]);
  set(H.lst_notes,'string',{},'value',[]);
  UpdateHistory;
  try
    set(H.edit_comp_parent,'visible','off');
    set(H.edit_comp_label,'visible','off');
    set(H.edit_comp_N,'visible','off');
    set(H.edit_comp_dynamics,'visible','off');
    set(H.edit_comp_mechs,'visible','off');
    set(H.btn_comp_copy,'visible','off');
    set(H.btn_comp_delete,'visible','off');
    set(H.btn_comp_edit,'visible','off');    
    set(H.p_comp_mechs,'visible','off');    
  end
  %refresh;
  return; 
elseif isempty(CURRSPEC) || isempty(CURRSPEC.cells)
  set(H.lst_comps,'string',{},'value',[]);
  set(H.lst_mechs,'string',{},'value',[]);
  return;
end
v=get(H.lst_comps,'value'); 
l=get(H.lst_comps,'string');
if ischar(CURRSPEC.cells(1).label)
  if isfield(CURRSPEC.cells,'parent') && ~isempty(CURRSPEC.cells(1).parent)
    l1={CURRSPEC.cells.label}; 
    %l2={CURRSPEC.cells.parent};
    %l=cellfun(@(x,y)[x '.' y],l2,l1,'uni',0);
    l=l1;
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

% cfg.pauseflag=1;
% DrawSimPlots;
% cfg.pauseflag=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawCellInfo(net)
global H
c=1.5; dy=-.07*c; ht=.1;
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
  if isempty(m)
    str='';
  else
    str=m{1}; for j=2:length(m), str=[str ', ' m{j}]; end
  end
  if ~isfield(H,'edit_comp_label') || length(H.edit_comp_label)<length(sel) || ~ishandle(H.edit_comp_label(i))
    H.edit_comp_parent(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.24 .8+dy*(i-1) .09 ht],'backgroundcolor','w','string',p{i},...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'parent'},'ButtonDownFcn',{@DrawUserParams,sel(i)},'visible','off');    
    H.edit_comp_label(i) = uicontrol('parent',H.p_net_select,'units','normalized','visible','off',...
      'style','edit','position',[.24 .8+dy*(i-1) .1 ht],'backgroundcolor','w','string',l{i},...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'label'},'ButtonDownFcn',{@DrawUserParams,sel(i)});
    H.btn_comp_delete(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','-','callback',{@DeleteCell,l{i}},...
      'position',[.205 .8+dy*(i-1) .03 ht],'TooltipString',l{i});%,'BackgroundColor','white');    
    H.edit_comp_N(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.24 .8+dy*(i-1) .06 ht],'backgroundcolor','w','string',N(i),...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'multiplicity'},'TooltipString',l{i});
    H.edit_comp_dynamics(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.3 .8+dy*(i-1) .28 ht],'backgroundcolor','w','string',[net.cells(sel(i)).dynamics{:}],...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'dynamics'},...
      'ButtonDownFcn',{@Display_Mech_Info,l{i},{},'cells'},'fontsize',9,'TooltipString',l{i});    
    H.edit_comp_mechs(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','edit','position',[.58 .8+dy*(i-1) .38 ht],'backgroundcolor','w','string',str,...
      'HorizontalAlignment','left','Callback',{@UpdateCells,l{i},'mechanisms'},...
      'ButtonDownFcn',{@Display_Mech_Info,l{i},{},'cells'},'fontsize',9,'TooltipString',l{i});
    H.p_comp_mechs(i) = uipanel('parent',H.p_net_select,'units','normalized',...
      'position',[.51 .8+dy*(i-1) .42 ht],'visible','off');    
    H.btn_comp_copy(i) = uicontrol('parent',H.p_net_select,'units','normalized',...
      'style','pushbutton','fontsize',10,'string','+','callback',{@CopyCell,l{i}},...
      'position',[.965 .8+dy*(i-1) .03 ht],'TooltipString',l{i});%,'BackgroundColor','white');    
    H.btn_comp_edit(i) = uicontrol('parent',H.p_net_select,'units','normalized','visible','off',...
      'style','pushbutton','fontsize',10,'string','...','callback',{@ShowClickMechList,i,'cells'},...%{@OpenCellModeler,l{i}},...
      'position',[.965 .8+dy*(i-1) .03 ht]);%,'BackgroundColor','white');            
  else
    % update properties
    set(H.edit_comp_parent(i),'string',p{i},'visible','off','Callback',{@UpdateCells,l{i},'parent'});
    set(H.edit_comp_label(i),'string',l{i},'visible','off','Callback',{@UpdateCells,l{i},'label'});
    set(H.edit_comp_dynamics(i),'string',[net.cells(sel(i)).dynamics{:}],'visible','on','Callback',{@UpdateCells,l{i},'dynamics'},'TooltipString',l{i});
    set(H.edit_comp_N(i),'string',N(i),'visible','on','Callback',{@UpdateCells,l{i},'multiplicity'},'TooltipString',l{i});
    set(H.edit_comp_mechs(i),'string',str,'visible','on','Callback',{@UpdateCells,l{i},'mechanisms'},'ButtonDownFcn',{@Display_Mech_Info,l{i},{},'cells'},'TooltipString',l{i});
    set(H.btn_comp_copy(i),'callback',{@CopyCell,l{i}},'visible','on','TooltipString',l{i});
    set(H.btn_comp_delete(i),'callback',{@DeleteCell,l{i}},'visible','on','TooltipString',l{i});
    set(H.btn_comp_edit(i),'callback',{@ShowClickMechList,i,'cells'},'visible','off');
    set(H.p_comp_mechs(i),'visible','off');    
  end
  if length(H.edit_comp_label)>i
    set(H.edit_comp_parent(i+1:end),'visible','off');
    set(H.edit_comp_label(i+1:end),'visible','off');
    set(H.edit_comp_dynamics(i+1:end),'visible','off');
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
dx=.15; x=.13; c=1.5; dy=-.07*c; ht=.1;
sel = get(H.lst_comps,'value');
net.cells = net.cells(sel);
net.connections = net.connections(sel,:);
net.connections = net.connections(:,sel);
l={net.cells.label}; 
for i=1:length(sel)
  for j=1:length(sel)
    m = net.connections(j,i).mechanisms;
    pos = [x+dx*(i-1) .8+dy*(j-1) .9*dx ht];
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
          'style','text','position',[x+dx*(j-1) .91 .11 ht],'string',l{j},...
          'callback',{@ShowClickMechList,this,'connections'});
      end
      if j==1 % from
        this=ones(1,max(sel));
        this(sel)=i;
        H.txt_from(i) = uicontrol('parent',H.p_net_connect,'units','normalized',...
          'style','text','position',[.01 .8+dy*(i-1) .11 ht],'string',l{i},...
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

% make list of connections
EL={CURRSPEC.cells.label};
CL={CURRSPEC.connections.label};
sel=find(~cellfun(@isempty,CL));
if isempty(sel) % no connections
  return; 
end
connlabels={};
for i=1:length(sel)
  m=CURRSPEC.connections(sel(i)).mechanisms;
  if ~iscell(m), m={m}; end
  for j=1:length(m)
    connlabels{end+1} = [CL{sel(i)} '.' m{j}];
  end  
end
connected=sel;

% make list of auxvars for the first connection
sel=strmatch([strrep(CL{sel(1)},'-','_') '_'],CURRSPEC.model.auxvars(:,1));
auxlist=CURRSPEC.model.auxvars(sel,1);
if isempty(auxlist), return; end
% get auxvar matrix 
a=CURRSPEC.model.auxvars(sel,:);
for i=1:size(a,1)
  key = a{i,1};
  try
    eval(sprintf('%s = %s;',a{i,1},a{i,2}));
    val = eval(a{i,2});
  catch
    val = nan;
  end
  userdata.(key).matrix = val;
  userdata.(key).equation = a{i,2};%[a{i,1} ' = ' a{i,2}];
end
auxeqn=[a{i,1} ' = ' a{i,2}];
auxmat=userdata.(auxlist{end}).matrix;
lims=[min(auxmat(:)) max(auxmat(:))];%[.9*min(auxmat(:)) 1.1*max(auxmat(:))];

% create controls
H.lst_conns = uicontrol('parent',H.p_net_kernel,'units','normalized','style','listbox','userdata',connected,...
  'position',[0 0 .2 1],'value',1,'string',connlabels,'BackgroundColor',[.9 .9 .9],'Max',1,'Callback',@UpdateAuxList);
% listbox to select which auxvar to plot
H.lst_auxvars = uicontrol('parent',H.p_net_kernel,'units','normalized','backgroundcolor',[.9 .9 .9],...
  'style','listbox','position',[.21 .5 .3 .5],'string',auxlist,'value',length(auxlist),'userdata',userdata,...
  'Callback',@UpdateAuxPlot); % 'backgroundcolor','w',
% edit box to modify defining expressions
H.txt_auxvar_eqn = uicontrol('parent',H.p_net_kernel,'units','normalized',...
  'style','text','position',[.21 .31 .3 .19],'backgroundcolor',[.9 .9 .9],'string',auxeqn,...
  'HorizontalAlignment','left','fontsize',8);%,'Callback',{@UpdateParams,'change','cells'});

% create button group for predefined adjacency matrices
if 0
  H.rad_adj = uibuttongroup('visible','off','units','normalized','Position',[.21 0 .3 .3],'parent',H.p_net_kernel);
  H.rad_adj_1 = uicontrol('Style','radiobutton','String','1-to-1','units','normalized',...
      'pos',[.05 .75 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
  H.rad_adj_2 = uicontrol('Style','radiobutton','String','all-to-all','units','normalized',...
      'pos',[.05 .45 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
  H.rad_adj_3 = uicontrol('Style','radiobutton','String','random','units','normalized',...
      'pos',[.05 .15 .8 .2],'parent',H.rad_adj,'HandleVisibility','off');
  set(H.rad_adj,'SelectionChangeFcn',@seladj);
  set(H.rad_adj,'SelectedObject',[]);  % No selection
  set(H.rad_adj,'Visible','on');
end

% plots
H.ax_conn_img = subplot('position',[.55 0 .45 1],'parent',H.p_net_kernel); 
H.img_connect = imagesc(auxmat); %axis xy;
if lims(2)>lims(1), caxis(lims); end
colorbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateAuxList(src,evnt)
global CURRSPEC H
if isempty(CURRSPEC.cells), return; end
sel=get(H.lst_conns,'userdata');
if isempty(sel), return; end
CL={CURRSPEC.connections.label};
% make list of auxvars for the first connection
sel=strmatch([strrep(CL{sel(get(H.lst_conns,'value'))},'-','_') '_'],CURRSPEC.model.auxvars(:,1));
auxlist=CURRSPEC.model.auxvars(sel,1);
if isempty(auxlist), return; end
% get auxvar matrix 
a=CURRSPEC.model.auxvars(sel,:);
for i=1:size(a,1)
  key = a{i,1};
  try
    eval(sprintf('%s = %s;',a{i,1},a{i,2}));
    val = eval(a{i,2});
  catch
    val = nan;
  end
  userdata.(key).matrix = val;
  userdata.(key).equation = a{i,2};% [a{i,1} ' = ' a{i,2}];
end
set(H.lst_auxvars,'string',auxlist,'value',length(auxlist),'userdata',userdata);
UpdateAuxPlot;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateAuxPlot(src,evnt)
global H
u=get(H.lst_auxvars,'userdata');
s=get(H.lst_auxvars,'string');
v=get(H.lst_auxvars,'value');
auxmat=u.(s{v}).matrix;
auxeqn=u.(s{v}).equation;
lims=[min(auxmat(:)) max(auxmat(:))];
set(H.img_connect,'cdata',auxmat); axis tight;
if lims(2)>lims(1), caxis(lims); end
set(H.txt_auxvar_eqn,'string',auxeqn);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulate(src,evnt,action)
global CURRSPEC H cfg eventdata t
clear global eventdata
if isequal(src,findobj('tag','start'))
  cfg.quitflag=-1;
end
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
  if cfg.publish~=0
    tmp_data = ts_matrix2data(cfg.record,'sfreq',1/cfg.dt);
    [tmp_data.sensor_info.label] = deal(allvars{:});
    assignin('base','sim_data',tmp_data);
    assignin('base','spec',CURRSPEC);
    clear tmp_data
    cfg.publish=0;
    fprintf('sim_data assigned to command line at t=%g\n',t);
    msgbox(sprintf('simulated data assigned to command line in variable ''sim_data'' at t=%g',t));
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
        list=get(H.lst_vars(k),'string');
        vind=get(H.lst_vars(k),'value');
        if var_flag(k)==1 % plot state var
          for j=1:length(inds)
            set(H.simdat_alltrace(k,j),'ydata',cfg.record(plotvars{k}(inds(j)),:));
          end
          %set(get(H.ax_state_plot(k),'title'),'string',sprintf('%s (n=%g/%g)',strrep(list{vind},'_','\_'),min(numcell,cfg.ncellshow),numcell));
        elseif var_flag(k)==0 % plot aux function
          % -----------------------------------------------------------
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
cfg.quitflag=1;
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
  %delete(H.ax_state_title(ishandle(H.ax_state_title)));
  delete(H.simdat_LFP_power(ishandle(H.simdat_LFP_power)));
  delete(H.lst_vars(ishandle(H.lst_vars)));
  %delete(H.ax_state_power(ishandle(H.ax_state_power)));
  H=rmfield(H,{'simdat_alltrace','simdat_LFP','ax_state_plot','simdat_LFP_power'});%,'ax_state_title','ax_state_power'});
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
  %H.ax_state_title(i) = title(titlestr);
  %title(strrep(vars{vind},'_','\_'));%[show{i} '.V']);  %ylabel([show{i} '.V']);  
  %H.ax_state_power(i) = subplot('position',[.7 .7+(i-1)*dy .25 -.8*dy],'parent',H.psims); 
  H.simdat_LFP_power(i)=line('color','k','LineStyle','-','erase','background','xdata',cfg.f,'ydata',zeros(size(cfg.f)),'zdata',[],'linewidth',2);
  %H.ax_state_img(i) = subplot('position',[.7 .7+(i-1)*dy .25 -.8*dy],'parent',H.psims); 
  %H.img_state = imagesc(zeros(10,40)); axis xy; axis square
end
% % slider control
if isempty(findobj('tag','speed'))
  uicontrol('Style','frame', 'Units','normalized', ...
            'Position',[0.1  0.05 0.41 0.05],'parent',H.psims,'visible','off');
  uicontrol('Style','text', 'Units','normalized',...
            'Position',[0.15  0.075 0.3 0.02],'parent',H.psims,'string','visualization speed','visible','off');
  uicontrol('Style','slider', 'Units','normalized', ...
            'Position',[0.15  0.055 0.3 0.015],'parent',H.psims,'visible','off',...
            'value',0.0, 'tag','speed');
end
if ~isfield(H,'edit_notes')
  H.edit_notes = uicontrol('style','edit','units','normalized','Position',[.02 .03 0.55 .13],'parent',H.psims,...
    'Callback',@RecordNotes,'string','[record observations here; review in history tab]','BackgroundColor','w','fontsize',10,'HorizontalAlignment','Left','Max',100);
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
if isempty(findobj('tag','publish'))
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.9 0.05 0.075 0.05],... % [0.85  0.11 0.04 0.05]
            'String','get sim_data','Callback','global cfg;cfg.publish=1;');
end
if isempty(findobj('tag','stop'))
  % btn: stop
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.85  0.05 0.04 0.05],...
            'String','stop','tag','stop','Callback','global cfg;cfg.quitflag=1;');
end     
if isempty(findobj('tag','dtlabel'))
  uicontrol('style','text','tag','dtlabel','string','dt','Units','normalized','position',[.75 .14 .04 .02]);
end
if isempty(findobj('tag','dt'))
  % btn: start <=> reset        
  uicontrol('Style','edit','Units','normalized', ...
            'Position',[0.75  0.11 0.04 0.03],...
            'String',num2str(cfg.dt),'tag','dt','Callback','global cfg; cfg.dt=str2num(get(gcbo,''string''));'); % start <=> pause
end        
if isempty(findobj('tag','bufferlabel'))
  uicontrol('style','text','tag','bufferlabel','string','#points','Units','normalized','position',[.8 .14 .04 .02]);
end
if isempty(findobj('tag','buffer'))
  % btn: start <=> reset        
  uicontrol('Style','edit', 'Units','normalized', ...
            'Position',[0.8  0.11 0.04 0.03],...
            'String',num2str(cfg.buffer),'tag','buffer','Callback','global cfg; cfg.buffer=str2num(get(gcbo,''string''));'); % start <=> pause
end     
% autoscale       
uicontrol('Style','pushbutton', 'Units','normalized', ...
          'Position',[0.9 0.11 0.075 0.05],...
          'String','autoscale','Callback',{@setlimits,'autoscale'});

% uicontrol('Style','pushbutton', 'Units','normalized', ...
%           'Position',[0.9 0.05 0.075 0.05],...
%           'String','email report','Callback',@emailreport);      

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
Display_Mech_Info;
DrawUserParams;
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
Display_Mech_Info;
DrawUserParams;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Display_Mech_Info(src,evnt,complabel,mechlabel,type)%mechpos,hmech,connname)
% purpose: display the mech model in a readable form
global H CURRSPEC cfg
if isempty(CURRSPEC.cells), return; end

umech=[]; cnt=1; mechlabels={};
for i=1:length(CURRSPEC.cells)
  for j=1:length(CURRSPEC.cells(i).mechanisms)
    umech(cnt).celllabel = CURRSPEC.cells(i).label;
    umech(cnt).mechlabel = CURRSPEC.cells(i).mechanisms{j};
    umech(cnt).type = 'cells';
    mechlabels{end+1}=sprintf('%s.%s',umech(cnt).celllabel,umech(cnt).mechlabel);
    cnt=cnt+1;
  end
  for j=1:length(CURRSPEC.cells)
    for k=1:length(CURRSPEC.connections(i,j).mechanisms)
      umech(cnt).celllabel = CURRSPEC.connections(i,j).label;
      umech(cnt).mechlabel = CURRSPEC.connections(i,j).mechanisms{k};
      umech(cnt).type = 'connections';
      mechlabels{end+1}=sprintf('%s.%s',umech(cnt).celllabel,umech(cnt).mechlabel);
      cnt=cnt+1;    
    end
  end
end
oldmechs=get(H.lst_mechs,'string');
oldvalue=get(H.lst_mechs,'value');
if numel(oldmechs)>=1 && numel(oldvalue)>=1
  oldmech=oldmechs{oldvalue};
else
  oldmech=[];
end
if ~isequal(oldmechs,mechlabels)
  if ~isempty(oldmech)
    newvalue=find(strcmp(oldmech,mechlabels));
  else
    newvalue=1;
  end
  newmechs=mechlabels;
else
  newvalue=oldvalue;
  newmechs=oldmechs;
end
if isempty(newvalue) && length(newmechs)>1
  newvalue=1;
end
% update mech list
set(H.lst_mechs,'string',newmechs,'value',newvalue,'userdata',umech);

% update mech text
if isempty(umech)
  return;
end
m=umech(newvalue);
if isempty(m), return; end
cfg.focus=find(cellfun(@(x)isequal(m.celllabel,x),{CURRSPEC.(m.type).label}));
cfg.focustype=m.type;
cfg.focusmech=find(strcmp(m.mechlabel,CURRSPEC.(m.type)(cfg.focus).mechanisms));
u.focustype=cfg.focustype;
u.focus=cfg.focus;
u.mechlabel=m.mechlabel;
mech = CURRSPEC.(cfg.focustype)(cfg.focus).mechs(cfg.focusmech);
set(H.txt_mech,'string',mech_spec2str(mech),'userdata',u);

% update mech plot
DrawAuxFunctions;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SaveMech(src,evnt)
global CURRSPEC cfg H allmechs
% purpose: write mech to new file
UpdateMech;
txt = get(H.txt_mech,'string');
if isempty(txt)
  msgbox('empty mechanism. nothing to write.');
  return;
end
% get file name
u=get(H.txt_mech,'userdata');
defaultname=[u.mechlabel '.txt'];
[filename,pathname] = uiputfile({'*.txt;'},'Save as',defaultname);
if isequal(filename,0) || isequal(pathname,0)
  return;
end
outfile = fullfile(pathname,filename);
fid = fopen(outfile,'wt');
for i=1:length(txt)
  fprintf(fid,[txt{i} '\n']);
end
fclose(fid);
fprintf('mechanism written to file: %s\n',outfile);
% update internal structures
thiscell=CURRSPEC.(u.focustype)(u.focus);
mind=strcmp(u.mechlabel,thiscell.mechanisms);
thismech=thiscell.mechs(mind);
thismech.file=outfile;
allmechs(strcmp(u.mechlabel,{allmechs.label}))=thismech;
cfg.allmechfiles{end+1}=outfile;
CURRSPEC.files{end+1}=outfile;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateMech(src,evnt)%,htxt,connname,mechname)
% purpose: apply user changes to the mech model
global CURRSPEC H cfg allmechs
u=get(H.txt_mech,'userdata');
txt = get(H.txt_mech,'string');
newmech = parse_mech_spec(txt);
newmech.label = u.mechlabel;
if ismember(newmech.label,cfg.newmechs)
  ind=find(strcmp(newmech.label,{allmechs.label})); % index into allmechs
  ind2=ismember(cfg.newmechs,newmech.label); % index into cfg.newmechs
  tmp=newmech;
  tmp.file=allmechs(ind).file;
  allmechs(ind)=tmp;
  cfg.newmechs(ind2)=[];
end
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
Display_Mech_Info;
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
      newmechlist = {};
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
  case 'dynamics'
    newspec.cells(this).(field) = get(src,'string'); % splitstr(get(src,'string'),',');
    updatemodel(newspec);
end
Display_Mech_Info;

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
Display_Mech_Info;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CURRSPEC = addmech(CURRSPEC,mechadded,type,index)
% purpose: add a mechanism to the compartment model
global allmechs cfg
if isempty(mechadded), return; end
if ~iscell(mechadded), mechadded={mechadded}; end
% addfile_flag=1;
for i=1:length(mechadded)
  newmech = mechadded{i};
  mechind = find(strcmp({allmechs.label},newmech),1,'last');
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
    fprintf('Creating new mechanism: %s\n',mechadded{i});
    this.params=[];
    this.auxvars={};
    this.functions={};
    this.statevars={};
    this.odes={};
    this.ic={};
    this.substitute={};
    this.inputvars={};
    this.label = mechadded{i};
    this.file = '';%'new';
%     addfile_flag=0;
    cfg.allmechfiles{end+1}=this.file;
    allmechs(end+1)=this;
    newmech=this;
    mechind=length(allmechs);
    cfg.newmechs{end+1}=this.label;
    %warndlg([mechadded{i} ' not found. Check spelling and case.']);
    %disp('known mechanisms include: '); disp(get_mechlist');
  end
%   else
    newmech = rmfield(newmech,'file');
    if ~isempty(CURRSPEC.(type)(index)) && isfield(CURRSPEC.(type)(index),'mechs') && isstruct(CURRSPEC.(type)(index).mechs)
      CURRSPEC.(type)(index).mechs(end+1)=newmech;
    else
      CURRSPEC.(type)(index).mechs = newmech;
    end
    CURRSPEC.(type)(index).mechanisms{end+1}=mechadded{i};
%     if addfile_flag
      CURRSPEC.files{end+1} = cfg.allmechfiles{mechind};
%     end
%   end
end
Display_Mech_Info;
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
Display_Mech_Info;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undo(src,evnt)
% revert to the last working model
global LASTSPEC
updatemodel(LASTSPEC);
refresh;%SelectCells;
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
%   global CURRSPEC
%   updatemodel(CURRSPEC);
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
DrawAuxFunctions;
DrawUserParams;%([],[],[],0);
DrawStudyInfo;
UpdateHistory;
Display_Mech_Info;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printmodel(src,evnt)
global CURRSPEC
buildmodel2(CURRSPEC,'verbose',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function txt = mech_spec2str(mech)
% Purpose: prepare text to display mech model parameters and equations
txt = {}; n=0;
if isempty(mech)
  return;
end
% print parameters
if ~isempty(mech.params)
  keys=fieldnames(mech.params);
  vals=struct2cell(mech.params);
  for i=1:length(keys)
    if i==1, n=n+1; txt{n}=sprintf('%% Parameters:'); end
    n=n+1; txt{n}=sprintf('%s = %s',keys{i},dat2str(vals{i}));
    if i==length(keys), n=n+1; txt{n}=sprintf(' '); end
  end
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
  val = ['[' num2str(val) ']'];
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
function DrawAuxFunctions(src,evnt,limflag)
global H cfg CURRSPEC
if nargin<3, limflag=0; end
if limflag
  l=str2num(get(H.edit_static_lims,'string')); 
  cfg.V=linspace(l(1),l(2),cfg.buffer);
end
% get list of functions in focuscomp
maxlhs = 20; % limit how much is shown in the listbox
maxlen = 150;
funcs={};
for i=1:length(CURRSPEC.cells)
  funcs=cat(1,funcs,CURRSPEC.cells(i).functions);
end
%funcs = CURRSPEC.cells(cfg.focuscomp).functions;
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
functions = funcs(sel,:);%CURRSPEC.cells(cfg.focuscomp).functions(sel,:);
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
if limflag
  set(H.ax_static_plot,'xlim',[min(cfg.V(:)) max(cfg.V(:))]); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawUserParams(src,evnt,jnk,movetabs)
if nargin<4, movetabs=0; end
global CURRSPEC cfg H
if movetabs
  % switch to cellview:
  set(findobj('tag','ptoggle'),'visible','off'); 
  set(findobj('tag','tab'),'backgroundcolor',[1 1 1]); 
  set(findobj('userdata','pcell'),'visible','on'); 
  set(H.bcell,'backgroundcolor',[.7 .7 .7]);
end  

uparm=[]; cnt=1; parmlabels={};
for i=1:length(CURRSPEC.cells)
  for j=1:length(CURRSPEC.cells(i).parameters)/2
    uparm(cnt).celllabel = CURRSPEC.cells(i).label;
    uparm(cnt).parmlabel = CURRSPEC.cells(i).parameters{1+2*(j-1)};
    uparm(cnt).parmvalue = CURRSPEC.cells(i).parameters{2*j};
    uparm(cnt).type = 'cells';
    parmlabels{end+1}=sprintf('%s.%s',uparm(cnt).celllabel,uparm(cnt).parmlabel);
    cnt=cnt+1;
  end
  for j=1:length(CURRSPEC.cells)
    for k=1:length(CURRSPEC.connections(i,j).parameters)/2
      uparm(cnt).celllabel = CURRSPEC.connections(i,j).label;
      uparm(cnt).parmlabel = CURRSPEC.connections(i,j).parameters{1+2*(k-1)};
      uparm(cnt).parmvalue = CURRSPEC.connections(i,j).parameters{2*k};
      uparm(cnt).type = 'connections';
      parmlabels{end+1}=sprintf('%s.%s',uparm(cnt).celllabel,uparm(cnt).parmlabel);
      cnt=cnt+1;    
    end
  end
end

if ~isempty(uparm)
  parmind=1;
  parmval=num2str(uparm(parmind).parmvalue);
else
  parmind=[];
  parmval='';
  parmlabels='';
end
if ~isfield(H,'lst_parms') || ~ishandle(H.lst_parms)
  H.lst_parms = uicontrol('parent',H.pcell,'units','normalized','BackgroundColor',[.9 .9 .9],'Max',1,'TooltipString','Right-click to edit parameter name. Hit ''d'' to delete.',...
    'position',[0 0 .2 1],'style','listbox','value',parmind,'string',parmlabels,...
    'ButtonDownFcn',@RenameParm,'Callback',{@UpdateParams,'show'},'KeyPressFcn',{@UpdateParams,'delete'},'userdata',uparm);
  uicontrol('style','text','parent',H.pcell,'units','normalized','string','add/update parameter:',...
    'position',[.3 .85 .6 .05],'HorizontalAlign','Left');
  H.edit_parmadd = uicontrol('parent',H.pcell,'units','normalized','style','edit','TooltipString','format: source-target.param = value',...
    'position',[.3 .8 .6 .05],'backgroundcolor','w','string','name = value',...
    'HorizontalAlignment','left','Callback',{@UpdateParams,'add'});
%   uicontrol('style','text','parent',H.pcell,'units','normalized','string','update value:',...
%     'position',[.3 .7 .6 .05],'HorizontalAlign','Left');
  H.edit_parmedit = uicontrol('parent',H.pcell,'units','normalized','style','edit','TooltipString','select parameter in list and enter new value here',...
    'position',[0 0 .2 .05],'backgroundcolor','w','string',parmval,...
    'HorizontalAlignment','left','Callback',{@UpdateParams,'change'},'visible','off');
else
  set(H.lst_parms,'value',parmind,'string',parmlabels,'userdata',uparm);
  set(H.edit_parmedit,'string',parmval);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateParams(src,evnt,action)
global H CURRSPEC
uparm=get(H.lst_parms,'userdata');
sel=get(H.lst_parms,'value');
switch action
  case 'show'
    s=get(H.lst_parms,'string');
    val=num2str(uparm(sel).parmvalue);
    set(H.edit_parmadd,'string',[s{sel} ' = ' val]);
    set(H.edit_parmedit,'string',val);
  case 'change'
    newspec=CURRSPEC;
    u=uparm(sel);
    cellind=cellfun(@(x)isequal(x,u.celllabel),{newspec.(u.type).label});
    keys=newspec.(u.type)(cellind).parameters(1:2:end);
    parmind = 2*find(strcmp(u.parmlabel,keys));
    newspec.(u.type)(cellind).parameters{parmind}=str2double(get(H.edit_parmedit,'string'));
    updatemodel(newspec);
    DrawUserParams;
  case 'add' % format: source-target.param = value
    newspec=CURRSPEC;
    str = get(H.edit_parmadd,'string');
    if isempty(str), return; end
    parts=splitstr(str,'=');
    value = str2double(strtrim(parts{2}));
    parts=splitstr(parts{1},'.');
    if numel(parts)==1 && length(CURRSPEC.cells)==1
      parts={CURRSPEC.cells(1).label,parts{1}};
    elseif numel(parts)==1
      parts={'__all__',parts{1}};
      %msgbox('improper syntax. use: compartment.parameter = value','syntax error','error');
      %return;
    end
    param= strtrim(parts{2});
    tmpparts=splitstr(parts{1},'-');
    if numel(tmpparts)==1
      type='cells';
      if strcmp(parts{1},'__all__')
        targets=1:length(CURRSPEC.cells);
      else
        target=find(strcmp(strtrim(parts{1}),{CURRSPEC.cells.label}));
        targets=target;
      end
    elseif numel(parts)==2
      type='connections';
      target=find(cellfun(@(x)isequal(strtrim(parts{1}),x),{CURRSPEC.connections.label}));
      if isempty(target)
        src=tmpparts{1}; srcind=find(strcmp(strtrim(src),{CURRSPEC.cells.label}));
        dst=tmpparts{2}; dstind=find(strcmp(strtrim(dst),{CURRSPEC.cells.label}));
        target=sub2ind(size(CURRSPEC.connections),srcind,dstind);
      end
      targets=target;
    end
    for i=1:length(targets)
      target=targets(i);
      tmp=newspec.(type)(target).parameters;
      if ~iscell(tmp)
        tmp={};
      end
      if ~isempty(tmp) && any(ismember(param,tmp(1:2:end))) % param already exists. just update it.
        ind=2*find(strcmp(param,tmp(1:2:end)));
        newspec.(type)(target).parameters{ind}=value;
      else % new parameter. add it.
        newspec.(type)(target).parameters{end+1}=param;
        newspec.(type)(target).parameters{end+1}=value;
      end
    end
    updatemodel(newspec);
    DrawUserParams;
  case 'delete' % on KeyPress 'd' or 'delete'
    if strcmp(evnt.Key,'d') || strcmp(evnt.Key,'delete')
      newspec=CURRSPEC;
      u=uparm(sel);
      cellind=cellfun(@(x)isequal(x,u.celllabel),{newspec.(u.type).label});
      keys=newspec.(u.type)(cellind).parameters(1:2:end);
      parmind = 2*find(strcmp(u.parmlabel,keys))-1;
      % remove parameter
      if ~isempty(parmind)
        newspec.(u.type)(cellind).parameters([parmind parmind+1])=[];
      else
        return;
      end
      updatemodel(newspec);
      DrawUserParams;      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DrawStudyInfo
global cfg H
if isfield(cfg,'study')
  study = cfg.study;
else
  study.scope = '';
  study.variable = '';
  study.values = '';  
  cfg.study = study;
end
% H.pbatchcontrols = pbatchcontrols;
% H.pbatchspace = pbatchspace;
% H.pbatchoutputs = pbatchoutputs;
if ~isfield(H,'text_scope') || ~ishandle(H.text_scope)
  % controls
  yshift=-.05; ht=.17;
  uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','text','position',[.05 .75+yshift .13 .2],'string','machine',...
    'HorizontalAlignment','left');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','text','position',[.5 .75+yshift .13 .2],'string','memlimit',...
    'HorizontalAlignment','right');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','text','position',[.05 .5+yshift .13 .2],'string','timelimits',...
    'HorizontalAlignment','left');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','text','position',[.05 .3+yshift .13 .2],'string','solver',...
    'HorizontalAlignment','left');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','text','position',[.3 .3+yshift .13 .2],'string','dt',...
    'HorizontalAlignment','right');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','text','position',[.05 .1+yshift .13 .2],'string','#repeats',...
    'HorizontalAlignment','left');%,'backgroundcolor','w'
  H.rad_machine=uibuttongroup('visible','off','units','normalized','Position',[.18 .8+yshift .3 .2],'parent',H.pbatchcontrols);
  H.rad_machine_1=uicontrol('style','radiobutton','string','local','parent',H.rad_machine,'HandleVisibility','off',...
    'units','normalized','pos',[0 0 .4 1]);
  H.rad_machine_2=uicontrol('style','radiobutton','string','cluster','parent',H.rad_machine,'HandleVisibility','off',...
    'units','normalized','pos',[.45 0 .5 1]);
  set(H.rad_machine,'SelectedObject',H.rad_machine_1);  % No selection
  set(H.rad_machine,'Visible','on');    
  H.edit_memlimit = uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','edit','position',[.65 .8+yshift .1 ht],'backgroundcolor','w','string','8G',...
    'HorizontalAlignment','left');
  H.edit_timelimits = uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','edit','position',[.18 .55+yshift .2 ht],'backgroundcolor','w','string','[0 40]',...
    'HorizontalAlignment','left');
  H.edit_SOLVER = uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','edit','position',[.18 .35+yshift .2 ht],'backgroundcolor','w','string','euler',...
    'HorizontalAlignment','left');
  H.edit_dt = uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','edit','position',[.45 .35+yshift .1 ht],'backgroundcolor','w','string','0.01',...
    'HorizontalAlignment','left');
  H.edit_repeats = uicontrol('parent',H.pbatchcontrols,'units','normalized',...
    'style','edit','position',[.18 .15+yshift .2 ht],'backgroundcolor','w','string','1',...
    'HorizontalAlignment','left');  
  % search space
  H.text_scope = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','text','position',[.1 .9 .1 .05],'string','scope',...
    'HorizontalAlignment','center');%,'backgroundcolor','w'
  H.text_variable = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','text','position',[.31 .9 .1 .05],'string','variable',...
    'HorizontalAlignment','center');%,'backgroundcolor','w'
  H.text_values = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','text','position',[.52 .9 .1 .05],'string','values',...
    'HorizontalAlignment','center'); %,'backgroundcolor','w'
  H.btn_batch_help = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','pushbutton','fontsize',10,'string','help','callback','web(''https://github.com/jsherfey/research/blob/master/modeling/biosim/readme'');',...
    'position',[.85 .92 .1 .06]); 
  % outputs
  uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','text','position',[.05 .85 .13 .1],'string','rootdir',...
    'HorizontalAlignment','left');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','text','position',[.8 .85 .1 .1],'string','dsfact',...
    'HorizontalAlignment','right');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','text','position',[.05 .7 .1 .1],'string','save',...
    'HorizontalAlignment','right');%,'backgroundcolor','w'
  uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','text','position',[.3 .7 .1 .1],'string','plot',...
    'HorizontalAlignment','right');%,'backgroundcolor','w'  
  H.edit_rootdir = uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','edit','position',[.15 .85 .65 .15],'backgroundcolor','w','string',pwd,...
    'HorizontalAlignment','left');
  H.edit_dsfact = uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','edit','position',[.92 .85 .05 .15],'backgroundcolor','w','string','1',...
    'HorizontalAlignment','left');  
  H.btn_run_simstudy = uicontrol('parent',H.pbatchoutputs,'units','normalized',...
    'style','pushbutton','fontsize',20,'string','submit!','callback',@RunSimStudy,...
    'position',[.67 .27 .3 .4]); 
  H.chk_savedata=uicontrol('style','checkbox','value',0,'parent',H.pbatchoutputs,'units','normalized','position'   ,[.13 .6 .15 .1],'string','data');
  H.chk_savesum=uicontrol('style','checkbox','value',0,'parent',H.pbatchoutputs,'units','normalized','position'    ,[.13 .5 .15 .1],'string','popavg');
  H.chk_savespikes=uicontrol('style','checkbox','value',0,'parent',H.pbatchoutputs,'units','normalized','position' ,[.13 .4 .15 .1],'string','spikes');
  H.chk_saveplots=uicontrol('style','checkbox','value',0,'parent',H.pbatchoutputs,'units','normalized','position'  ,[.13 .3 .15 .1],'string','plots');
  H.chk_plottraces=uicontrol('style','checkbox','value',1,'parent',H.pbatchoutputs,'units','normalized','position' ,[.38 .6 .15 .1],'string','state vars');
  H.chk_plotrates=uicontrol('style','checkbox','value',0,'parent',H.pbatchoutputs,'units','normalized','position'  ,[.38 .5 .17 .1],'string','spike rates');
  H.chk_plotspectra=uicontrol('style','checkbox','value',0,'parent',H.pbatchoutputs,'units','normalized','position',[.38 .4 .15 .1],'string','spectrum');
end
if isfield(H,'edit_scope') 
  if ishandle(H.edit_scope)
    delete(H.edit_scope);
    delete(H.edit_variable);
    delete(H.edit_values);
    delete(H.btn_simset_delete);
    delete(H.btn_simset_copy);
  end
  H = rmfield(H,{'edit_scope','edit_variable','edit_values','btn_simset_delete','btn_simset_copy'});
end 
for i=1:length(study)
  H.edit_scope(i) = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','edit','position',[.1 .8-.1*(i-1) .2 .08],'backgroundcolor','w','string',study(i).scope,...
    'HorizontalAlignment','left','Callback',sprintf('global cfg; cfg.study(%g).scope=get(gcbo,''string'');',i));
  H.edit_variable(i) = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','edit','position',[.31 .8-.1*(i-1) .2 .08],'backgroundcolor','w','string',study(i).variable,...
    'HorizontalAlignment','left','Callback',sprintf('global cfg; cfg.study(%g).variable=get(gcbo,''string'');',i));
  H.edit_values(i) = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','edit','position',[.52 .8-.1*(i-1) .4 .08],'backgroundcolor','w','string',study(i).values,...
    'HorizontalAlignment','left','Callback',sprintf('global cfg; cfg.study(%g).values=get(gcbo,''string'');',i)); 
  H.btn_simset_delete(i) = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','pushbutton','fontsize',10,'string','-','callback',{@DeleteSimSet,i},...
    'position',[.06 .8-.1*(i-1) .03 .08]);
  H.btn_simset_copy(i) = uicontrol('parent',H.pbatchspace,'units','normalized',...
    'style','pushbutton','fontsize',10,'string','+','callback',{@CopySimSet,i},...
    'position',[.93 .8-.1*(i-1) .03 .08]);
end

function DeleteSimSet(src,evnt,index)
global cfg
cfg.study(index) = [];
DrawStudyInfo;

function CopySimSet(src,evnt,index)
global cfg
cfg.study(end+1) = cfg.study(index);
DrawStudyInfo;

function RunSimStudy(src,evnt)
global cfg CURRSPEC H
if isempty(CURRSPEC.cells) || isempty([cfg.study.scope]), return; end
scope = {cfg.study.scope};
variable = {cfg.study.variable};
values = {cfg.study.values};
dir=get(H.edit_rootdir,'string');
mem=get(H.edit_memlimit,'string');
dt=str2num(get(H.edit_dt,'string'));
lims=str2num(get(H.edit_timelimits,'string'));
dsfact=str2num(get(H.edit_dsfact,'string'));

machine=get(get(H.rad_machine,'SelectedObject'),'String');
if strcmp(machine,'local')
  clusterflag = 0;
elseif strcmp(machine,'cluster')
  clusterflag = 1;
end
nrepeats=str2num(get(H.edit_repeats,'string'));
for i=1:nrepeats
  if nrepeats>1
    fprintf('Submitting batch iteration %g of %g...\n',i,nrepeats);
  end
  % record note
  if ~isfield(CURRSPEC,'history') || isempty(CURRSPEC.history)
    id=1; 
  else
    id=max([CURRSPEC.history.id])+1;
  end
  note.id=id;
  timestamp=datestr(now,'yyyymmdd-HHMMSS');
  note.date=timestamp;
  if strcmp(machine,'local')
    note.text=sprintf('BATCH: machine=%s. ',machine);
  else
    note.text=sprintf('BATCH: machine=%s. ',machine);% rootdir=%s',machine,dir);
  end
  tmp=CURRSPEC;
  if isfield(tmp.model,'eval')
    tmp.model=rmfield(tmp.model,'eval');
  end
  if isfield(tmp,'history')
    tmp = rmfield(tmp,'history');
  end
  if id>1 && isequal(CURRSPEC.history(end).spec.model,CURRSPEC.model)
    note.id=id-1;
    note.spec=tmp;
    note.changes={};
  else
    note.spec=tmp;
    note.changes={'changes made'};
  end
  note.isbatch=1;  
  note.batch.space=cfg.study;
  note.batch.rootdir = dir;
  note.batch.machine = machine;
  % update controls
  s=get(H.lst_notes,'string');
  v=get(H.lst_notes,'value');
  s={s{:} num2str(note.id)};
  v=[v length(s)];
  set(H.lst_notes,'string',s,'value',v);
  set(H.edit_notes,'string','');
  if id==1
    CURRSPEC.history = note;
  else
    CURRSPEC.history(end+1) = note;
  end
  UpdateHistory;
  % submit to simstudy
  [allspecs,timestamp]=simstudy(CURRSPEC,scope,variable,values,'dt',dt,'rootdir',dir,'memlimit',mem,...
    'timelimits',lims,'dsfact',dsfact,'sim_cluster_flag',clusterflag,'timestamp',timestamp,...
    'savedata_flag',get(H.chk_savedata,'value'),'savepopavg_flag',get(H.chk_savesum,'value'),'savespikes_flag',get(H.chk_savespikes,'value'),'saveplot_flag',get(H.chk_saveplots,'value'),...
    'plotvars_flag',get(H.chk_plottraces,'value'),'plotrates_flag',get(H.chk_plotrates,'value'),'plotpower_flag',get(H.chk_plotspectra,'value'));
end

% notes(2).id='model2';
% notes(2).date='yyyymmdd-hhmmss';
% notes(2).text='batch sim note...';
% notes(2).changes{1}='E.n: 1 => 2';
% notes(2).changes{2}='+E2: {iNa,iK}';
% notes(2).isbatch = 1;
% notes(2).batch.space(1).scope = '(E,I)';
% notes(2).batch.space(1).variables = 'N';
% notes(2).batch.space(1).values = '[1 2 3]';
% notes(2).model = CURRSPEC;
% CURRSPEC.history(end+1)=notes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateHistory(src,evnt)
global CURRSPEC H

if ~isfield(H,'lst_notes') || ~ishandle(H.lst_notes)
  notes=[]; ids={};
  H.lst_notes = uicontrol('units','normalized','position',[0 0 .2 1],'parent',H.phistory,'BackgroundColor',[.9 .9 .9],'ToolTipString','Select notes and hit ''d'' to delete.',...
    'style','listbox','value',1:min(3,length(notes)),'string',ids,'Max',100,'Callback',@UpdateHistory,'KeyPressFcn',@NoteKeyPress);
  H.edit_comparison = uicontrol('parent',H.pcomparison,'style','edit','units','normalized','tag','modelcomparison',...
  'position',[0 0 1 .85],'string','','FontName','Courier','FontSize',9,'HorizontalAlignment','Left','Max',100,'BackgroundColor',[.9 .9 .9]);
  jEdit = findjobj(H.edit_comparison);
  jEditbox = jEdit.getViewport().getComponent(0);
  jEditbox.setEditable(false);                % non-editable  
  H.btn_compare = uicontrol('parent',H.pcomparison,'units','normalized',...
    'style','pushbutton','fontsize',10,'string','compare','callback',@CompareModels,...
    'position',[.05 .85 .2 .15]);  
  H.btn_report = uicontrol('parent',H.pcomparison,'style','pushbutton','fontsize',10,'string','gen report','callback',@GenerateReport,...
    'units','normalized','position',[.75 .85 .2 .15]);%[.25 .935 .1 .03]);
end

if ~isfield(CURRSPEC,'history') || isempty(CURRSPEC.history)
  return;
else
  notes=CURRSPEC.history;
end

ids = cellfun(@num2str,{notes.id},'uni',0);
if isempty(ids)
  return;
end
str = get(H.lst_notes,'string');
if ~isequal(ids,str)
  set(H.lst_notes,'string',ids);
end
sel = get(H.lst_notes,'value'); 
sel = sel(1:min(length(sel),length(ids)));
if numel(sel)<1
  sel = 1;
elseif any(sel)>length(notes)
  sel(sel>length(notes))=[];
end
set(H.lst_notes,'value',sel);
notes = notes(sel);

delete(findobj('tag','note'));

ypos=.9; ht=.05;
for i=1:length(notes)
  if notes(i).isbatch==1 && strcmp(notes(i).batch.machine,'cluster')
    saved_flag=1; fontcolor='k';
  else
    saved_flag=0; fontcolor='k';
  end
  if length(notes(i).changes)>=1
    changed_flag=1; fontweight='bold';
  else
    changed_flag=0; fontweight='normal';
  end
  if changed_flag
    H.chk_notes(i) = uicontrol('style','checkbox','value',0,'parent',H.pnotes,'units','normalized','userdata',notes(i),...
      'position',[.05 ypos .7 ht],'string',sprintf('Model %s (%s)',ids{sel(i)},notes(i).date),'visible','on','foregroundcolor',fontcolor,'tag','note','fontweight',fontweight); % ['Model ' ids{sel(i)} ' (' notes(i).date ')']
  else
    H.chk_notes(i) = uicontrol('style','checkbox','value',0,'parent',H.pnotes,'units','normalized','userdata',notes(i),...
      'position',[.05 ypos .7 ht],'string',sprintf('Model %s (%s)',ids{sel(i)},notes(i).date),'visible','on','foregroundcolor',fontcolor,'tag','note','fontweight',fontweight); % ['Model ' ids{sel(i)} ' (' notes(i).date ')']    
  end
  H.btn_revert(i) = uicontrol('parent',H.pnotes,'units','normalized','style','pushbutton','fontsize',10,'visible','on','tag','note',...
    'string','<--','position',[.87 ypos .1 .04],'callback',{@updatemodel,notes(i).spec,notes(i).id});%%sprintf('global CURRSPEC; updatemodel(CURRSPEC.history(%g).spec); refresh;',sel(i)));
  ypos=ypos-ht;
  H.edit_note(i) = uicontrol('style','edit','units','normalized','HorizontalAlignment','left','parent',H.pnotes,'BackgroundColor','w','visible','on',...
    'string',notes(i).text,'position',[.15 ypos .82 ht],'tag','note','callback',sprintf('global CURRSPEC; CURRSPEC.history(%g).text=get(gcbo,''string'');',sel(i)));
  if notes(i).isbatch==0
    ypos=ypos-ht;
  else
    if saved_flag
      H.btn_batchmanager(i) = uicontrol('parent',H.pnotes,'units','normalized','style','pushbutton','fontsize',10,...
        'string','open','tag','note','position',[.87 ypos-ht+.01 .1 .04],'callback',{@Load_File,notes(i).batch.rootdir},'visible','on');  
    end
    for j=1:length(notes(i).batch.space)
      b=notes(i).batch.space(j);
      ypos=ypos-ht;
      H.txt_searchspace = uicontrol('style','text','string',sprintf('(%s).(%s)=%s',b.scope,b.variable,b.values),'fontsize',10,'fontweight','bold',...
        'parent',H.pnotes,'tag','note','units','normalized','position',[.15 ypos .72 ht],'HorizontalAlignment','left','visible','on','foregroundcolor','b');      
    end
    ypos=ypos-ht/1.1;% 1.5
  end
end

% notes(1).id='model1';
% notes(1).date='yyyymmdd-hhmmss';
% notes(1).text='interactive sim note...';
% notes(1).changes{1}='E.N: 1 => 2';
% notes(1).changes{2}='+E2: {iNa,iK}';
% notes(1).changes{3}='E2: +iM';
% notes(1).changes{4}='E2.iM.gM: 10 => 20';
% notes(1).isbatch = 0;
% notes(1).batch = [];
% notes(1).model = CURRSPEC;
% 
% notes(2).id='model2';
% notes(2).date='yyyymmdd-hhmmss';
% notes(2).text='batch sim note...';
% notes(2).changes{1}='E.n: 1 => 2';
% notes(2).changes{2}='+E2: {iNa,iK}';
% notes(2).isbatch = 1;
% notes(2).batch.space(1).scope = '(E,I)';
% notes(2).batch.space(1).variables = 'N';
% notes(2).batch.space(1).values = '[1 2 3]';
% notes(2).model = CURRSPEC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NoteKeyPress(src,evnt)
switch evnt.Key
  case {'delete','d'}
    global CURRSPEC;
    newspec=CURRSPEC;
    s=get(src,'string');
    v=get(src,'value');
    newspec.history(v)=[];
    s(v)=[];
    set(src,'string',s,'value',[]);
    updatemodel(newspec);
    UpdateHistory;
    %refresh;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GenerateReport(src,evnt)
global CURRSPEC cfg
notes=CURRSPEC.history;
txt={};
ids=[notes.id];
uids=unique(ids);
lastmodel=[];
for i=1:length(uids) % loop over models
  id=uids(i);
  these=notes([notes.id]==id);
  thismodel=these(1).spec;
  if ~isempty(lastmodel)
    txt{end+1}=sprintf('diff(Model%g,Model%g): changes in Model%g compared to Model%g',uids(i-1),id,id,uids(i-1));
    txt=cat(2,txt,modeldiff(lastmodel,thismodel));
    txt{end+1}='';
  end  
  txt{end+1}=sprintf('Model %g notes',id);
  for j=1:length(these)
    note=these(j);
    txt{end+1}=sprintf('%s: %s',note.date,note.text);
    if note.isbatch && ~isempty(note.batch.space)
      for k=1:length(note.batch.space)
        b=note.batch.space(k);
        txt{end+1}=sprintf('\t(%s).(%s)=%s',b.scope,b.variable,b.values);
      end    
    end
  end
  if i < length(uids)
    txt{end+1}='-----------------------------------------------------------------';    
  else
    txt{end+1}='';
  end
  lastmodel=thismodel;
  end
txt{end+1}=cfg.modeltext;

h=figure('position',[70 120 930 580]);
uicontrol('parent',h,'style','edit','units','normalized','position',[0 0 1 .85],'tag','report',...
  'string',txt,'FontName','Courier','FontSize',9,'HorizontalAlignment','Left','Max',100,'BackgroundColor','w'); % Courier Monospaced
uicontrol('parent',h,'style','pushbutton','fontsize',10,'string','write','callback',@WriteReport,...
    'units','normalized','position',[0 .9 .1 .05]);
uicontrol('parent',h,'style','pushbutton','fontsize',10,'string','email','callback',@EmailReport,...
    'units','normalized','position',[.1 .9 .1 .05]);

function WriteReport(src,evnt)
timestamp=datestr(now,'yyyymmdd-HHMMSS');
[filename,pathname] = uiputfile({'*.txt;'},'Save as',['report_' timestamp '.txt']);
if isequal(filename,0) || isequal(pathname,0)
  return;
end
txt=get(findobj('tag','report'),'string');
outfile = fullfile(pathname,filename);
fid = fopen(outfile,'wt');
for i=1:length(txt)
  fprintf(fid,[txt{i} '\n']);
end
fclose(fid);
fprintf('report written to file: %s\n',outfile);

function EmailReport(src,evnt)
global CURRSPEC
% prepare report text
txt=get(findobj('tag','report'),'string');
if isempty(txt), return; end
if iscell(txt{1}), txt=txt{end}; end
str='';
for i=1:length(txt)
  str=sprintf('%s%s\n',str,txt{i});
end
% get email address
emailaddress=inputdlg('Email address:','Enter email address');
if isempty(emailaddress)
  return; 
else
  emailaddress=emailaddress{1};
end
% add system info to email content
sysinfo='';
[r, username] = system('echo $USER');
if ~r
  username=regexp(username,'\w+','match'); % remove new line characters
  if ~isempty(username), username=username{1}; end
  sysinfo=sprintf('%sUser: %s\n',sysinfo,username); 
end
[r, computername] = system('hostname');           % Uses linux system command to get the machine name of host. 
if ~r, sysinfo=sprintf('%sHost: %s',sysinfo,computername); end
sysinfo=sprintf('%sMatlab version: %s\n',sysinfo,version);
if exist('BIOSIMROOT','var') % get unique id for this version of the git repo
  cwd=pwd;
  cd(BIOSIMROOT);
  [r,m]=system('git log -1');
  cd(cwd);
else
  [r,m]=system('git log -1');
end
if ~r
  githash=strrep(regexp(m,'commit\s[\w\d]+','match','once'),'commit ','');
  sysinfo=sprintf('%sDNSim Git hash: %s',sysinfo,githash);
end
str=sprintf('Report generated at %s.\n%s\n\n%s',datestr(now,31),sysinfo,str);
timestamp=datestr(now,'yyyymmdd-HHMMSS');
% create temporary file for attachment
tmpfile=sprintf('model_%s.mat',timestamp);
spec=CURRSPEC;
save(tmpfile,'spec');
% send email
setpref('Internet','SMTP_Server','127.0.0.1'); % Sets the outgoing mail server - often the default 127.0.0.1
setpref('Internet','E_mail','sherfey@bu.edu'); % Sets the email FROM/reply address for all outgoing email reports.
sendmail(emailaddress,['DNSim Report: ' timestamp],str,tmpfile);
fprintf('report emailed to: %s\n',emailaddress);
% remove temporary file
delete(tmpfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CompareModels(src,evnt)
global H
if ~isfield(H,'chk_notes'), return; end
H.chk_notes(~ishandle(H.chk_notes))=[];
if isempty(H.chk_notes), return; end
sel=get(H.chk_notes,'value');
if isempty(sel),return; end
sel=[sel{:}]==1;
if length(find(sel))<2, return; end
notes = get(H.chk_notes(sel),'userdata');
basemodel=notes{1}.spec;
othermodels=cellfun(@(x)x.spec,notes(2:end),'uni',0);
txt={};
for i=1:length(othermodels)
  txt{end+1}=sprintf('diff(Model%g,Model%g): changes in Model%g compared to Model%g',notes{1}.id,notes{i+1}.id,notes{i+1}.id,notes{1}.id);
  txt=cat(2,txt,modeldiff(basemodel,othermodels{i}));
  txt{end+1}='------------------------------------------------';
end
set(H.edit_comparison,'string',txt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RenameComponent(src,evnt)
global CURRSPEC
s=get(src,'string');
v=get(src,'value');
if length(v)>1, return; end
newname=inputdlg(['Rename Compartment: ' s{v}],'New name');
if isempty(newname), return; end
newname=newname{1};
newspec=CURRSPEC;
I=find(strcmp(s{v},{newspec.cells.label}));
newspec.cells(I).label=newname;
for i=1:length(s)
  l=newspec.connections(v,i).label;
  if ~isempty(l)
    newspec.connections(v,i).label=strrep(l,[s{v} '-'],[newname '-']);
  end
  l=newspec.connections(i,v).label;
  if ~isempty(l)
    newspec.connections(i,v).label=strrep(l,['-' s{v}],['-' newname]);
  end
end
s{v}=newname;
set(src,'string',s);
updatemodel(newspec);
SelectCells;
Display_Mech_Info;
DrawAuxView;
DrawUserParams;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RenameMech(src,evnt)
global CURRSPEC
ud=get(src,'userdata');
v=get(src,'value');
s=get(src,'string');
if length(v)>1, return; end
newname=inputdlg(['Rename Mechanism: ' ud(v).mechlabel ' (in ' ud(v).celllabel ')'],'New name');
if isempty(newname), return; end
newname=newname{1};
newspec=CURRSPEC;
u=ud(v);
ud(v).mechlabel=newname;
comp=u.celllabel;
mech=u.mechlabel;
i=find(cellfun(@(x)isequal(comp,x),{newspec.(u.type).label}));
j=find(strcmp(mech,newspec.(u.type)(i).mechanisms));
newspec.(u.type)(i).mechanisms{j}=newname;
newspec.(u.type)(i).mechs(j).label=newname;
s{v}=[comp '.' newname];
set(src,'string',s,'userdata',ud);
updatemodel(newspec);
SelectCells;
Display_Mech_Info;
DrawAuxView;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function RecordNotes(src,evnt)
global CURRSPEC H cfg t
if ~isfield(CURRSPEC,'history') || isempty(CURRSPEC.history)
  id=1; 
else
  id=max([CURRSPEC.history.id])+1;
end
note.id=id;
note.date=datestr(now,'yyyymmdd-HHMMSS');
note.text=get(H.edit_notes,'string');
if cfg.quitflag<0
  note.text = sprintf('SIM(t=%g): %s',t,note.text);
end
tmp=CURRSPEC;
if isfield(tmp.model,'eval')
  tmp.model=rmfield(tmp.model,'eval');
end
if isfield(tmp,'history')
  tmp = rmfield(tmp,'history');
end
if id>1 && isequal(CURRSPEC.history(end).spec.cells,CURRSPEC.cells) && isequal(CURRSPEC.history(end).spec.connections,CURRSPEC.connections) % isequal(CURRSPEC.history(end).spec.model,CURRSPEC.model)
  note.id=id-1;
  note.spec=tmp;
  note.changes={};
else
  note.spec=tmp;
  note.changes={'changes made'};
end
note.isbatch=0;
note.batch=[];
if id==1
  CURRSPEC.history = note;
else
  CURRSPEC.history(end+1) = note;
end
s=get(H.lst_notes,'string');
v=get(H.lst_notes,'value');
s={s{:} num2str(note.id)};
v=[v length(s)];
set(H.lst_notes,'string',s,'value',v);
set(H.edit_notes,'string','');
UpdateHistory;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatemodel(newspec,varargin) % maybe use same specification as for "override"
% purpose: update the integratable model after each change to its specification
if nargin>2, newspec=varargin{2}; end % get spec from history:revert
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
if nargin>2
  CURRSPEC.history = LASTSPEC.history; % hold onto history since the reverted model
  refresh; 
  if nargin>3
    fprintf('reverting to model %g\n',varargin{3});
  else
    fprintf('reverting to selected model\n');
  end
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
