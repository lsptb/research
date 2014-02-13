function cellmodeler(spec)
% Purpose: GUI for building multi-compartment biophysical cell models
% Calls the following external functions:
%   get_mechlist(): get list of known biophysical mechanisms from database
%   buildmodel2(): build ODE system from high-level model specification
%   parse_mech_spec(): read model equations from ascii text file
%   findjobj(): see Matlab File Exchange (http://www.mathworks.com/matlabcentral/fileexchange/14317-findjobj-find-java-handles-of-matlab-graphic-objects)
% Refs:
%   - http://undocumentedmatlab.com/blog/undocumented-mouse-pointer-functions/
%   - http://undocumentedmatlab.com/blog/findjobj-find-underlying-java-object/
%   - http://undocumentedmatlab.com/blog/uicontrol-callbacks/

% TODO (after finishing GUI dev steps 1-4):
% 5. add static plots (auxfuncs)
% 6. add representative sim plots (10ms; state vars and interface funcs)
% 7. add real-time simulation w/ sim controls (start,stop,pause)
% 8. add model load/save capability

global cfg H
if ~isfield(spec,'cells'), tmp.cells=spec; spec=tmp; clear tmp; end
if isfield(spec.cells,'files'), spec.files = spec.cells(1).files; end

% get list of all known mechs (stored in DB)
% TODO: read list fom MySQL DB (see http://introdeebee.wordpress.com/2013/02/22/connecting-matlab-to-mysql-database-using-odbc-open-database-connector-for-windows-7/)
DBPATH = '/space/mdeh3/9/halgdev/projects/jsherfey/code/modeler/database';
if ~exist(DBPATH,'dir')
  DBPATH = 'C:\Users\jsherfey\Desktop\My World\Code\modelers\database';
end
[allmechlist,allmechfiles]=get_mechlist(DBPATH);
% use stored mechs if user did not provide list of mech files
if ~isfield(spec,'files') || isempty(spec.files)
  selmechlist=allmechlist;
  selmechfiles=allmechfiles;
elseif ischar(spec.files) && exist(spec.files,'dir')
  % user provided a directory that contains the user mech files
  d=dir(spec.files);
  selmechlist = {d(cellfun(@(x)any(regexp(x,'.txt$')),{d.name})).name};
  allmechlist = {selmechlist{:} allmechlist{:}}; % prepend user mechs
  selmechfiles = cellfun(@(x)fullfile(DBPATH,x),spec.files,'unif',0);
  allmechfiles = {selmechfiles{:} allmechfiles{:}};
end

% only keep mech selection matching those included in the cell
cellmechs={spec.cells.mechanisms}; cellmechs=[cellmechs{:}];
if isfield(spec,'connections')
  connmechs={spec.connections.mechanisms}; 
  connmechs=connmechs(~cellfun(@isempty,connmechs));
else
  connmechs={};
end
cellmechs={cellmechs{:} connmechs{:}};
selmechlist = cellfun(@(x)strrep(x,'.txt',''),selmechlist,'unif',0);
sel=cellfun(@(x)any(strmatch(x,cellmechs,'exact')),selmechlist);
selmechlist = selmechlist(sel);
selmechfiles = selmechfiles(sel);
spec.files = selmechfiles;

%cfg = [];
for i=1:length(spec.cells), cfg.userparams{i} = spec.cells(i).parameters; end
cfg.selmechlist = selmechlist;
cfg.allmechlist = allmechlist;
cfg.allmechfiles = allmechfiles;
cfg.focuscolor = [.7 .7 .7];
cfg.focuscomp = 1; % index to component-of-focus in spec.cells (start w/ root)
cfg.focusmech = min(1,length(selmechlist));
cfg.pauseflag = -1;
cfg.quitflag = -1;
cfg.tlast=-inf; 
cfg.buffer = 10000;
cfg.dt = .01;

[model,IC,functions,auxvars,spec] = buildmodel2(spec);
try spec.cells=spec.entities; spec=rmfield(spec,'entities'); end

global currspec lastspec
currspec = spec;
lastspec = spec;

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

%% set up figure
%H = [];
if isfield(currspec.cells,'parent') && ischar(currspec.cells(1).parent)
  cname=currspec.cells(1).parent;
else
  cname='cell';
end
if exist('H','var') && isfield(H,'f_cell') && ishandle(H.f_cell)
  close(H.f_cell);
end
sz = get(0,'ScreenSize'); 
H.f_cell = figure('position',[.005*sz(3) .01*sz(4) .94*sz(3) .88*sz(4)],...%'position',[125 100 1400 800],...
  'WindowScrollWheelFcn',@ZoomFunction,'CloseRequestFcn','delete(gcf);');
% H.f_cell = figure('position',[125 150 1400 800],...
%   'WindowScrollWheelFcn',@ZoomFunction,'CloseRequestFcn','delete(gcf);');%'global currspec; currspec=[]; delete(gcf);');
jobj=findjobj(H.f_cell); 
set(jobj,'MouseEnteredCallback','global H; figure(H.f_cell)');
set(jobj,'MouseClickedCallback',{@Update,'controls'});
H.p_cell_morph = uipanel('parent',H.f_cell,'Position',[.02 .65 .2 .25],'BackgroundColor','white','BorderWidth',.2,'BorderType','line'); % cell morphology
H.p_cell_spec = uipanel('parent',H.f_cell,'Position',[.263 .65 .47 .35],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','rightclick mechanism for model','fontangle','italic'); % cell specification
H.txt_cellname = uicontrol('parent',H.f_cell,'units','normalized','position',[.02+.5*.2-.01-.06-.03 .95 .05 .03],'BackgroundColor','white','style','edit','string',cname,'fontsize',14,'callback',@changename);
H.btn_undo = uicontrol('parent',H.f_cell,'units','normalized','position',[.02+.5*.2+.05-.06 .95 .04 .03],'style','pushbutton','string','undo','fontsize',10,'callback',@undo);
H.btn_print = uicontrol('parent',H.f_cell,'units','normalized','position',[.02+.5*.2+.05-.01 .95 .04 .03],'style','pushbutton','string','print','fontsize',10,'callback',@printmodel);
H.btn_close = uicontrol('parent',H.f_cell,'units','normalized','position',[.02+.5*.2+.05+.04 .95 .04 .03],'style','pushbutton','string','close','fontsize',10,'callback','global currspec; currspec=[]; close;');
H.p_cell_parms = uipanel('parent',H.f_cell,'Position',[.75 .65 .23 .35],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','Parameters'); % cell specification
H.p_static_plots = uipanel('parent',H.f_cell,'Position',[.02 .05 .35 .55],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','static'); % cell specification
H.p_simulation = uipanel('parent',H.f_cell,'Position',[.4 .05 .57 .55],'BackgroundColor','white','BorderWidth',.2,'BorderType','line','title','dynamic'); % cell specification

DrawControls;

%% Callback Functions
%{
------------------
Callback functions
------------------
Display_Mech_Info(mechpos,uitype): print model equations for mechanism w/ focus
  calls: mech_spec2str, mech_str2spec
Update(what,compname,mechname): switch calls UpdateControls and UpdatePlots.UpdateControls().UpdatePlots().
UpdateMechList(): DONE
AddMech(): DONE
RemoveMech(): DONE
UpdateCellMech(): need to rethink mech model handling
AddComp()
RemoveComp()
?Revise_Params()?
-------------------
Helper functions
-------------------
updatecellmodel(): change the model structure (spec) and pass through buildmodel2()
get_neighb(n,x,y): finds linear identifiers of neighbors of (x,y)
get_connected(ind): finds linear indices in spec.connections() connected to spec.cells(ind)
mech_spec2str(): convert spec.cells(#).mechs(#) to txt{:}
mech_str2spec(): convert txt{:} to updated spec.cells(#).mechs(#)
maybe: data2str(): convert all data classes to char
%}

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Display_Mech_Info(src,evnt,mechpos,uitype,hmech,hcomp)
% purpose: display the mech model in a readable form
try
global H currspec cfg
try
  H.ui_mechlist_edit=H.ui_mechlist_edit(ishandle(H.ui_mechlist_edit));
end
if any(strcmp('on',get(H.ui_mechlist_edit,'visible'))), return; end
cfg.focuscomp = strmatch(get(hcomp,'string'),{currspec.cells.label},'exact');
spec = currspec.cells(cfg.focuscomp);
cfg.focusmech = strmatch(get(hmech,'string'),spec.mechanisms,'exact');
mech = spec.mechs(cfg.focusmech);
compname = spec.label;
mechname = spec.mechanisms{cfg.focusmech};
showname = [compname '.' mechname ' (click text to edit)'];
Update(src,evnt,'controls',compname,mechname);
fpos=get(H.f_cell,'position');
panpos=get(H.p_cell_spec,'position'); lo=.5;
pos(1)=fpos(1)+fpos(3)*panpos(1);
pos(2)=fpos(2)+fpos(4)*(panpos(2)-lo);
pos(3)=fpos(3)*panpos(3);
pos(4)=fpos(4)*lo+.86*fpos(4)*mechpos(2)*panpos(4);
% pos(3)=.75*pos(3);
% pos(4)=.75*pos(4);
% pos(1)=pos(1)+.5*pos(3);
% pos(2)=pos(2)+pos(4);
if isfield(H,'f_hover') && ishandle(H.f_hover)
  figure(H.f_hover); 
  set(gcf,'name',showname,'position',pos);
else
  H.f_hover = figure('MenuBar','none','name',showname,'NumberTitle','off','position',pos);  
end
% -------------------------------------
% CHANGELOG: now forcing uitype to be editable. remove this line to restore click-to-edit
uitype = 'edit';
% -------------------------------------
% show mech info in uicontrol determined by uitype
clf
fontname = 'courier'; % fixedwidth, courier, helvetica (default)
H.mech_txt = uicontrol('style',uitype,'string','','FontName',fontname,'FontSize',10,...
  'BackgroundColor',[.9 .9 .9],'Max',50,'HorizontalAlignment','Left','unit','normalized','Position',[0 0 1 1]);
try
  jobj=findjobj(H.mech_txt);
  set(jobj,'MouseClickedCallback',{@Display_Mech_Info,mechpos,'edit',hmech,hcomp});
end
if strcmp(uitype,'edit')
  set(H.f_hover,'MenuBar','none'); 
  mech_apply = uimenu(H.f_hover,'Label','Apply','Callback',{@UpdateCellMech,H.mech_txt,compname,mechname});
  file_write = uimenu(H.f_hover,'Label','Write','Callback',{@SaveCellMech,H.mech_txt,compname,mechname});
  set(H.mech_txt,'BackgroundColor','w');
end
txt = mech_spec2str(mech);
set(H.mech_txt,'string',txt,'position',[0 0 1 1]);
catch
  %disp('ERROR!');
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateMechList(src,evnt,action)
% apply user changes to the mechs included in a compartment model
global H cfg currspec
mechlist = currspec.cells(cfg.focuscomp).mechanisms;
switch action
  case 'toggle'
    if strcmp(get(H.ui_mechlist_edit(cfg.focuscomp),'visible'),'off')
      str=''; for j=1:length(mechlist),str=[str mechlist{j} ', ']; end; str=str(1:end-2);
      set(H.ui_mechlist_edit(cfg.focuscomp),'visible','on','string',str);%,'units','normalized','position',[.1 .1 .2 .2]);
      set(H.ui_comp_mechs(cfg.focuscomp,ishandle(H.ui_comp_mechs(cfg.focuscomp,:))),'visible','off');
    else
      set(H.ui_mechlist_edit(cfg.focuscomp),'visible','off');
      set(H.ui_comp_mechs(cfg.focuscomp,ishandle(H.ui_comp_mechs(cfg.focuscomp,:))),'visible','on');
    end
  case 'update'
    if strcmp(evnt.Key,'return')
      pause(.5); 
      newspec = currspec;
      set(H.ui_mechlist_edit(cfg.focuscomp),'visible','off');
      str = get(H.ui_mechlist_edit(cfg.focuscomp),'string');
      if ~isempty(str)
        newmechlist = strtrim(splitstr(str,','));
        mechadded = setdiff(newmechlist,mechlist);
        mechremoved = setdiff(mechlist,newmechlist);
        if ~isempty(mechadded) || ~isempty(mechremoved)
          newspec = removemech(newspec,mechremoved);
          newspec = addmech(newspec,mechadded);
          updatecellmodel(newspec);
        end
        DrawControls;
      end
    else
      % do nothing
    end
end

function SaveCellMech(src,evnt,htxt,compname,mechname)
global currspec
% purpose: write mech to new file
outfile = '';
txt = get(htxt,'string');
%thiscell = strmatch(compname,{currspec.cells.label},'exact');
%comp = currspec.cells(thiscell);
%thismech = strmatch(mechname,comp.mechanisms,'exact');

% - popup dialog box
% - use fprintf to write each line of txt{:}
% - store new mech name in spec and update the model

% currspec.files{end+1} = outfile;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateCellMech(src,evnt,htxt,compname,mechname)
% purpose: apply user changes to the mech model
global currspec
spec=currspec;
txt = get(htxt,'string');
newmech = parse_mech_spec(txt);
newmech.label = mechname;
thiscell = strmatch(compname,{spec.cells.label},'exact');
comp = spec.cells(thiscell);
thismech = strmatch(mechname,comp.mechanisms,'exact');
if ~isempty(thismech)
  comp.mechs(thismech) = newmech;
else
  comp.mechs(end+1) = newmech;
  comp.mechanisms{end+1} = mechname;
end
spec.cells(thiscell) = comp;
updatecellmodel(spec);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateCell(src,evnt,xi,yi)
% purpose: add a new compartment (or remove an existing compartment)
CONNECTION_MECHANISM = 'iCOM';
global cfg H currspec
newspec = currspec;
nrows = sqrt(numel(H.ui_cellcomp));
% determine if a compartment was added or removed
contents = get(H.ui_cellcomp,'string');
gridlabels = contents(~cellfun(@isempty,contents));
lastlabels = {currspec.cells.label};
if ~isempty(setdiff(gridlabels,lastlabels))
  action = 'add';
elseif ~isempty(setdiff(lastlabels,gridlabels))
  action = 'remove';
else
  return;
end
switch action
  case 'add' % add a compartment
    % add new compartment
    newlabel = get(src,'string'); % contents{xi,yi}
    newcomp = currspec.cells(cfg.focuscomp);
    newspec.cells(end+1) = newcomp;
    newspec.cells(end).label = newlabel;
    thisind = length(newspec.cells);
    % add connections (ohmic currents) to/from the compartment:
    ind = find(~cellfun(@isempty,contents));
    [X,Y] = ind2sub([cfg.maxcomp cfg.maxcomp],ind);
    idx = abs(X-xi)<=1 & abs(Y-yi)<=1;
    if length(find(idx))>1 % will always be at least 1 b/c of self-self adjacency
      % add connections
      sel = sub2ind([cfg.maxcomp cfg.maxcomp],X(idx),Y(idx));
      labels = get(H.ui_cellcomp(sel),'string');
      others = setdiff(labels,newlabel);
      for i=1:length(others)
        % add bidirectional connections (iCOM) to each
        thatind = find(strcmp(others{i},{newspec.cells.label}));
        newspec.connections(thisind,thatind).label = [newlabel '-' others{i}];
        newspec.connections(thisind,thatind).mechanisms = CONNECTION_MECHANISM;
        newspec.connections(thatind,thisind).label = [others{i} '-' newlabel];
        newspec.connections(thatind,thisind).mechanisms = CONNECTION_MECHANISM;
      end
    end
    % update callbacks
    jobj=findjobj(H.ui_cellcomp(xi,yi));
    set(jobj,'MouseEnteredCallback',{@Update,'controls',newlabel});
  case 'remove' % remove a compartment
    %deleted = setdiff(lastlabels,gridlabels);
    %thatind = strmatch(deleted,{newspec.cells.label});
    newspec.cells(cfg.focuscomp) = [];
    newspec.connections(cfg.focuscomp,:) = [];
    newspec.connections(:,cfg.focuscomp) = [];
    cfg.focuscomp = 1;
    cfg.focusmech = 1;
    % update callbacks
    %jobj=findjobj(H.ui_cellcomp(xi,yi));
    %set(jobj,'MouseEnteredCallback',[]);
end
updatecellmodel(newspec);
DrawControls;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function updatecellmodel(newspec) % maybe use same specification as for "override"
% purpose: update the integratable model after each change to its specification
% note: should be called after any change to the multiscale model
% purpose to consider: input override matrix and update currspec in this function...
global currspec lastspec
lastspec = currspec;
currspec = newspec;
% remove all connections for which gCOM=0 before building the model
% ...
[model,IC,functions,auxvars,currspec] = buildmodel2(currspec,'verbose',0);
try currspec.cells=currspec.entities; currspec=rmfield(currspec,'entities'); end
fprintf('model updated successfully\n');

% redraw whatever's appropriate (ex. plots; TBD)
% ...

% How to update the model:
% - update global entity param in model: change param in spec.cells(1).parameters = {key,val} and pass through buildmodel2()
% - update mech param in model: change param in spec.cells(1).mechs(1).params.(key) and pass through buildmodel2()
% - update mech list in model: 
% 	add: add mech label to spec.cells(1).mechanisms{end+1} and file to spec.files{end+1}, then pass through buildmodel2()
% 	remove: remove mech label from spec.cells(1).mechanisms and model from spec.cells(1).mechs
% - update compartments/cells in model:
% 	add: add spec.cells(end+1) and give a unique spec.cells(#).label; pass through buildmodel2()
% 	remove: spec.cells(#)=[], spec.connections(:,#)=[], spec.connections(#,:)=[]

% % override parameters in specification (this works for label,N,mechlist,Pg)
% if ~isempty(parms.override)
%   o = parms.override;
%   [nover,ncols] = size(o);
%   for k = 1:nover
%     l = o{k,1}; f = o{k,2}; v = o{k,3}; 
%     if ncols>3, a = o{k,4}; else a = []; end
%     if ~ischar(l) || ~ischar(f), continue; end
%     if ismember(l,Elabels), type='entities';
%     elseif ismember(l,Clabels), type='connections';
%     else continue; 
%     end
%     n = strmatch(l,{spec.(type).label},'exact');
%     if isequal(f,'parameters')
%       if isempty(spec.(type)(n).(f)), continue; end
%       matched = find(cellfun(@(x)isequal(v,x),spec.(type)(n).(f)));
%       if ~isempty(matched)
%         spec.(type)(n).(f){matched+1} = a;
%       else
%         spec.(type)(n).(f){end+1} = v;
%         spec.(type)(n).(f){end+1} = a;
%       end
%     else
%       spec.(type)(n).(f) = v;
%     end
%   end
%   clear o nover ncols l f v a n
% end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function undo(src,evnt)
% revert to the last working model
global lastspec
updatecellmodel(lastspec);
DrawControls;

function printmodel(src,evnt)
global currspec
buildmodel2(currspec,'verbose',1);

function changename(src,evnt)
global currspec
name=get(src,'string');
[currspec.cells.parent]=deal(name);
%(1:length(currspec.cells))

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Update(src,evnt,what,compname,mechname)
% purpose: handle updating graphics
global H cfg currspec
if cfg.focuscomp>length(currspec.cells), cfg.focuscomp=1; end
if cfg.focusmech>length(currspec.cells(cfg.focuscomp).mechanisms), cfg.focusmech=1; end
if nargin < 4, compname=currspec.cells(cfg.focuscomp).label; end
if nargin < 5
  if isempty(cfg.focusmech) || cfg.focusmech==0
    mechname='';
  else
    mechname=currspec.cells(cfg.focuscomp).mechanisms{cfg.focusmech}; 
  end
end
cfg.focuscomp = find(strcmp({currspec.cells.label},compname));
cfg.focusmech = find(strcmp(currspec.cells(cfg.focuscomp).mechanisms,mechname));
switch what
  case 'controls'
    set(H.ui_cellcomp,'BackgroundColor','w');
    set(H.ui_cellcomp(strcmp(get(H.ui_cellcomp,'string'),compname)),'BackgroundColor',cfg.focuscolor);
    set(H.ax_comp_area(ishandle(H.ax_comp_area)),'linewidth',.5,'color','w');
    set(H.ax_comp_area(cfg.focuscomp),'linewidth',3,'color',cfg.focuscolor);
    set(H.ax_comp_hover(ishandle(H.ax_comp_hover)),'BackgroundColor','w');
    set(H.ax_comp_hover(cfg.focuscomp),'BackgroundColor',cfg.focuscolor);
    %set(H.ui_mechlist_edit(cfg.focuscomp),'visible','off');
    %set(H.ui_comp_mechs(cfg.focuscomp,ishandle(H.ui_comp_mechs(cfg.focuscomp,:))),'visible','on');
    set(H.p_cell_parms,'Title',currspec.cells(cfg.focuscomp).label,'FontWeight','bold');
    % update parameter controls
    % intrinsic parameters
    p = currspec.cells(cfg.focuscomp).parameters;
    if ~isempty(p)
      intparms = p(1:2:end); valind=1; val=num2str(p{2});
    else
      intparms = ''; valind=[]; val='';
    end
    set(H.ui_cells_paramlist,'string',intparms,'value',valind);
    set(H.ui_cells_paramedit,'string',val);
    % connection parameters
    conninds=get_connected(currspec,cfg.focuscomp);
    connparms = {}; val=''; valind=[];
    for i=1:length(conninds)
      this = currspec.connections(conninds(i),cfg.focuscomp);
      if isempty(this.parameters), continue; end
      valind=1; val=num2str(this.parameters{2});
      thisp = this.parameters(1:2:end);
      for j=1:length(thisp), thisp{j}=[this.label '.' thisp{j}]; end
      connparms = {connparms{:} thisp{:}};
    end
    set(H.ui_connections_paramlist,'string',connparms,'value',valind);
    set(H.ui_connections_paramedit,'string',val);   
end
% update static plots
% TODO: coordinate (+consolidate) this with the matching code in DrawControls & plotfunctions
    maxlhs = 20; % limit how much is shown in the listbox
    maxlen = 150;
    funcs = currspec.cells(cfg.focuscomp).functions;
    len = min(maxlhs,max(cellfun(@length,funcs(:,1))));
    str = {};
    for i=1:size(funcs,1)
      str{i} = sprintf(['%-' num2str(len) 's  = %s'],funcs{i,1},strrep(funcs{i,2},' ',''));
      if length(str{i})>maxlen, str{i}=str{i}(1:maxlen); end
    end
    val=get(H.lst_static_funcs,'value');
    set(H.lst_static_funcs,'string',str);
    set(H.lst_static_funcs,'value',val(val<=length(str)));
    plotfunctions;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function currspec = addmech(currspec,mechadded)
% purpose: add a mechanism to the compartment model
global allmechs cfg
if isempty(mechadded), return; end
if ~iscell(mechadded), mechadded={mechadded}; end
for i=1:length(mechadded)
  newmech = mechadded{i};
  newmech = find(strcmp({allmechs.label},newmech)); mechind=newmech;
  newmech = allmechs(newmech);
  if isempty(newmech)
    warndlg([mechadded{i} ' not found. Check spelling and case.']);
    disp('known mechanisms include: '); disp(get_mechlist');
  else
    %newfld = fieldnames(newmech);
    %oldfld = fieldnames(currspec.cells(cfg.focuscomp).mechs);
    newmech = rmfield(newmech,'file');
    currspec.cells(cfg.focuscomp).mechs(end+1)=newmech;
    currspec.cells(cfg.focuscomp).mechanisms{end+1}=mechadded{i};
    currspec.files{end+1} = cfg.allmechfiles{mechind};
  end
end
cfg.focusmech = length(currspec.cells(cfg.focuscomp).mechs);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function currspec = removemech(currspec,mechremoved)
% purpose: remove a mechanism from the compartment model
global cfg
if isempty(mechremoved), return; end
if ~iscell(mechremoved), mechremoved={mechremoved}; end
for i=1:length(mechremoved)
  oldmech = mechremoved{i};
  oldmech = find(strcmp(currspec.cells(cfg.focuscomp).mechanisms,oldmech));
  if isempty(oldmech)
    warndlg([mechremoved{i} ' not found. Unknown error.']);
    fprintf('Current mechanisms for %s include: ',currspec.cells(cfg.focuscomp).label)
    currspec.cells(cfg.focuscomp).mechanisms
  else
    currspec.cells(cfg.focuscomp).mechs(oldmech)=[];
    currspec.cells(cfg.focuscomp).mechanisms(oldmech)=[];
  end
end
cfg.focusmech = min(1,length(currspec.cells(cfg.focuscomp).mechs));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function DrawControls
% purpose: draw/redraw uicontrols with properties reflecting the most recent model (currspec)
global H cfg currspec
% cleanup
if isfield(H,'f_hover')
  try close(H.f_hover); end
  H=rmfield(H,'f_hover');
end
delete(get(H.p_cell_morph,'Children'));
delete(get(H.p_cell_spec,'Children'));
try
fld='ui_comp_mechs';    delete(H.(fld)(ishandle(H.(fld)))); H=rmfield(H,fld);
fld='ui_comp_edit';     delete(H.(fld)(ishandle(H.(fld)))); H=rmfield(H,fld);
fld='ui_mechlist_edit'; delete(H.(fld)(ishandle(H.(fld)))); H=rmfield(H,fld);
fld='ui_comp_txt';      delete(H.(fld)(ishandle(H.(fld)))); H=rmfield(H,fld);
fld='ax_comp_area';     delete(H.(fld)(ishandle(H.(fld)))); H=rmfield(H,fld);
fld='ax_comp_hover';    delete(H.(fld)(ishandle(H.(fld)))); H=rmfield(H,fld);
end

spec=currspec;
compnames={spec.cells.label}; % list root first then children nodes
ncomp = numel(compnames);
maxmechs = numel(cfg.selmechlist);
% add morphology boxes
maxcomp = 9; cfg.maxcomp=maxcomp;
w = 1/(maxcomp-1);
x = 0:w:1; cnt=0;
for i=1:maxcomp
  for j=1:maxcomp
    cnt=cnt+1;
    if cnt==(maxcomp^2+1)/2
      str = '';
      bg = cfg.focuscolor;
    else
      str = '';
      bg = 'white';
    end
    H.ui_cellcomp(i,j) = uicontrol('parent',H.p_cell_morph,'style','edit','string','',...
      'units','normalized','position',[x(i) x(j) w w],'BackgroundColor',bg,...
      'Callback',{@UpdateCell,i,j});
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

compcnt=1; compind=1; linked{1}=compnames{1};
gridnums=flipud(reshape(1:maxcomp^2,[maxcomp maxcomp]));
xi=round((maxcomp+1)/2); yi=xi; allxi=xi; allyi=yi;
set(H.ui_cellcomp(xi,yi),'string',compnames{compcnt});
jobj=findjobj(H.ui_cellcomp(xi,yi));
set(jobj,'MouseEnteredCallback',{@Update,'controls',compnames{compcnt}});
while compcnt<ncomp
  if compcnt>1
    xi=allxi(compcnt);
    yi=allyi(compcnt);
  end
  % what is around this compartment?
  neighbs = get_neighb(maxcomp^2,xi,yi);
  % who does it connect to?
  conninds = get_connected(spec,compind);
  if isempty(conninds)
    jj=ceil(compcnt/2/maxcomp);
    tmp=find(cellfun(@isempty,get(H.ui_cellcomp(:,jj),'string'))==1);
    ii=tmp(2*compcnt);
    set(H.ui_cellcomp(ii,jj),'string',compnames{compind});
    allxi=[allxi ii]; allyi=[allyi jj];
  end
  % where do they go?
  for k=1:length(conninds) % loop over connections to this compartment
    if ismember(compnames{conninds(k)},get(H.ui_cellcomp,'string'))
      continue;
    end
    for i=1:length(neighbs) % loop over grid neighb of this compartment
      [ii,jj]=find(gridnums==neighbs(i));
      if ~isempty(get(H.ui_cellcomp(ii,jj),'string')), continue; end
      targetneighbs = get_neighb(maxcomp^2,ii,jj);
      targetconns = get_connected(spec,conninds(k));
      tmp=cellfun(@(x)find(x==gridnums),num2cell(targetneighbs));
      [II,JJ]=ind2sub(size(gridnums),tmp);
      lbls=cellfun(@(a,b)(get(H.ui_cellcomp(a,b),'string')),num2cell(II),num2cell(JJ),'unif',0);
      lbls=lbls(~cellfun(@isempty,lbls));
      if all(ismember(lbls,compnames(targetconns)))
        set(H.ui_cellcomp(ii,jj),'string',compnames{conninds(k)});
        jobj=findjobj(H.ui_cellcomp(ii,jj));
        set(jobj,'MouseEnteredCallback',{@Update,'controls',compnames{conninds(k)}});        
        linked{end+1}=compnames{conninds(k)};
        allxi=[allxi ii]; allyi=[allyi jj];
        break;
      end
    end
  end
  compcnt=compcnt+1;
  compind=compcnt;
end   

%% add cell specification area

l=.5; mw=l/(maxmechs+3); mx=0:mw:l;
for i=1:ncomp
  comp=spec.cells(i);
  % border around compartment specification
  if i==1
    H.ax_comp_area(i) = subplot('position',[.01 x(end-i)-w/10 .98 w/1.3],...
      'parent',H.p_cell_spec,'linewidth',3,'color',cfg.focuscolor); set(gca,'xtick',[],'ytick',[]); box on 
    H.ax_comp_hover(i) = uicontrol('parent',H.p_cell_spec,'style','text','BackgroundColor',cfg.focuscolor,...
      'units','normalized','position',[.01 x(end-i)-w/10 .98 w/1.3],'string','');
  else
    H.ax_comp_area(i) = subplot('position',[.01 x(end-i)-w/10 .98 w/1.3],...
      'parent',H.p_cell_spec,'linewidth',.5,'color','w'); set(gca,'xtick',[],'ytick',[]); box on 
    H.ax_comp_hover(i) = uicontrol('parent',H.p_cell_spec,'style','text','BackgroundColor','w',...
      'units','normalized','position',[.01 x(end-i)-w/10 .98 w/1.3],'string','');
  end
  jobj=findjobj(H.ax_comp_hover(i)); 
  set(jobj,'MouseEnteredCallback',{@Update,'controls',comp.label});
  %set(jobj,'PropertyChangeCallback',{@UpdateMechList,'update'});
  %set(jobj,'FocusLostCallback',{@UpdateMechList,'toggle'});
  % hidden edit control for updating the comma-separated list of mechs for this compartment
  H.ui_mechlist_edit(i) = uicontrol('parent',H.p_cell_spec,'style','edit','string','','units','normalized',...
    'visible','off','position',[.25+mx(1) x(end-i)-w/10 .98-.5 w/1.3],'BackgroundColor','w','HorizontalAlignment','left',...
    'KeyPressFcn',{@UpdateMechList,'update'});
  % compartment label in specification area
  H.ui_comp_txt(i) = uicontrol('parent',H.p_cell_spec,'style','text','string',comp.label,...
    'units','normalized','position',[.025 x(end-i) .15 w/2],'tooltip','left-click to edit mechanism list','ButtonDownFcn',{@UpdateMechList,'toggle'});
  jobj=findjobj(H.ui_comp_txt(i)); 
  %set(jobj,'MouseEnteredCallback',{@Update,'controls',comp.label});
  set(jobj,'MouseClickedCallback',{@UpdateMechList,'toggle'});
  % edit control for setting compartment dynamics
  H.ui_comp_edit(i) = uicontrol('parent',H.p_cell_spec,'style','edit','string',[comp.dynamics{:}],...
    'units','normalized','position',[.8 x(end-i) .18 w/2],'BackgroundColor','w','HorizontalAlignment','left','tooltipstring',[spec.cells(i).ode_labels{1} ''' = ' spec.cells(i).odes{1}]);%'Mechanism functions will be substituted here.');
  %jobj=findjobj(H.ui_comp_edit(i)); set(jobj,'MouseEnteredCallback',{@Update,'controls',comp.label});
  % list of mechanisms for this compartment
  for j=1:length(comp.mechs) 
    mech=comp.mechs(j);
    pos=[.25+mx(j) x(end-i) .9*mw w/2];
    pos=[.25+mx(j) x(end-i) .9*mw w/2];
    if ~isempty(mech.substitute)
      tooltip='';
      for k=1:size(mech.substitute,1)
        tooltip=[tooltip mech.substitute{k,1} ' => ' mech.substitute{k,2} ' ']; 
      end
    else
      tooltip = 'WARNING: This mechanism is disconnected from the compartmental dynamics. Need to add substitution term "==>".';
    end
    H.ui_comp_mechs(i,j) = uicontrol('parent',H.p_cell_spec,'style','text','string',spec.cells(i).mechanisms{j},...
      'units','normalized','position',pos,'TooltipString',tooltip,'ForegroundColor','b');
    set(H.ui_comp_mechs(i,j),'ButtonDownFcn',{@Display_Mech_Info,pos,'text',H.ui_comp_mechs(i,j),H.ui_comp_txt(i)});
    %jobj=findjobj(H.ui_comp_mechs(i,j));
    %try jobj.setCursor(java.awt.Cursor.getPredefinedCursor(java.awt.Cursor.HAND_CURSOR)); end
    %set(jobj,'MouseEnteredCallback',{@Display_Mech_Info,pos,'text',H.ui_comp_mechs(i,j),H.ui_comp_txt(i)});% sprintf('disp(%g)',j)); % display mech info in new fig
    %set(jobj,'MouseClickedCallback',{@Display_Mech_Info,pos,'text',H.ui_comp_mechs(i,j),H.ui_comp_txt(i)});
    %set(jobj,'MouseClickedCallback',{@UpdateMechList,'toggle'}); % open editable textbox with complete mech file
    % note: if new fig save button is pressed, a new mech file should be written and the modeler figure updated to include the new info throughout
  end
end

%% add parameter boxes for global intrinsic and connection parameters
set(H.p_cell_parms,'Title',currspec.cells(cfg.focuscomp).label,'FontWeight','bold');

% intrinsic parameters
p = currspec.cells(cfg.focuscomp).parameters;
if ~isempty(p)
  intparms = p(1:2:end); valind=1; val=num2str(p{2});
else
  intparms = ''; valind=[]; val='';
end
dy=-.015;
H.ui_cells_paramlabel = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','text','position',[.03 .94+dy .25 .08],'backgroundcolor','w','string','intrinsic','HorizontalAlignment','left');
H.ui_cells_paramlist = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','listbox','position',[.03 .18+dy .25 .78],'value',valind,'string',intparms,...
  'backgroundcolor','w','Max',1,'Min',0,'Callback',{@UpdateParams,'show','cells'});
H.ui_cells_paramedit = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[.03 .11+dy .25 .065],'backgroundcolor','w','string',val,...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'change','cells'});
H.ui_cells_paramadd = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[.03 .03+dy .25 .065],'backgroundcolor','w','string','key = value',...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'add','cells'});

% connection parameters
conninds=get_connected(currspec,cfg.focuscomp);
connparms = {}; val=''; valind=[];
for i=1:length(conninds)
  this = currspec.connections(conninds(i),cfg.focuscomp);
  if isempty(this.parameters), continue; end
  if i==1, valind=1; val=num2str(this.parameters{2}); end
  thisp = this.parameters(1:2:end);
  for j=1:length(thisp), thisp{j}=[this.label '.' thisp{j}]; end
  connparms = {connparms{:} thisp{:}};
end
H.ui_connections_paramlabel = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','text','position',[.33 .94+dy .6 .08],'backgroundcolor','w','string','connections','HorizontalAlignment','left');
H.ui_connections_paramlist = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','listbox','position',[.33 .18+dy .6 .78],'value',valind,'string',connparms,...
  'backgroundcolor','w','Max',1,'Min',0,'Callback',{@UpdateParams,'show','connections'});
H.ui_connections_paramedit = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[.33 .11+dy .6 .065],'backgroundcolor','w','string',val,...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'change','connections'});
H.ui_connections_paramadd = uicontrol('parent',H.p_cell_parms,'units','normalized',...
  'style','edit','position',[.33 .03+dy .6 .065],'backgroundcolor','w','string','key = value',...
  'HorizontalAlignment','left','Callback',{@UpdateParams,'add','connections'});

%% create plots
cfg.T=(0:cfg.buffer-1)*cfg.dt; % ms
cfg.V=linspace(-100,100,cfg.buffer);%(-100:.01:100); % mV
cfg.colors  = 'kbrgmy';
cfg.lntype  = {'-',':','-.','--'};
%usage: lnstyle = [cfg.colors(mod(outer,length(cfg.colors))) cfg.lntype{mod(inner,length(cfg.colors))}];

% dynamic plots (simulated data)
H.ax_sim_plot = subplot('position',[.05 .25 .7 .7],...
  'parent',H.p_simulation,'linewidth',3,'color','w'); 
box on; % set(gca,'xtick',[],'ytick',[]);
legstr = {};
for k=1:ncomp
  H.simdat_alltrace(k)=line('color',cfg.colors(max(1,mod(k,length(cfg.colors)))),'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',cfg.T,'ydata',zeros(1,cfg.buffer),'zdata',[]);
  legstr{k} = currspec.cells(k).label;
end
axis([cfg.T(1) cfg.T(end) -100 30]); xlabel('time (ms)');
h=legend(legstr{:});

% TODO: add listbox of statevars to plot

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
% if isempty(findobj('tag','start'))
%   % btn: start <=> reset        
%   uicontrol('Style','pushbutton', 'Units','normalized', ...
%             'Position',[0.85  0.19 0.1 0.05],...
%             'String','start','tag','start','Callback',{@simulate,'restart'}); % start <=> pause
% end
% if isempty(findobj('tag','pause'))
%   % btn: pause <=> resume
%   uicontrol('Style','pushbutton', 'Units','normalized', ...
%             'Position',[0.85  0.12 0.1 0.05],...
%             'String','pause','tag','pause','Callback','global cfg;cfg.pauseflag;cfg.pauseflag=-cfg.pauseflag;');
% end
% if isempty(findobj('tag','stop'))
%   % btn: stop
%   uicontrol('Style','pushbutton', 'Units','normalized', ...
%             'Position',[0.85  0.05 0.1 0.05],...
%             'String','stop','tag','stop','Callback','global cfg;cfg.quitflag=1;');
% end
% % autoscale       
% uicontrol('Style','pushbutton', 'Units','normalized', ...
%           'Position',[0.85 0.5 0.1 0.05],...
%           'String','autoscale','Callback',{@setlimits,'autoscale'});
% 
% end
if isempty(findobj('tag','start'))
    uicontrol('Style','pushbutton', 'Units','normalized', ...
      'Position',[0.75  0.05 0.05 0.05],...
      'String','start','tag','start','Callback',{@simulate,'restart'}); % start <=> pause
end
if isempty(findobj('tag','pause'))
  % btn: pause <=> resume
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.80  0.05 0.05 0.05],...
            'String','pause','tag','pause','Callback','global cfg;cfg.pauseflag;cfg.pauseflag=-cfg.pauseflag;');
end
if isempty(findobj('tag','stop'))
  % btn: stop
  uicontrol('Style','pushbutton', 'Units','normalized', ...
            'Position',[0.85  0.05 0.05 0.05],...
            'String','stop','tag','stop','Callback','global cfg;cfg.quitflag=1;');
end
% autoscale       
uicontrol('Style','pushbutton', 'Units','normalized', ...
          'Position',[0.9 0.05 0.075 0.05],...
          'String','autoscale','Callback',{@setlimits,'autoscale'});

H.ax_static_plot = subplot('position',[.04 .45 .9 .5],'parent',H.p_static_plots,'linewidth',3,'color','w'); box on; 
H.lst_static_funcs = uicontrol('units','normalized','position',[.04 .02 .9 .35],'parent',H.p_static_plots,...
  'style','listbox','value',1:5,'string',str,'Max',50,'Callback',@plotfunctions);
uicontrol('Style','edit', 'Units','normalized','Position',[0.75 0.45 0.25 0.05],...
          'String',sprintf('[%g,%g]',min(cfg.V),max(cfg.V)),'Callback',[],'parent',H.p_static_plots);
uicontrol('Style','pushbutton', 'Units','normalized','Position',[0.75 0.9 0.1 0.05],...
          'String','fit?','Callback',[],'parent',H.p_static_plots);

plotfunctions;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotfunctions(src,evnt)
global H cfg currspec
if isfield(H,'static_traces')
  try delete(H.static_traces); end
  H=rmfield(H,'static_traces'); 
  cla(H.ax_static_plot);
end
% maxlhs = 20; % limit how much is shown in the listbox
% maxlen = 150;
% funcs = currspec.cells(cfg.focuscomp).functions;
% len = min(maxlhs,max(cellfun(@length,funcs(:,1))));
% str = {};
% for i=1:size(funcs,1)
%   str{i} = sprintf(['%-' num2str(len) 's  = %s'],funcs{i,1},strrep(funcs{i,2},' ',''));
%   if length(str{i})>maxlen, str{i}=str{i}(1:maxlen); end
% end
% set(H.lst_static_funcs,'string',str);
sel = get(H.lst_static_funcs,'value');
list = get(H.lst_static_funcs,'string');
functions = currspec.cells(cfg.focuscomp).functions(sel,:);
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
  Y = f(X);
  H.static_traces(k)=line('parent',H.ax_static_plot,'color',cfg.colors(max(1,mod(k,length(cfg.colors)))),...
    'LineStyle',cfg.lntype{max(1,mod(k,length(cfg.lntype)))},'erase','background','xdata',X,'ydata',Y,'zdata',[]);
end
h=legend(functions(:,1)); set(h,'fontsize',6,'location','EastOutside');
axis tight; set(gca,'fontsize',6);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setlimits(src,evnt,action)
global cfg H
switch action
  case 'autoscale'
    ymin = min(cfg.record(:));
    ymax = max(cfg.record(:));
    if ymin~=ymax
      set(H.ax_sim_plot,'ylim',[ymin ymax]);
    end
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function simulate(src,evnt,action)
global currspec cfg H
functions = currspec.model.functions;
auxvars = currspec.model.auxvars;
ode = currspec.model.ode;
IC = currspec.model.IC;
dt = cfg.dt;
buffer = cfg.buffer;
T = cfg.T;

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
X=IC; t=0; cnt=0; cfg.record=zeros(length(IC),buffer);
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
  X = X+dt*F(t,X);
  t = t + dt;
  if cnt<=buffer
    cfg.record(:,cnt)=X;
  else
    cfg.record = [cfg.record(:,2:end) X];
  end
 
  % get output to plot
  for k=1:length(currspec.cells)
    if k>length(H.simdat_alltrace)
      break;
    end
    var = currspec.cells(k).ode_labels{1};
    var = find(cellfun(@(x)isequal(x,var),currspec.cells(k).var_list));
    set(H.simdat_alltrace(k),'ydata',cfg.record(var,:));
  end
  drawnow
end
cfg.quitflag=-1;
p=findobj('tag','pause');
set(p,'string','pause'); 
set(findobj('tag','start'),'string','restart');

function UpdateParams(src,evnt,action,type) 
% note: type=connections|cells, and is used w/ srcnum to get params from
% cells or connections fields using the same statements.
global cfg currspec H
h = H.(['ui_' type '_paramlist']);
g = H.(['ui_' type '_paramedit']);
f = H.(['ui_' type '_paramadd']);
newspec = currspec;
list = get(h,'string');
if strcmp(type,'connections')
  if strcmp(action,'show') || strcmp(action,'change')
    label = splitstr(list{get(h,'value')},'.');
  elseif strcmp(action,'add')
    label = regexp(get(f,'string'),'^.+=','match');
    label = splitstr(strtrim(label{1}(1:end-1)),'.');
  end
  srcnum = find(cellfun(@(x)isequal(x,label{1}),{currspec.(type)(:,cfg.focuscomp).label}));
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
    updatecellmodel(newspec);
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
    updatecellmodel(newspec);
end

