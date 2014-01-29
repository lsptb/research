function [T_mri_head, T_head_mri] = ts_pointreg(hptsfile,surffile,varargin);
%function [T_mri_head, T_head_mri] = ts_pointreg(hptsfile,surffile,varargin);
%
% Purpose: Manual registration between MEG and MRI reference frames
%   Plots head points in relation to scalp surface vertices and allows for
%   rotation and translation in 3 axes
%
% Usage:
%   [T_mri_head, T_head_mri] = ts_pointreg(hptsfile,surffile,'key1', value1,...);
%
% Required input:
%  hptsfile - text file containing head point coordinates
%  surffile - scalp surface file in "tri" format
%
% Optional parameters:
%  outtransfile - file name for output 4x4 transformation matrix
%    {default = ts_pointreg.trans}
%  intransfile - file name for input 4x4 transformation matrix
%    {default = use identity matrix}
%  intrans - 4x4 transformation matrix
%    {default = use identity matrix}
%  trans_type - ['head2mri' or 'mri2head']
%     if 'head2mri' - input and output transformation matrices are head2mri
%       i.e. head  points registered to surface (derived from MRI)
%     if 'mri2head' - input and output transformation matrices are mri2head
%    {default = 'mri2head'}
%
%  colorval - value between 1 and 100 that determines color of surface mesh
%    {default: 10} 
%  face_alpha - value between 0 and 1 that determines transparency of triangles
%    {default: 0.5} 
%  edge_alpha - value between 0 and 1 that determines transparency of edges
%    {default: 0.1} 
%
%
% example of hpts file format:
%cardinal 001 -0.070086 0.000000 -0.000000
%cardinal 002 -0.000000 0.099511 -0.000000
%cardinal 003 0.064181 0.000000 -0.000000
%hpi      001 0.056898 -0.018247 -0.036789
%hpi      002 0.038970 0.098632 0.046132
%hpi      003 -0.018337 0.102662 0.058440
%hpi      004 -0.060982 -0.012347 -0.029710
%extra    001 0.022850 0.103267 0.048321
%extra    002 0.019972 0.101263 0.061493
%extra    003 0.018250 0.096633 0.074459
%extra    004 0.017615 0.089393 0.085852
%
% the hptsfile can be obtained by running fiff2ascii -p on Neuromag fif file
%
% example of trans file format:
%1.000000 0.000000 0.000000 0.000000 
%0.000000 1.000000 0.000000 0.000000 
%0.000000 0.000000 1.000000 0.000000 
%0.000000 0.000000 0.000000 1.000000 
%
% units of the transformation matrix should be meters
%   if old .trans file created with tkmedit is used, the transformation
%   matrix should be loaded first with read_transfile like this:
%   intrans = ts_read_transfile(fname,0.001);
%   to convert units from millimeters to meters.
%   Then run:
%     ts_pointreg(hptsfile,surffile,'intrans',intrans)
%
%  created:       05/02/06   by Don Hagler
%  last modified: 08/07/09   by Don Hagler
%

%% todo: use Uutela's hpipoints.m

if (~mmil_check_nargs(nargin,2)) return; end;

DEFAULT_AXES_POSITION = [0.03 0.12 0.94  0.85];
DEFAULT_FIG_COLOR = [0.6 0.6 0.6];
BUTTON_COLOR =[0.8 0.8 0.8];
DEFAULT_CARD_POINT_COLOR = 'b';
DEFAULT_POINT_COLOR = 'r';

DEFAULT_TRANSLATION = 1;
DEFAULT_ROTATION = 1;

DEFAULT_VIEW_AZ = -150;
DEFAULT_VIEW_EL = 20;

DEFAULT_COLORVAL = 10;
DEFAULT_FACE_ALPHA = 0.5;
DEFAULT_EDGE_ALPHA = 0.1;

DEFAULT_OUTTRANSFILE = 'ts_pointreg.trans';
DEFAULT_TRANS_TYPE = 'mri2head';

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), data=struct(options{:}); 
  else data = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

try, data.outtransfile; catch, data.outtransfile = DEFAULT_OUTTRANSFILE; end;
try, data.intransfile;  catch, data.intransfile  = []; end;
try, data.intrans;      catch, data.intrans = eye(4); end;
try, data.colorval;     catch, data.colorval = DEFAULT_COLORVAL; end;
try, data.face_alpha;   catch, data.face_alpha = DEFAULT_FACE_ALPHA; end;
try, data.edge_alpha;   catch, data.edge_alpha = DEFAULT_EDGE_ALPHA; end;
try, data.trans_type;   catch, data.trans_type = DEFAULT_TRANS_TYPE; end;

datafields = fieldnames(data);
for index=1:length(datafields)
   switch datafields{index}
   case {'outtransfile' 'intransfile' 'intrans' 'trans_type'...
         'colorval' 'face_alpha' 'edge_alpha'...
   },;
   otherwise, error(['ts_pointreg: unrecognized option: ''' datafields{index} '''' ]);
   end;
end;

if ~strcmp(data.trans_type,'head2mri') & ~strcmp(data.trans_type,'mri2head')
  fprintf('%s: trans_type must be head2mri or mri2head\n',mfilename);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare figure and axes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figh = figure('Color',DEFAULT_FIG_COLOR, 'name', mfilename, ...
   'numbertitle', 'off', 'visible', 'on', 'toolbar','figure');

ax = axes('tag','UI_plot','parent',figh,...
   'Position',DEFAULT_AXES_POSITION); 
axis off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up uicontrols
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lower left corners of uicontrol groups
ui_transx      = [0.01,   0.02];
ui_transy      = [0.155,  0.02];
ui_transz      = [0.30,   0.02];
ui_rotx        = [0.51,   0.02];
ui_roty        = [0.655,  0.02];
ui_rotz        = [0.80,   0.02];

% define commands
trans_cmd = @translate_points;
rot_cmd = @rotate_points';

% group labels

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.06      0.13     0.05   ]; n=n+1; % translate
uipos(n,:)  = [ 0.51       0.06      0.13     0.05   ]; n=n+1; % rotate

% create uicontrols
n=1;
% transx +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',DEFAULT_FIG_COLOR,'Style','text', ...
 'Position',uipos(n,:), ...
 'Tag','UI_trans_label','string','Translation');
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',DEFAULT_FIG_COLOR,'Style','text', ...
 'Position',uipos(n,:), ...
 'Tag','UI_rot_label','string','Rotation');
n=n+1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transx group
ui_offset = ui_transx;

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.00      0.04     0.05   ]; n=n+1; % trans -
uipos(n,:)  = [ 0.102      0.00      0.04     0.05   ]; n=n+1; % trans +
uipos(n,:)  = [ 0.051      0.001     0.05     0.05   ]; n=n+1; % trans txt
uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

% create uicontrols
n=1;
% transx +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_transx_neg','string','-X',...
 'TooltipString','negative X translation',...
 'callback', trans_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_transx_pos','string','+X',...
 'TooltipString','positive X translation',...
 'callback', trans_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',[1 1 1], ...
 'Position',uipos(n,:), ...
 'Style','edit', ...
 'Tag','UI_transx_txt',...
 'TooltipString','specify X translation',...
 'string', DEFAULT_TRANSLATION);
n=n+1;

clear uipos u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transy group
ui_offset = ui_transy;

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.00      0.04     0.05   ]; n=n+1; % trans -
uipos(n,:)  = [ 0.102      0.00      0.04     0.05   ]; n=n+1; % trans +
uipos(n,:)  = [ 0.051      0.001     0.05     0.05   ]; n=n+1; % trans txt
uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

% create uicontrols
n=1;
% transy +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_transy_neg','string','-Y',...
 'TooltipString','negative Y translation',...
 'callback', trans_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_transy_pos','string','+Y',...
 'TooltipString','positive Y translation',...
 'callback', trans_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',[1 1 1], ...
 'Position',uipos(n,:), ...
 'Style','edit', ...
 'Tag','UI_transy_txt',...
 'TooltipString','specify Y translation',...
 'string', DEFAULT_TRANSLATION);
n=n+1;

clear uipos u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% transz group
ui_offset = ui_transz;

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.00      0.04     0.05   ]; n=n+1; % trans -
uipos(n,:)  = [ 0.102      0.00      0.04     0.05   ]; n=n+1; % trans +
uipos(n,:)  = [ 0.051      0.001     0.05     0.05   ]; n=n+1; % trans txt
uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

% create uicontrols
n=1;
% transz +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_transz_neg','string','-Z',...
 'TooltipString','negative Z translation',...
 'callback', trans_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_transz_pos','string','+Z',...
 'TooltipString','positive Z translation',...
 'callback', trans_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',[1 1 1], ...
 'Position',uipos(n,:), ...
 'Style','edit', ...
 'Tag','UI_transz_txt',...
 'TooltipString','specify Z translation',...
 'string', DEFAULT_TRANSLATION);
n=n+1;

clear uipos u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotx group
ui_offset = ui_rotx;

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.00      0.04     0.05   ]; n=n+1; % trans -
uipos(n,:)  = [ 0.102      0.00      0.04     0.05   ]; n=n+1; % trans +
uipos(n,:)  = [ 0.051      0.001     0.05     0.05   ]; n=n+1; % trans txt
uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

% create uicontrols
n=1;
% rotx +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_rotx_neg','string','-X',...
 'TooltipString','negative X rotation',...
 'callback', rot_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_rotx_pos','string','+X',...
 'TooltipString','positive X rotation',...
 'callback', rot_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',[1 1 1], ...
 'Position',uipos(n,:), ...
 'Style','edit', ...
 'Tag','UI_rotx_txt',...
 'TooltipString','specify X rotation',...
 'string', DEFAULT_ROTATION);
n=n+1;

clear uipos u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% roty group
ui_offset = ui_roty;

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.00      0.04     0.05   ]; n=n+1; % trans -
uipos(n,:)  = [ 0.102      0.00      0.04     0.05   ]; n=n+1; % trans +
uipos(n,:)  = [ 0.051      0.001     0.05     0.05   ]; n=n+1; % trans txt
uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

% create uicontrols
n=1;
% roty +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_roty_neg','string','-Y',...
 'TooltipString','negative Y rotation',...
 'callback', rot_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_roty_pos','string','+Y',...
 'TooltipString','positive Y rotation',...
 'callback', rot_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',[1 1 1], ...
 'Position',uipos(n,:), ...
 'Style','edit', ...
 'Tag','UI_roty_txt',...
 'TooltipString','specify Y rotation',...
 'string', DEFAULT_ROTATION);
n=n+1;

clear uipos u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rotz group
ui_offset = ui_rotz;

% set positions
n=1;
uipos(n,:)  = [ 0.01       0.00      0.04     0.05   ]; n=n+1; % trans -
uipos(n,:)  = [ 0.102      0.00      0.04     0.05   ]; n=n+1; % trans +
uipos(n,:)  = [ 0.051      0.001     0.05     0.05   ]; n=n+1; % trans txt
uipos(1:n-1,1) = uipos(1:n-1,1) + ui_offset(1);
uipos(1:n-1,2) = uipos(1:n-1,2) + ui_offset(2);

% create uicontrols
n=1;
% rotz +- buttons, text box
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_rotz_neg','string','-Z',...
 'TooltipString','negative Z rotation',...
 'callback', rot_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'Position',uipos(n,:), ...
 'Tag','UI_rotz_pos','string','+Z',...
 'TooltipString','positive Z rotation',...
 'callback', rot_cmd);
n=n+1;
u(n) = uicontrol('Parent',figh,'Units', 'normalized',...
 'BackgroundColor',[1 1 1], ...
 'Position',uipos(n,:), ...
 'Style','edit', ...
 'Tag','UI_rotz_txt',...
 'TooltipString','specify Z rotation',...
 'string', DEFAULT_ROTATION);
n=n+1;

clear uipos u;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get points and apply coordinate transform
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load trans file
if ~isempty(data.intransfile)
  if ~exist(data.intransfile,'file')
    fprintf('%s: input transfile not found, using identity matrix\n',...
      mfilename);
    data.intransfile = [];
  end;
end;

if ~isempty(data.intransfile)
  fprintf('%s: reading transfile %s\n',mfilename,data.intransfile);
  data.T = ts_read_transfile(data.intransfile,1000);
else
  fprintf('%s: using intrans\n',mfilename);
  data.T = data.intrans;
  data.T(1,4) = data.T(1,4)*1000;
  data.T(2,4) = data.T(2,4)*1000;
  data.T(3,4) = data.T(3,4)*1000;
end;
% internal trans_type is always head2mri (applied to head points)
if strcmp(data.trans_type,'mri2head')
  data.T = inv(data.T);
end;

% load hpts file
points = [];
fid=fopen(hptsfile,'rt');
if fid==-1
  error('failed to open file %s',hptsfile);
end;
p = 1;
while (~feof(fid))
  temp=fgetl(fid);
  if strcmp(temp(1),'#'), continue; end;
  spaces=regexp(temp,' +');
  strings = [];
  k = 1;
  for j=1:length(spaces)
    strings{j} = temp(k:spaces(j)-1);
    k = spaces(j) + 1;
  end;
  strings{j+1} = temp(k:length(temp));
  points(p).name =  strings{1};
  points(p).num = str2num(strings{2});
  % convert points to millimeters since surface is in those units
  points(p).x = str2double(strings{3})*1000;
  points(p).y = str2double(strings{4})*1000;
  points(p).z = str2double(strings{5})*1000;
  p=p+1;
end
fclose(fid);

% reformat point coords
data.numpoints = length(points);
data.coords = [];
data.coords(:,1) = cell2mat({points.x});
data.coords(:,2) = cell2mat({points.y});
data.coords(:,3) = cell2mat({points.z});
data.coords(:,4) = ones(data.numpoints,1);
data.pointcolor = DEFAULT_POINT_COLOR;
data.card_pointcolor = DEFAULT_CARD_POINT_COLOR;

% load tri file
data.surf = fs_read_trisurf(surffile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot head model and points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf, 'userdata', data);
redraw;
view([DEFAULT_VIEW_AZ,DEFAULT_VIEW_EL]);

return;

function translate_points(hobj,ed)
  transx = 0;
  transy = 0;
  transz = 0;
  data = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  if strcmp(tmp_UI,'UI_transx_neg')
    transx = -str2num(get(findobj(gcbf,'tag','UI_transx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transx_pos')
    transx = str2num(get(findobj(gcbf,'tag','UI_transx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transy_neg')
    transy = -str2num(get(findobj(gcbf,'tag','UI_transy_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transy_pos')
    transy = str2num(get(findobj(gcbf,'tag','UI_transy_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transz_neg')
    transz = -str2num(get(findobj(gcbf,'tag','UI_transz_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_transz_pos')
    transz = str2num(get(findobj(gcbf,'tag','UI_transz_txt'),'string'));
  end;
  tmp_T = eye(4);
  tmp_T(1,4) = transx;
  tmp_T(2,4) = transy;
  tmp_T(3,4) = transz;
  data.T = tmp_T*data.T;
  set(gcbf,'userdata',data);
  tmp_T = data.T;
  if strcmp(data.trans_type,'mri2head')
    tmp_T = inv(tmp_T);
  end;
  ts_write_transfile(data.outtransfile,tmp_T,0.001);
  redraw;
return;


function rotate_points(hobj,ed)
  rotx = 0;
  roty = 0;
  rotz = 0;
  data = get(gcbf, 'userdata');
  tmp_UI = get(gcbo,'tag');
  if strcmp(tmp_UI,'UI_rotx_neg')
    rotx = -str2num(get(findobj(gcbf,'tag','UI_rotx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_rotx_pos')
    rotx = str2num(get(findobj(gcbf,'tag','UI_rotx_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_roty_neg')
    roty = -str2num(get(findobj(gcbf,'tag','UI_roty_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_roty_pos')
    roty = str2num(get(findobj(gcbf,'tag','UI_roty_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_rotz_neg')
    rotz = -str2num(get(findobj(gcbf,'tag','UI_rotz_txt'),'string'));
  elseif strcmp(tmp_UI,'UI_rotz_pos')
    rotz = str2num(get(findobj(gcbf,'tag','UI_rotz_txt'),'string'));
  end;
  tmp_T = eye(4);
  if (rotx)
    tmp_T(2,2) = cos(rotx*pi/180);
    tmp_T(2,3) = sin(rotx*pi/180);
    tmp_T(3,2) = -sin(rotx*pi/180);
    tmp_T(3,3) = cos(rotx*pi/180);
  elseif (roty)
    tmp_T(1,1) = cos(roty*pi/180);
    tmp_T(1,3) = -sin(roty*pi/180);
    tmp_T(3,1) = sin(roty*pi/180);
    tmp_T(3,3) = cos(roty*pi/180);
  elseif (rotz)  
    tmp_T(1,1) = cos(rotz*pi/180);
    tmp_T(1,2) = sin(rotz*pi/180);
    tmp_T(2,1) = -sin(rotz*pi/180);
    tmp_T(2,2) = cos(rotz*pi/180);
  end;
  data.T = tmp_T*data.T;
  set(gcbf,'userdata',data);
  tmp_T = data.T;
  if strcmp(data.trans_type,'mri2head')
    tmp_T = inv(tmp_T);
  end;
  ts_write_transfile(data.outtransfile,tmp_T,0.001);
  redraw;
return;

function redraw()
  data = get(gcf, 'userdata');

  % plot surf
  cla;
  hold on;

  colors = data.colorval*ones(data.surf.nverts,1);
  h=trisurf(data.surf.faces,...
            data.surf.coords(:,1),...
            data.surf.coords(:,2),...
            data.surf.coords(:,3),...
            colors);
  colormap hsv;
  caxis([0 100]);

  set(h,'EdgeAlpha',data.edge_alpha,'FaceAlpha',data.face_alpha);

  % plot points
  tmp_coords = (data.T*data.coords')';
  if (data.numpoints<3)  
    plot3(tmp_coords(1:data.numpoints,1),tmp_coords(1:data.numpoints,2),...
      tmp_coords(1:data.numpoints,3),...
      [data.card_pointcolor '.'],'MarkerSize',16);
  else
    plot3(tmp_coords(1:3,1),tmp_coords(1:3,2),tmp_coords(1:3,3),...
      [data.card_pointcolor '.'],'MarkerSize',16);
  end;
  if (data.numpoints>3)
    plot3(tmp_coords(4:end,1),tmp_coords(4:end,2),tmp_coords(4:end,3),...
      [data.pointcolor '.'],'MarkerSize',16);
  end;

  axis image;
  axis off;
  xlabel('x');
  ylabel('y');
  zlabel('z');
return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

