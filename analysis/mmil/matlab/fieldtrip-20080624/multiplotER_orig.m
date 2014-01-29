function multiplotER(cfg, varargin)

% multiplotER plots the event-related fields or potentials versus time
% or of oscillatory activity (power or coherence) versus frequency. Multiple
% datasets can be overlayed.  The plots are arranged according to their
% location specified in the layout.
%
% Use as:
%   multiplotER(cfg, data)
%   multiplotER(cfg, data, data2, ..., dataN)
%
% The data can be an ERP/ERF produced by TIMELOCKANALYSIS, a powerspectrum 
% produced by FREQANALYSIS or a coherencespectrum produced by FREQDESCRIPTIVES. 
% If you specify multiple datasets they must contain the same channels, etc.
%
% The configuration can have the following parameters:
% cfg.xparam        = field to be plotted on x-axis (default depends on data.dimord)
%                     'time' or 'freq' 
% cfg.zparam        = field to be plotted on y-axis (default depends on data.dimord)
%                     'avg', 'powspctrm' or 'cohspctrm' 
% cfg.maskparameter = field in the first dataset to be used for marking significant data
% cfg.xlim          = 'maxmin' or [xmin xmax] (default = 'maxmin')
% cfg.ylim          = 'maxmin' or [ymin ymax] (default = 'maxmin')
% cfg.cohrefchannel = Name of reference-channel, only for visualizing coherence 
% cfg.baseline      = 'yes','no' or [time1 time2] (default = 'no'), see TIMELOCKBASELINE or FREQBASELINE
% cfg.baselinetype  = 'absolute' or 'relative' (default = 'absolute')
% cfg.trials        = 'all' or a selection given as a 1xN vector (default = 'all')
% cfg.axes          = 'yes', 'no' (default = 'yes')
%                     Draw x- and y-axes for each graph
% cfg.box           = 'yes', 'no' (default = 'no')
%                     Draw a box around each graph
% cfg.comment       = string of text (default = date + colors)
%                     Add 'comment' to graph (according to COMNT in the layout)
% cfg.showlabels    = 'yes', 'no' (default = 'no')
% cfg.fontsize      = font size of comment and labels (if present) (default = 8)
% cfg.interactive   = Interactive plot 'yes' or 'no' (default = 'no')
%                     In a interactive plot you can select areas and produce a new
%                     interactive plot when a selected area is clicked. Multiple areas 
%                     can be selected by holding down the SHIFT key.
% cfg.renderer      = 'painters', 'zbuffer',' opengl' or 'none' (default = 'opengl')
%
% cfg.layout        = specify the channel layout for plotting using one of 
%                     the following ways:
%
% The layout defines how the channels are arranged and what the size of each
% subplot is. You can specify the layout in a variety of ways:
%  - you can provide a pre-computed layout structure (see prepare_layout)
%  - you can give the name of an ascii layout file with extension *.lay
%  - you can give the name of an electrode file
%  - you can give an electrode definition, i.e. "elec" structure
%  - you can give a gradiometer definition, i.e. "grad" structure
% If you do not specify any of these and the data structure contains an
% electrode or gradiometer structure, that will be used for creating a
% layout. If you want to have more fine-grained control over the layout
% of the subplots, you should create your own layout file.
%
% See also:
%   multiplotTFR, singleplotER, singleplotTFR, topoplotER, topoplotTFR,
%   prepare_layout.

% Undocumented local options:
% cfg.graphcolor
% cfg.layoutname
%
% This function depends on TIMELOCKBASELINE which has the following options:
% cfg.baseline, documented
% cfg.channel
% cfg.blcwindow
% cfg.previous
% cfg.version
%
% This function depends on FREQBASELINE which has the following options:
% cfg.baseline, documented
% cfg.baselinetype

% Copyright (C) 2003-2006, Ole Jensen
%
% $Log: multiplotER.m,v $
% Revision 1.38  2008/01/29 19:43:33  sashae
% added option for trial selection; plot functions now also accept data with
% repetitions (either trials or subjects), the avg is computed and plotted
% removed some old code
%
% Revision 1.37  2007/07/05 08:36:19  roboos
% fixed bug related to layout and varargin
%
% Revision 1.36  2007/06/19 13:57:07  ingnie
% axis wider if cfg.box = 'yes', changed color box and mask-patch
%
% Revision 1.35  2007/06/07 14:32:09  ingnie
% fixed error: baselining should depend on xparam not yparam
%
% Revision 1.34  2007/06/05 16:14:23  ingnie
% added cfg.maskparameter
%
% Revision 1.33  2007/04/25 17:23:46  ingnie
% *** empty log message ***
%
% Revision 1.32  2007/04/19 10:26:11  roboos
% Added a warning to the "Apply baseline correction" sections. If a user doesn't set yparam, baselining is not applied. (thanks to Doug)
%
% Revision 1.31  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.30  2007/03/27 11:05:19  ingnie
% changed call to fixdimord in call to checkinput
%
% Revision 1.29  2007/03/21 15:56:18  chrhes
% updated documentation regarding the fact that cfg.layout can also contain a
% layout structure obtained using the function prepare_layout.m
%
% Revision 1.28  2007/03/14 08:43:12  roboos
% replaced call to createlayout to prepare_layout, made some small changes to 
% the lay structure
%
% Revision 1.27  2007/01/09 10:46:31  roboos
% fixed accidental typo in function definition
%
% Revision 1.26  2007/01/09 10:41:43  roboos
% Added the option cfg.renderer, default is opengl. Interactive plotting
% on linux/VNC sometimes does not work, using cfg.renderer='painters'
% seems to fix it.
%
% Revision 1.25  2006/07/27 15:33:29  roboos
% autodetect params for timelocked data with keeptrials
%
% Revision 1.24  2006/06/19 11:11:37  roboos
% fixed small bug in the conversion of coherence data: first select labels for 
% the channels, then for the channelcombinations
%
% Revision 1.23  2006/05/30 14:19:18  ingnie
% if axis is 'yes' axis are always drawn and not only when data crosses zero,
% updated documentation
%
% Revision 1.22  2006/05/26 12:47:28  ingnie
% added error when labels in layout and labels in data do not match and therefore
% no data is selected to be plotted
%
% Revision 1.21  2006/05/11 14:53:52  ingnie
% when cfg.ylim is 'maxmin', the scaling (ymin ymax) only determined by channels
% that are present in the layout.
%
% Revision 1.20  2006/05/03 08:12:51  ingnie
% updated documentation
%
% Revision 1.19  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.18  2006/03/17 14:47:27  denpas
% Updated documentation.
%
% Revision 1.17  2006/03/17 14:43:59  denpas
% Fixed cfg.yparam / cfg.zparam bug. Either param may now be used to specify
% the y-axis in case of 2D data.
%
% Revision 1.16  2006/03/13 14:08:47  denpas
% Fixed a bug concerning dynamic fieldnames in combination with (:), which doesn't work.
% Removed the (:) and added min(min(...)) and max(max(...)).
%
% Revision 1.15  2006/03/02 13:54:47  jansch
% fixed multiple small bugs
%
% Revision 1.14  2006/02/28 12:43:15  roboos
% made plotting of coherence consistent between all xxxplotER functions
% made baselining consistent, use cfg.xparam to decide between freqbaseline 
% and timelockbaseline
%
% Revision 1.13  2006/02/27 15:03:03  denpas
% many changes, most important is added interactive functionality
% made data selection consistent between different plot functions
% changed dimord for consistency
%
% Revision 1.11  2005/08/18 12:15:41  jansch
% added support to plot subfields
%
% Revision 1.10  2005/05/17 17:50:37  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer 
% correspondence with Matlab documentation on shortcircuited evaluation of 
% sequential boolean constructs
%
% Revision 1.9  2005/04/29 12:46:44  roboos
% small change in help
%
% Revision 1.8  2005/04/29 12:40:43  roboos
% cleaned up and updated the help
%
% Revision 1.7  2005/04/12 10:06:57  olejen
% clf added
%
% Revision 1.6  2005/04/06 07:40:56  jansch
% included option cfg.cohrefchannel. updated help.
%
% Revision 1.5  2005/02/07 17:12:00  roboos
% changed handling of layout files (using new function createlayout), now also 
% supports automatic layout creation based on gradiometer/electrode definition 
% in data, updated help, cleaned up indentation
%
% Revision 1.4  2005/01/27 09:31:49  roboos
% applied autoindentation on code, removed many empty lines and spaces,
% replaced layoutfile reading with read_lay, applied doudavs code to all input arguments
% implemented automatic detection of arguments to plot,
% updated help, changed input from p1, p2, p3... to varargin
%
% Revision 1.3  2004/09/24 15:54:54  roboos
% included the suggested improvements by Doug Davidson: added option cfg.cohtargetchannel
% and updated the help
%
% Revision 1.2  2004/09/01 17:59:28  roboos
% added copyright statements to all filed
% added cfg.version to all functions that give configuration in their output
% added cfg.previous to all functions with input data containing configuration details

clf

% For backward compatibility with old data structures:
for i=1:length(varargin)
  varargin{i} = checkdata(varargin{i});
end

% set the defaults:
if ~isfield(cfg,'baseline'),    cfg.baseline    = 'no';                        end
if ~isfield(cfg,'trials'),      cfg.trials      = 'all';                       end
if ~isfield(cfg,'xlim'),        cfg.xlim        = 'maxmin';                    end
if ~isfield(cfg,'ylim'),        cfg.ylim        = 'maxmin';                    end
if ~isfield(cfg,'comment'),     cfg.comment     = strcat([date '\n']);         end
if ~isfield(cfg,'axes'),        cfg.axes        = 'yes';                       end
if ~isfield(cfg,'showlabels'),  cfg.showlabels  = 'no';                        end
if ~isfield(cfg,'box'),         cfg.box         = 'no';                        end
if ~isfield(cfg,'fontsize'),    cfg.fontsize    = 8;                           end
if ~isfield(cfg,'graphcolor')   cfg.graphcolor  = ['brgkywrgbkywrgbkywrgbkyw'];end
if ~isfield(cfg,'interactive'), cfg.interactive = 'no';                        end
if ~isfield(cfg,'renderer'),    cfg.renderer    = 'opengl';                    end
if ~isfield(cfg,'maskparameter'),cfg.maskparameter = [];                       end

GRAPHCOLOR = ['k' cfg.graphcolor ];

% Set x/y/zparam defaults according to varargin{1}.dimord value:
if strcmp(varargin{1}.dimord, 'chan_time')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';                   end
elseif strcmp(varargin{1}.dimord, 'chan_freq')
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
elseif strcmp(varargin{1}.dimord, 'subj_chan_time') || strcmp(varargin{1}.dimord, 'rpt_chan_time')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  for i=1:(nargin-1)
    varargin{i} = timelockanalysis(tmpcfg, varargin{i});
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='time';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='avg';                   end
elseif strcmp(varargin{1}.dimord, 'subj_chan_freq') || strcmp(varargin{1}.dimord, 'rpt_chan_freq')
  tmpcfg = [];
  tmpcfg.trials = cfg.trials;
  tmpcfg.jackknife = 'no';
  for i=1:(nargin-1)
    if isfield(varargin{i}, 'crsspctrm'), varargin{i} = rmfield(varargin{i}, 'crsspctrm'); end % on the fly computation of coherence spectrum is not supported
    varargin{i} = freqdescriptives(tmpcfg, varargin{i});
  end
  if ~isfield(cfg, 'xparam'),      cfg.xparam='freq';                  end
  if ~isfield(cfg, 'zparam'),      cfg.zparam='powspctrm';             end
end

% Make sure cfg.yparam and cfg.zparam become equivalent if only one is defined:
if (isfield(cfg, 'yparam')) & (~isfield(cfg, 'zparam'))
  cfg.zparam = cfg.yparam;
elseif (~isfield(cfg, 'yparam')) & (isfield(cfg, 'zparam'))
  cfg.yparam = cfg.zparam;
end

% Old style coherence plotting with cohtargetchannel is no longer supported:
if isfield(cfg,'cohtargetchannel'), 
  error('cfg.cohtargetchannel is obsolete, check the documentation for help about coherence plotting.');
end

for k=1:length(varargin)
  % Check for unconverted coherence spectrum data:
  if (strcmp(cfg.zparam,'cohspctrm')) & (isfield(varargin{k}, 'labelcmb'))
    % A reference channel is required:
    if ~isfield(cfg,'cohrefchannel'),
      error('no reference channel specified');
    end
    % Convert 2-dimensional channel matrix to a single dimension:
    sel1                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,2));
    sel2                  = strmatch(cfg.cohrefchannel, varargin{k}.labelcmb(:,1));
    fprintf('selected %d channels for coherence\n', length(sel1)+length(sel2));
    varargin{k}.cohspctrm = varargin{k}.cohspctrm([sel1;sel2],:,:);
    varargin{k}.label     = [varargin{k}.labelcmb(sel1,1);varargin{k}.labelcmb(sel2,2)];
    varargin{k}.labelcmb  = varargin{k}.labelcmb([sel1;sel2],:);
    varargin{k}           = rmfield(varargin{k}, 'labelcmb');
  end

  % Apply baseline correction:
  if ~strcmp(cfg.baseline, 'no')
    if strcmp(cfg.xparam, 'time')
      varargin{k} = timelockbaseline(cfg, varargin{k});
    elseif strcmp(cfg.xparam, 'freq')
      varargin{k} = freqbaseline(cfg, varargin{k});
    else 
      warning('Baseline not applied, please set cfg.xparam');
    end
  end
end

% Get physical x-axis range:
if strcmp(cfg.xlim,'maxmin')
  % Find maxmin throughout all varargins:
  xmin = [];
  xmax = [];
  for i=1:length(varargin)
    xmin = min([xmin varargin{i}.(cfg.xparam)]);
    xmax = max([xmax varargin{i}.(cfg.xparam)]);
  end
else
  xmin = cfg.xlim(1);
  xmax = cfg.xlim(2);
end

% Find corresponding x-axis bins:
xidc = find(varargin{1}.(cfg.xparam) >= xmin & varargin{1}.(cfg.xparam) <= xmax);

% Align physical x-axis range to the array bins:
xmin = varargin{1}.(cfg.xparam)(xidc(1));
xmax = varargin{1}.(cfg.xparam)(xidc(end));

% Read or create the layout that will be used for plotting:
lay = prepare_layout(cfg, varargin{1});

% Get physical y-axis range (ylim / zparam):
if strcmp(cfg.ylim,'maxmin')
  % Find maxmin throughout all varargins:
  ymin = [];
  ymax = [];
  for i=1:length(varargin)
    % Select the channels in the data that match with the layout:
    dat = [];  
    dat = getsubfield(varargin{i}, cfg.zparam);
    [seldat, sellay] = match_str(varargin{k}.label, lay.label);
    if isempty(seldat)
      error('labels in data and labels in layout do not match'); 
    end
    data = dat(seldat,:);
    ymin = min([ymin min(min(min(data)))]);
    ymax = max([ymax max(max(max(data)))]);
  end
else
  ymin = cfg.ylim(1);
  ymax = cfg.ylim(2);
end

% convert the layout to Ole's style of variable names
X      = lay.pos(:,1);
Y      = lay.pos(:,2);
Width  = lay.width;
Height = lay.height;
Lbl    = lay.label;

% Create empty channel coordinates and labels arrays:
chanX(1:length(Lbl)) = NaN;
chanY(1:length(Lbl)) = NaN;
chanLabels = cell(1,length(Lbl));

hold on;
colorLabels = [];

% Plot each data set:
for k=1:length(varargin)
  P          = getsubfield(varargin{k}, cfg.zparam);
  Labels     = getfield(varargin{k}, 'label');

  if length(varargin) > 1
    colorLabels = [colorLabels inputname(k+1) '=' GRAPHCOLOR(k+1) '\n'];
  end

  style = GRAPHCOLOR(k+1);
  
  for m=1:length(Lbl)
    l = cellstrmatch(Lbl(m),Labels);
    if ~isempty(l)
      if ~isempty(cfg.maskparameter)
        mask = varargin{1}.(cfg.maskparameter)(l,:);
      else
        mask = [];
      end
      % Plot ER:  
      plotWnd(varargin{k}.(cfg.xparam),P(l,:),xidc,[xmin xmax],[ymin ymax], ...
        X(m), ...
        Y(m), ...
        Width(m), ...
        Height(m), ...
        Lbl(m), ...
        cfg,style, mask);

      % Keep ER plot coordinates (at centre of ER plot), and channel labels (will be stored in the figure's UserData struct):
      chanX(m) = X(m) + 0.5 * Width(m);
      chanY(m) = Y(m) + 0.5 * Height(m);
      chanLabels{m} = Lbl{m};
    end
  end
end

% Add the colors of the different datasets to the comment:
cfg.comment = [cfg.comment colorLabels];

% Write comment text:
l = cellstrmatch('COMNT',Lbl);
if ~isempty(l)
  text(X(l),Y(l),sprintf(cfg.comment),'Fontsize',cfg.fontsize);
end

% Plot scales:
l = cellstrmatch('SCALE',Lbl);
if ~isempty(l)
  plotScales([xmin xmax],[ymin ymax],X(l),Y(l),Width(1),Height(1),cfg)
end

% Make the figure interactive:
if strcmp(cfg.interactive, 'yes')
  userData.hFigure = gcf;
  userData.hAxes = gca;
  for i=1:10
    userData.hSelection{i} = plot(mean(X), mean(Y));
    set(userData.hSelection{i}, 'XData', [mean(X)]);
    set(userData.hSelection{i}, 'YData', [mean(Y)]);
    set(userData.hSelection{i}, 'Color', [0 0 0]);
    set(userData.hSelection{i}, 'EraseMode', 'xor');
    set(userData.hSelection{i}, 'LineStyle', '--');
    set(userData.hSelection{i}, 'LineWidth', 1.5);
    set(userData.hSelection{i}, 'Visible', 'on');
    userData.range{i} = [];
  end
  userData.iSelection = 0;
  userData.plotType = 'multiplot';
  userData.selecting = 0;
  userData.selectionType = '';
  userData.selectAxes = 'z';
  userData.lastClick = [];
  userData.cfg = cfg;
  userData.data = varargin;
  userData.chanX = chanX;
  userData.chanY = chanY;
  userData.chanLabels = chanLabels;
  tag = sprintf('%.5f', 10000 * rand(1));
  set(gcf, 'Renderer', cfg.renderer);
  set(gcf, 'Tag', tag);
  set(gcf, 'UserData', userData);
  set(gcf, 'WindowButtonMotionFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 0);']);
  set(gcf, 'WindowButtonDownFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 1);']);
  set(gcf, 'WindowButtonUpFcn', ['plotSelection(get(findobj(''Tag'', ''' tag '''), ''UserData''), 2);']);
end

axis tight
axis off
if strcmp(cfg.box, 'yes')
  abc = axis;
  axis(abc + [-1 +1 -1 +1]*mean(abs(abc))/10)
end
orient landscape
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotScales(xlim,ylim,xpos,ypos,width,height,cfg)
x1 =  xpos;
x2 =  xpos+width;
y1 =  ypos;
y2 =  ypos+width;
plot([xpos xpos+width xpos+width xpos xpos],[ypos ypos ypos+height ypos+height ypos],'b');
if xlim(1) <=  0 && xlim(2) >= 0
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  plot(xs,ys,'b');
end

if ylim(1) <= 0 && ylim(2) >= 0
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  plot(xs,ys,'b');
end

text( x1,y1,num2str(xlim(1),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
text( x2,y1,num2str(xlim(2),3),'rotation',90,'HorizontalAlignment','Right','VerticalAlignment','middle','Fontsize',cfg.fontsize);
text( x2,y1,num2str(ylim(1),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);
text( x2,y2,num2str(ylim(2),3),'HorizontalAlignment','Left','VerticalAlignment','bottom','Fontsize',cfg.fontsize);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotWnd(x,y,xidc,xlim,ylim,xpos,ypos,width,height,label,cfg,style,mask)
set(gca,'FontSize',cfg.fontsize);

x = x(xidc);
y = y(xidc);

% Clip out of bounds y values:
y(find(y > ylim(2))) = ylim(2);
y(find(y < ylim(1))) = ylim(1);

xs = xpos+width*(x-xlim(1))/(xlim(2)-xlim(1));
ys = ypos+height*(y-ylim(1))/(ylim(2)-ylim(1));
plot(xs,ys,style)

if strcmp(cfg.showlabels,'yes')
  text(xpos,ypos+1.0*height,label,'Fontsize',cfg.fontsize)
end

% Draw axes:
if strcmp(cfg.axes,'yes')
  % Draw y axis
  xs =  xpos+width*([0 0]-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*(ylim-ylim(1))/(ylim(2)-ylim(1));
  plot(xs,ys,'k');
  % Draw x axis
  xs =  xpos+width*(xlim-xlim(1))/(xlim(2)-xlim(1));
  ys =  ypos+height*([0 0]-ylim(1))/(ylim(2)-ylim(1));
  plot(xs,ys,'k');
end

% Draw box around plot:
if strcmp(cfg.box,'yes')
  plot([xpos xpos+width xpos+width xpos xpos],[ypos ypos ypos+height ypos+height ypos],'k');
end

% Add mask patch
if ~isempty(mask)
  % determine how many boxes
  foundbeg = 0;
  foundend = 0;
  beg  = [];
  eind = [];
  for i = 1:length(mask)
    if ~foundbeg  && mask(i) == 1 
      beg(length(beg)+1) = i;
      foundbeg = 1;
      foundend = 0;
    elseif ~foundbeg  && mask(i) == 0
      %next
    elseif ~foundend  && mask(i) == 1
      %next
    elseif ~foundend  && mask(i) == 0
      eind(length(eind)+1) = i-1;
      foundend = 1;
      foundbeg = 0;
    end
  end
  if length(eind) == length(beg)-1
    eind(length(eind)+1) = length(mask);
  end
  numbox = length(beg);
  for i = 1:numbox
    xmaskmin = xpos+width*(x(beg(i))-xlim(1))/(xlim(2)-xlim(1));
    xmaskmax = xpos+width*(x(eind(i))-xlim(1))/(xlim(2)-xlim(1));
    %plot([xmaskmin xmaskmax xmaskmax xmaskmin xmaskmin],[ypos ypos ypos+height ypos+height ypos],'r');
    hs = patch([xmaskmin xmaskmax xmaskmax xmaskmin xmaskmin],[ypos ypos ypos+height ypos+height ypos], [.6 .6 .6]);
    set(hs, 'EdgeColor', 'none');
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function l = cellstrmatch(str,strlist)
l = [];
for k=1:length(strlist)
  if strcmp(char(str),char(strlist(k)))
    l = [l k];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = pwrspctm2cohspctrm(cfg, freq)
% This creates a copy of freq with new entries in the freq.pwrspctrm. The
% TFR corresponding to cfg.cohtargetchannel is set to ones, but all other
% channels are replaced with the coherence TFRs from freq.cohspctrm.
%
% [z] = pwrspctrm2cohspctrm(cfg, freq);
%
% freq should be organised in a structure as obtained from
% the freqdescriptives function.
%
% doudav; 24.09.04
z = freq;
[a a1] = match_str( cfg.cohtargetchannel, z.label);
[a b1] = match_str( cfg.cohtargetchannel, z.labelcmb(:,1));
[a b2] = match_str( cfg.cohtargetchannel, z.labelcmb(:,2));
[a c1] = match_str( z.labelcmb(b1,2), z.label);
[a c2] = match_str( z.labelcmb(b2,1), z.label);
targets  = [c1; c2];
targets2 = [b1; b2];
z.powspctrm(a1,:,:) = ones(size(z.powspctrm(a1,:,:)));
for i=1:size(targets,1),
  z.powspctrm(targets(i),:,:) = squeeze(z.cohspctrm(targets2(i),:,:));
end
clear a a1 b1 b2 c1 c2;
