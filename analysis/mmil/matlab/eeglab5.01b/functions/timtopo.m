% timtopo()   - plot all channels of a data epoch on the same axis 
%               and map its scalp map(s) at selected latencies.
% Usage:
%  >> timtopo(data,'chan_locs');
%  >> timtopo(data,'chan_locs',[limits],[plottimes]','title',[plotchans], ...
%                                                     [voffsets], 'key', 'val', ...);
% Inputs:
%  data       = (channels,frames) single-epoch data matrix
%  chan_locs  = channel location file or EEG.chanlocs structure. 
%               See >> topoplot example for file format.
%
% Optional ordered inputs:
%  [limits]   = [minms maxms minval maxval] data limits for latency (in ms) and y-values
%                (assumes uV) {default|0 -> use [0 npts-1 data_min data_max]; 
%               else [minms maxms] or [minms maxms 0 0] -> use
%                [minms maxms data_min data_max]
%  plottimes  = [vector] latencies (in ms) at which to plot scalp maps 
%                {default|NaN -> latency of maximum variance}
% 'title'     = [string] plot title {default|0 -> none}
%  plotchans  = vector of data channel(s) to plot {default|0 -> all}
%  voffsets   = vector of (plotting-unit) distances vertical lines should extend 
%                above the data (in special cases) {default -> all = standard}
%
% Optional keyword, arg pair inputs (must come after the above):
% 'topokey','val' = optional topoplot() scalp map plotting arguments. See >> help topoplot 
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1-10-98 
%
% See also: envtopo(), topoplot()

% Copyright (C) 1-10-98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: timtopo.m,v $
% Revision 1.72  2005/05/16 20:27:46  scott
% adjusted several plotting details, oblique line height, cbar, text, etc.
%
% Revision 1.71  2005/03/05 02:28:36  arno
% debug argument passing to topoplot
%
% Revision 1.70  2004/10/20 15:20:27  scott
% still bug at line 551 -sm
%
% Revision 1.69  2004/10/20 15:19:03  scott
% nothing
%
% Revision 1.68  2004/04/25 20:07:28  scott
% typo
%
% Revision 1.67  2004/04/25 20:05:36  scott
% help message modernize. adjust text, add white vertical underlines, made limits input more flexible
%
% Revision 1.66  2003/11/26 18:20:44  scott
% Time -> Latency
%
% Revision 1.65  2003/08/06 00:27:10  arno
% remove postp (made matlab 5.3 crash
%
% Revision 1.64  2003/03/05 03:22:40  scott
% cleanup
%
% Revision 1.63  2003/03/05 03:19:23  scott
% topoleft
%
% Revision 1.62  2003/03/05 03:17:01  scott
% typo
%
% Revision 1.61  2003/03/05 03:16:41  scott
% typo
%
% Revision 1.60  2003/03/05 03:16:10  scott
% topoleft
%
% Revision 1.59  2003/03/05 03:06:37  scott
% cleanup
%
% Revision 1.58  2003/03/05 03:04:52  scott
% head_sep
%
% Revision 1.57  2003/03/05 03:03:43  scott
% topowidth
%
% Revision 1.56  2003/03/05 02:54:49  scott
% topoargs
%
% Revision 1.55  2003/03/05 02:53:26  scott
% same
%
% Revision 1.54  2003/03/05 02:53:00  scott
% topostring
%
% Revision 1.53  2003/03/05 02:46:01  scott
% topowdith
%
% Revision 1.52  2003/03/05 02:44:04  scott
% topowidth
% .,
%
% Revision 1.51  2003/03/05 02:41:49  scott
% topowidth
%
% Revision 1.50  2003/03/05 02:40:58  scott
% same
%
% Revision 1.49  2003/03/05 02:40:00  scott
% topowidth
%
% Revision 1.48  2003/03/05 02:36:31  scott
% same
%
% Revision 1.47  2003/03/05 02:36:01  scott
% eval topoplot
%
% Revision 1.46  2003/03/05 02:34:32  scott
% topoplot
%
% Revision 1.45  2003/03/05 02:33:01  scott
% topoplot
%
% Revision 1.44  2003/03/05 02:30:37  scott
% emarkersize
%
% Revision 1.43  2003/03/05 02:28:35  scott
% same
%
% Revision 1.42  2003/03/05 02:27:54  scott
% emarkersize
%
% Revision 1.41  2003/03/05 02:04:58  scott
% emarkersize
%
% Revision 1.40  2003/03/05 02:04:07  scott
% emarkersize -sm
%
% Revision 1.39  2003/03/05 01:56:34  scott
% topowidth -sm
%
% Revision 1.38  2003/03/05 01:55:20  scott
% topowidth -sm
%
% Revision 1.37  2003/03/05 01:51:57  scott
% topowidth -sm
%
% Revision 1.36  2003/03/04 21:17:52  scott
% axfont -sm
%
% Revision 1.35  2003/03/04 21:10:50  scott
% titlefont -sm
%
% Revision 1.34  2003/03/04 18:52:41  scott
% cleaning up -sm
%
% Revision 1.33  2003/03/04 18:48:25  scott
% test size -sm
%
% Revision 1.32  2003/03/04 18:44:58  scott
% title text debug -sm
%
% Revision 1.31  2003/03/04 18:43:45  scott
% title text -sm
%
% Revision 1.30  2003/03/04 18:41:09  scott
% title text -sm
%
% Revision 1.29  2003/03/04 18:40:34  scott
% title text -sm
%
% Revision 1.28  2003/03/04 18:38:57  scott
% title text
%
% Revision 1.27  2003/03/04 18:35:02  scott
% cbar text -sm
%
% Revision 1.26  2003/03/04 18:34:07  scott
% cbar text
%
% Revision 1.25  2003/03/04 18:32:59  scott
% cbar text -sm
%
% Revision 1.24  2003/03/04 18:28:38  scott
% cbar text -sm
%
% Revision 1.23  2003/03/04 18:24:36  scott
% final debugs? -sm
%
% Revision 1.22  2003/03/04 18:22:30  scott
% debug last -sm
%
% Revision 1.21  2003/03/04 18:20:01  scott
% debug last -sm
%
% Revision 1.20  2003/03/04 18:15:06  scott
% debug last -sm
%
% Revision 1.19  2003/03/04 18:13:58  scott
% debug last -sm
%
% Revision 1.18  2003/03/04 18:11:36  scott
% debug last -sm
%
% Revision 1.17  2003/03/04 18:09:45  scott
% debug last -sm
%
% Revision 1.16  2003/03/04 18:07:29  scott
% debug last -sm
%
% Revision 1.15  2003/03/04 18:06:21  scott
% using changeunits -sm
%
% Revision 1.14  2003/03/04 17:49:36  scott
% debug oblique lines -sm
%
% Revision 1.13  2003/03/04 17:41:55  scott
% debug last -sm
%
% Revision 1.12  2003/03/04 17:40:51  scott
% debug head size for sbplots -sm
%
% Revision 1.11  2003/03/04 17:37:47  scott
% debug last -sm
%
% Revision 1.10  2003/03/04 17:36:26  scott
% debug last -sm
%
% Revision 1.9  2003/03/04 17:21:43  scott
% debug last -sm
%
% Revision 1.8  2003/03/04 17:20:23  scott
% debug last -sm
%
% Revision 1.7  2003/03/04 17:16:01  scott
% debug last -sm
%
% Revision 1.6  2003/03/04 17:11:51  scott
% edit to work in subplot axes -sm
%
% Revision 1.5  2002/11/15 03:07:45  arno
% header for web
%
% Revision 1.4  2002/08/28 00:52:30  arno
% allow to plot NaN with other latencies
%
% Revision 1.3  2002/08/27 00:20:54  arno
% debugging colorbar->cbar (for menus)
%
% Revision 1.2  2002/08/12 23:48:05  arno
% debug absmax
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 5-31-00 added o-time line and possibility of plotting 1 channel -sm & mw
% 11-02-99 added maplimits arg -sm
% 01-22-01 added to help message -sm
% 01-25-02 reformated help & license, added link -ad 
% 03-15-02 add all topoplot options -ad

function M = timtopo(data,chan_locs,limits,plottimes,titl,plotchans,voffsets, varargin)

MAX_TOPOS = 24;

if nargin < 1
   help timtopo;
   return
end

[chans,frames] = size(data);
icadefs;   

if nargin < 7 | voffsets == 0
  voffsets = zeros(1,MAX_TOPOS);
end

if nargin < 6
   plotchans = 0;
end

if plotchans==0
   plotchans = 1:chans;
end

if nargin < 5,
   titl = '';     % DEFAULT NO TITLE
end

plottimes_set=1;   % flag variable
if nargin< 4 | isempty(plottimes) | any(isnan(plottimes))
   plottimes_set = 0;
end

limitset = 0;
if nargin < 3,
    limits = 0;
elseif length(limits)>1
    limitset = 1;
end

if nargin < 2
    chan_locs = 'chan.locs';  % DEFAULT CHAN_FILE
end
if isnumeric(chan_locs) & chan_locs == 0,
    chan_locs = 'chan.locs';  % DEFAULT CHAN_FILE
end

  %
  %%%%%%%%%%%%%%%%%%%%%%% Read and adjust limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % defaults: limits == 0 or [0 0 0 0]
  if ( length(limits) == 1 & limits==0) | (length(limits)==4 & ~any(limits))  
    xmin=0;
    xmax=frames-1;
    ymin=min(min(data));
    ymax=max(max(data));
  elseif length(limits) == 2  % [minms maxms] only
    ymin=min(min(data));
    ymax=max(max(data));
    xmin = limits(1);
    xmax = limits(2);
 elseif length(limits) == 4
    xmin = limits(1);
    xmax = limits(2);
    if any(limits([3 4]))
      ymin = limits(3);
      ymax = limits(4);
    else % both 0
      ymin=min(min(data));
      ymax=max(max(data));
    end
  else
    fprintf('timtopo(): limits format not correct. See >> help timtopo.\n');
    return
  end;

  if xmax == 0 & xmin == 0,
    x = (0:1:frames-1);
    xmin = 0;
    xmax = frames-1;
  else
    dx = (xmax-xmin)/(frames-1);
    x=xmin*ones(1,frames)+dx*(0:frames-1); % compute x-values
  end;
  if xmax<=xmin,
      fprintf('timtopo() - in limits, maxms must be > minms.\n')
      return
  end

  if ymax == 0 & ymin == 0,
      ymax=max(max(data));
      ymin=min(min(data));
  end
  if ymax<=ymin,
      fprintf('timtopo() - in limits, maxval must be > minmval.\n')
      return
  end

sampint = (xmax-xmin)/(frames-1); % sampling interval = 1000/srate;
x = xmin:sampint:xmax;   % make vector of x-values

%
%%%%%%%%%%%%%%% Compute plot times/frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if plottimes_set == 0
  [mx plotframes] = max(sum(data.*data)); 
                  % default plotting frame has max variance
  if nargin< 4 | isempty(plottimes)
	  plottimes = x(plotframes);
  else
	  plottimes(find(isnan(plottimes))) = x(plotframes);
  end;
  plottimes_set = 1;
end;

if plottimes_set == 1
  ntopos = length(plottimes);
  if ntopos > MAX_TOPOS
    fprintf('timtopo(): too many plottimes - only first %d will be shown!\n',MAX_TOPOS);
    plottimes = plottimes(1:MAX_TOPOS);
    ntopos = MAX_TOPOS;
  end

  if max(plottimes) > xmax | min(plottimes)< xmin
    fprintf(...
'timtopo(): at least one plottimes value outside of epoch latency range - cannot plot.\n');
    return
  end

  plottimes = sort(plottimes); % put map latencies in ascending order, 
                               % else map lines would cross.
  xshift = [x(2:frames) xmax];
  plotframes = ones(size(plottimes));
  for t = 1:ntopos
    time = plottimes(t);
    plotframes(t) = find(time>=x & time < xshift);
  end
end

vlen = length(voffsets); % extend voffsets if necessary
i=1;
while vlen< ntopos
        voffsets = [voffsets voffsets(i)];
        i=i+1;
        vlen=vlen+1;
end

%
%%%%%%%%%%%%%%%%  Compute title and axes font sizes %%%%%%%%%%%%%%%
%
pos = get(gca,'Position');
axis('off')
cla % clear the current axes
if pos(4)>0.70
   titlefont= 16;
   axfont = 16;
elseif pos(4)>0.40
   titlefont= 14;
   axfont = 14;
elseif pos(4)>0.30
   titlefont= 12;
   axfont = 12;
elseif pos(4)>0.22
   titlefont= 10;
   axfont = 10;
else
   titlefont= 8;
   axfont = 8;
end

%
%%%%%%%%%%%%%%%% Compute topoplot head width and separation %%%%%%%%%%%%%%%
%
head_sep = 0.2;
topowidth = pos(3)/((6*ntopos-1)/5); % width of each topoplot
if topowidth> 0.25*pos(4) % dont make too large (more than 1/4 of axes width)!
  topowidth = 0.25*pos(4);
end

halfn = floor(ntopos/2);
if rem(ntopos,2) == 1  % odd number of topos
   topoleft = pos(3)/2 - (ntopos/2+halfn*head_sep)*topowidth;
else % even number of topos
   topoleft = pos(3)/2 - ((halfn)+(halfn-1)*head_sep)*topowidth;
end
topoleft = topoleft - 0.01; % adjust left a bit for colorbar

if max(plotframes) > frames |  min(plotframes) < 1
    fprintf('Requested map frame %d is outside data range (1-%d)\n',max(plotframes),frames);
    return
end

%
%%%%%%%%%%%%%%%%%%%% Print times and frames %%%%%%%%%%%%%%%%%%%%%%%%%%
%

fprintf('Scalp maps will show latencies: ');
for t=1:ntopos
  fprintf('%4.0f ',plottimes(t));
end
fprintf('\n');
fprintf('                     at frames: ');
for t=1:ntopos
  fprintf('%4d ',plotframes(t));
end
fprintf('\n');

%
%%%%%%%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%% site the plot at bottom of the figure %%%%%%%%%%%%%%%%%%
%
axdata = axes('Units','Normalized','Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],'FontSize',axfont);
set(axdata,'Color',BACKCOLOR);

limits = get(axdata,'Ylim');
set(axdata,'GridLineStyle',':')
set(axdata,'Xgrid','off')
set(axdata,'Ygrid','on')
axes(axdata)
axcolor = get(gcf,'Color');
set(axdata,'Color',BACKCOLOR);
pl=plot(x,data(plotchans,:));    % plot the data
if length(plotchans)==1
  set(pl,'color','k');
  set(pl,'linewidth',2);
end
xl= xlabel('Latency (ms)');
set(xl,'FontSize',axfont);
yl=ylabel('Potential (\muV)');
set(yl,'FontSize',axfont,'FontAngle','normal');
axis([xmin xmax ymin ymax]);
hold on

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot zero time line %%%%%%%%%%%%%%%%%%%%%%%%%%
%

if xmin<0 & xmax>0
   plot([0 0],[ymin ymax],'k:','linewidth',1.5);
else
  fprintf('xmin %g and xmax %g do not cross time 0.\n',xmin,xmax)
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw vertical lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
width  = xmax-xmin;
height = ymax-ymin;
lwidth = 1.5;  % increment line thickness

for t=1:ntopos % dfraw vertical lines through the data at topoplot frames
 if length(plotchans)>1 | voffsets(t)
  l1 = plot([plottimes(t) plottimes(t)],...
       [min(data(plotchans,plotframes(t))) ...
       voffsets(t) + max(data(plotchans,plotframes(t)))],'w'); % white underline behind
  set(l1,'linewidth',2);
  l1 = plot([plottimes(t) plottimes(t)],...
       [min(data(plotchans,plotframes(t))) ...
       voffsets(t) + max(data(plotchans,plotframes(t)))],'b'); % blue line
  set(l1,'linewidth',lwidth);
 end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw oblique lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
axall = axes('Position',pos,...
             'Visible','Off','FontSize',axfont);   % whole-gca invisible axes
axes(axall)
set(axall,'Color',BACKCOLOR);
axis([0 1 0 1])
  axes(axall)
  axis([0 1 0 1]);
  set(gca,'Visible','off'); % make whole-figure axes invisible

for t=1:ntopos % draw oblique lines through to the topoplots 
  maxdata = max(data(:,plotframes(t))); % max data value at plotframe
  axtp = axes('Units','Normalized','Position',...
       [pos(1)+topoleft+(t-1)*(1+head_sep)*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth ...
                       topowidth*(1+head_sep)]); % this will be the topoplot axes
                       % topowidth]); % this will be the topoplot axes
  axis([-1 1 -1 1]);

  from = changeunits([plottimes(t),maxdata],axdata,axall); % data axes
  to   = changeunits([0,-0.74],axtp,axall);                % topoplot axes
  delete(axtp);
  axes(axall);                                             % whole figure axes
  l1 = plot([from(1) to(1)],[from(2) to(2)]);
  set(l1,'linewidth',lwidth);

  hold on
  set(axall,'Visible','off');
  axis([0 1 0 1]);
end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%
topoaxes = zeros(1,ntopos);
for t=1:ntopos
       % [pos(3)*topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
  axtp = axes('Units','Normalized','Position',...
       [topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
              pos(2)+0.66*pos(4) ...
                  topowidth topowidth*(1+head_sep)]);
  axes(axtp)                             % topoplot axes
  topoaxes(t) = axtp; % save axes handles
  cla

  if ~isempty(varargin)
    topoargs = varargin;
  else
    topoargs = {};
  end
  if topowidth<0.12
      topoargs = { topoargs{:} 'electrodes' 'off' };
  end
  topoplot( data(:,plotframes(t)),chan_locs, topoargs{:});

  % Else make a 3-D headplot
  %
  % headplot(data(:,plotframes(t)),'chan.spline'); 
  
  timetext = [num2str(plottimes(t),'%4.0f')];
  % timetext = [num2str(plottimes(t),'%4.0f') ' ms']; % add ' ms'
  text(0.00,0.80,timetext,'FontSize',axfont-3,'HorizontalAlignment','Center'); % ,'fontweight','bold');
end

%
%%%%%%%%%%%%%%%%%%% Plot a topoplot colorbar %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
axcb = axes('Position',[pos(1)+pos(3)*0.995 pos(2)+0.62*pos(4) pos(3)*0.02 pos(4)*0.09]);
h=cbar(axcb);                        % colorbar axes
pos_cb = get(axcb,'Position');
set(h,'Ytick',[]);

axes(axall)
set(axall,'Color',axcolor);
%
%%%%%%%%%%%%%%%%%%%%% Plot the color bar '+' and '-' %%%%%%%%%%%%%%%%%%%%%%%%%%
%
text(0.986,0.695,'+','FontSize',axfont,'HorizontalAlignment','Center');
text(0.986,0.625,'-','FontSize',axfont,'HorizontalAlignment','Center');

%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the plot title if any %%%%%%%%%%%%%%%%%%%%%%%%%%
%
% plot title between data panel and topoplots (to avoid crowding at top of figure), on the left
ttl = text(0.03,0.635,titl,'FontSize',titlefont,'HorizontalAlignment','left'); % 'FontWeight','Bold');

% textent = get(ttl,'extent');
% titlwidth = textent(3);
% ttlpos = get(ttl,'position');
% set(ttl,'position',[     ttlpos(2), ttlpos(3)]);

axes(axall)
set(axall,'layer','top'); % bring component lines to top
for t = 1:ntopos
  set(topoaxes(t),'layer','top'); % bring topoplots to very top
end

  if ~isempty(varargin)
    try,
		if ~isempty( strmatch( 'absmax', varargin))
			text(0.86,0.624,'0','FontSize',axfont,'HorizontalAlignment','Center');
		end;
	catch, end;
  end

%
% Turn on axcopy()
%
axcopy(gcf);
