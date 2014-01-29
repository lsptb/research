function ts_plot_avg(avg_data,varargin);
% ts_plot_avg - plots MEG/EEG sensor waveforms for single channel
%
% Usage:
%   ts_plot_avg(avg_data,'key1', value1,...);
%
% Required input:
%  avg_data - average data structure (see ts_avg_fif_data)
%
% Optional parameters:
%  'conditions' - vector of condition numbers (not event code)
%    used to index avg_data.averages
%    {default: [1]}
%  'legend_flag' - [0|1] toggle display legend
%    {default: 1}
%  'legend_loc' - legend location
%    {default: 'EastOutside'}
%  'axis_flag' - [0|1] toggle display default matlab axis
%    {default: 1}
%  'vert_lines' - vector of time points to make vertical lines
%    {default: [0]}
%  'horz lines' - vector of y-values make horizontal lines
%    {default: [0]}
%  'condnames' - cell array of condition names for legend
%    {default: []}
%  'channame' - channel name
%    supply either channame or channum (channame takes precedence)
%    {default: []}
%  'channum' - channel number
%    {default: 1}
%  'scale_range' - min and max values for scale (in uVolts, fT, or fT/cm)
%    {default: [] -> will use max/min}
%  'time_range' - min and max time points (in msec)
%    {default: [] -> will use max/min}
%  'label_flag' - [0|1] toggle display x and y axis labels and title
%    {default: 1}
%  'fontsize'
%    {default: 12}
%  'fontname'
%    {default: 'Helvetica'}
%  'linewidth' - trace line width
%    {default: 1.5}
%  'grid_linewidth' - line width for horizontal and vertical lines
%    {default: 1}
%  'color_order' - cell array containing matlab color characters
%    {default: {'b' 'g' 'r' 'c' 'm' 'y' 'k'}}
%
%  created:       10/11/07   by Don Hagler
%  last modified: 08/23/08   by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
  help(mfilename);
  return;
end;

nconds = length(avg_data.averages);
nchans = avg_data.num_sensors;
time = avg_data.averages(1).time*1000;
ntpoints = length(time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parms = mmil_args2parms( varargin, {...
  'conditions',     [1:nconds],  [],...
  'condnames',      [],          [],...
  'channum',        1,           [1 nchans],...
  'channame',       [],          [],...
  'scale_range',    [],          [],...
  'time_range',     [],          [],...
  'linewidth',      1.5,         [0.1 100],...
  'grid_linewidth', 1,           [0.1 100],...
  'fontsize',       12,          [1 100],...
  'fontname',       'Helvetica', [],...
  'label_flag',     1,           [0 1],...
  'legend_flag',    1,           [0 1],...
  'legend_loc',    'EastOutside', [],...
  'axis_flag',      1,           [0 1],...
  'vert_lines',     [0],         [],...
  'horz_lines',     [0],         [],...
  'color_order'  {'b' 'g' 'r' 'c' 'm' 'y' 'k'}, [],...
  'bgcolor'         'w',         []...
},1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if isempty(parms.conditions)
  parms.conditions = [1:nconds];
end;
if ~isempty(parms.condnames)
  if ~iscell(parms.condnames), parms.condnames = {parms.condnames}; end;
  if length(parms.condnames) ~= length(parms.conditions)
    error('number of condnames (%d) does not match number of conditions (%d)',...
      length(parms.condnames),length(parms.conditions));
  end;
end;

if isempty(parms.color_order)
  error('color_order is empty');
end;
if ~iscell(parms.color_order), parms.color_order = {parms.color_order}; end;
for i=1:length(parms.color_order)
  if ~ischar(parms.color_order{i})
    error('color_order must contain color strings (e.g. ''b'', ''b*'')')
  end;
end;

if ~isempty(parms.channame)
  parms.channum = find(strcmp(parms.channame,{avg_data.sensor_info.label}));
  if isempty(parms.channum)
    error('channame %s not found in avg_data.sensor_info.label',...
      parms.channame);
  end;
  chan_info = avg_data.sensor_info(parms.channum);
else
  chan_info = avg_data.sensor_info(parms.channum);
  parms.channame = chan_info.label;
end;
chantype = chan_info.typestring;

switch chantype
case {'grad1' 'grad2','grad'}
  scalefact = 10^13;
  units = 'fT/m';
case 'mag'
  scalefact = 10^15;
  units = 'fT';
case 'eeg'
  scalefact = 10^6;
  units = 'uV';
otherwise
  scalefact = 1;
  units = 'unknown';
end; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot data

cla;
whitebg(parms.bgcolor);
hold on;

data = zeros(length(parms.conditions),ntpoints);
j=1;
legstrs = [];
for i=1:length(parms.conditions)
  if (j>length(parms.color_order)) j=1; end;
  color = parms.color_order{j}; j=j+1;
  c = parms.conditions(i);
  if c<1 | c>nconds
    error('bad condition number %d -- avg_data has %d conditions',...
      c,nconds);
  end;
  tmp_data = squeeze(avg_data.averages(c).data(parms.channum,:))*scalefact;
  plot(time,tmp_data,color,'LineWidth',parms.linewidth);
  if parms.legend_flag
    if isempty(parms.condnames)
      legstrs{i} = sprintf('%d: Event Code = %d',...
        c,avg_data.averages(c).event_code);
    else
      legstrs{i} = parms.condnames{i};
    end;
  end;
end;

axis tight;
if ~isempty(parms.scale_range)
  set(gca,'YLim',parms.scale_range);
end;
if ~isempty(parms.time_range)
  set(gca,'XLim',parms.time_range);
end;
set(gca,'FontName',parms.fontname);
if parms.label_flag
  xlabel('Time (ms)','FontSize',parms.fontsize,'FontName',parms.fontname);
  ylabel(units,'FontSize',parms.fontsize,'FontName',parms.fontname);
  title(parms.channame,'FontSize',parms.fontsize,'FontName',parms.fontname,...
       'Color','k');
end;
if parms.legend_flag
  legend(legstrs,'Location',parms.legend_loc);
end;
if ~parms.axis_flag, axis off; end;

for i=1:length(parms.vert_lines)
  x = parms.vert_lines(i);
  plot([x,x],get(gca,'YLim'),'k','LineWidth',parms.grid_linewidth);
end;
for i=1:length(parms.horz_lines)
  y = parms.horz_lines(i);
  plot(get(gca,'XLim'),[y,y],'k','LineWidth',parms.grid_linewidth);
end;

return;

