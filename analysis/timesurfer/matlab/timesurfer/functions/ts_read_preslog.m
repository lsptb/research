function evnts = ts_read_preslog(fname,varargin);
% function evnts = ts_read_preslog(fname,[options])
%
% Required Input:
%   fname - input presentation log file
%
% Optional Input:
%   'reject_events' - cell array of event names to exclude
%      (will use regexp, so partial names are ok)
%     {default = []}
%
% created:         08/22/08    by Don Hagler
% last modified:   08/22/08    by Don Hagler
%

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'reject_events',[],[],...
});

evnts = [];

if ~isempty(parms.reject_events) && ~iscell(parms.reject_events)
  parms.reject_events = {parms.reject_events};
end;

reject_types = {'Quit'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(fname);
if fid==-1
  error('failed to open file %s',fname);
end;
[trials] = ...
  textscan(fid,'%d%s%s%d%d%d%d%d%[^\n]','delimiter','\t','headerlines',5);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(trials{1})
  type = trials{2}{i};
  if isempty(type) || ismember(type,reject_types), continue; end;

  code = trials{3}{i};
  skipflag = 0;
  for j=1:length(parms.reject_events)
    if ~isempty(regexp(lower(code),lower(parms.reject_events{j})))
      skipflag=1;
      break;
    end;
  end;
  if skipflag, continue; end;

  if isnan(str2double(code))
    name = code;
    n = regexp(code,'^(?<condition>\d+)\s','names');
    if isempty(n)
      condition = 0;
    else
      condition = str2double(n.condition);
    end;
  else
    name = [];
    condition = str2double(code);
  end;

  if isempty(name) && strcmp(lower(type),'response')
    name = 'response';
  end;

  latency = round(double(trials{4}(i))/10);

  evnts(end+1).type = 'trigger';
  evnts(end  ).duration = 1;
  evnts(end  ).condition = condition;
  evnts(end  ).name = name;
  evnts(end  ).latency = latency;
end;

