function outdata = ts_trials2avg(data,varargin)
% function outdata = ts_trials2avg(data)
% Purpose: Averages epochs or power across trials and returns the appropriate
% data structure.
%
% Input:  timefreq_data or epoch_data with trials
% Output: timefreq_data or avg_data with averages
%
% Created by:       Jason Sherfey 02-Feb-2009
% Last Modified by: Jason Sherfey 05-May-2009
if nargin < 1, help(mfilename); end
data  = ts_checkdata_header(data);
parms = mmil_args2parms(varargin,...
						{'events',[],[],...
						 'conditions',[],[],...
						 'toilim',[],[],...
             'toi',[],[],...
             'keep_event_codes',1,[],...
             'stdev','data',[],...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...                
						},false);

[datatype,datafield,dataparam]   = ts_object_info(data,varargin{:});          

% extract events
if isempty(parms.events) && isempty(parms.conditions)
  conds = 1:length(data.(datafield));
elseif ~isempty(parms.conditions)
  conds = parms.conditions;
elseif ~isempty(parms.events)
  if iscell(parms.events)
    parms.events = unique([parms.events{:}]);
  end
  conds = find(ismember([data.(datafield).event_code],parms.events));
end
data.(datafield) = data.(datafield)(conds);

% abort if data is an average
if strcmp(datatype,'avg_data') || (strcmp(datatype,'timefreq_data') && ndims(data.(datafield)(1).(dataparam{1})) == 3)
  warning('%s: aborting trial averaging: data provided is already averaged over trials',mfilename);
  return;
end

% add stdev field to indata
if ischar(parms.stdev) && any(strcmp(parms.stdev,dataparam))
  [data.(datafield).stdev] = deal([]);
end

%% Create average
outdata = rmfield(data,datafield);
for c = 1:length(data.(datafield))
  outdata.averages(c) = data.(datafield)(c);
  for k = 1:length(dataparam)
    parm = dataparam{k};
    ndim = ndims(data.(datafield)(c).(parm));
    if ndim < 3
      outdata.averages(c).(parm) = data.(datafield)(c).(parm);
      continue;
    end
    outdata.averages(c).(parm) = [];
    outdata.averages(c).(parm) = mean(data.(datafield)(c).(parm),ndim);    
    if ~isempty(parms.stdev) && strcmp(parms.stdev,parm)
      outdata.averages(c).stdev = nanstd(double(data.(datafield)(c).(parm)),0,ndim);
    end
    data.(datafield)(c).(parm) = [];
  end
end