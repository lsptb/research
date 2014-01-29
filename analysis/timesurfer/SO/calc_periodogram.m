function psd = calc_periodogram(data,foi,windowspec)
% calculate psd over complete time series

% windows:
windows = {'barthannwin','bartlett','blackman','blackmanharris','bohmanwin','chebwin',...
'flattopwin','gausswin','hamming','hann','kaiser','nuttallwin','parzenwin',...
'rectwin','taylorwin','triang','tukeywin'};

% sampling rate
Fs = data.sfreq;

% freqs
if nargin < 2, foi    = 5:5:round(Fs)/2;  end

% data descriptors
[datatype,datafield,dataparam] = ts_object_info(data);

% sizes
ncond = length(data.(datafield));
nchan = data.num_sensors;
npts  = size(data.(datafield)(1).(dataparam{1}),2);
ntrl  = size(data.(datafield)(1).(dataparam{1}),3);
nfrq  = length(foi);

% window
try
  if nargin < 3
    w = []; 
  elseif ischar(windowspec) && ismember(windowspec,windows)
    w = eval([windowspec '(npts);']);
  %   w = window(str2func(windowspec),npts);
  elseif ischar(windowspec) && ismember(windowspec(2:end),windows) % if @ is in string
    w = eval([windowspec(2:end) '(npts);']);
  %   w = window(str2func(windowspec(2:end)),npts);
  elseif strcmp(class(windowspec),'function_handle') && ismember(func2str(windowspec),windows)
    w = eval([func2str(windowspec) '(npts);']);
  %   w = window(windowspec,npts);
  elseif isnumeric(windowspec) && length(windowspec)==npts
    w = windowspec;
  else
    fprintf('Warning! Window specification not recognized; using default rectangular window.\n');
    w = [];
  end
catch
  fprintf('Warning! Failed to define window function; using default rectangular window.\n');
  w = [];
end

% initialization & preallocation of psd data structure
psd                             = rmfield(data,datafield);
psd.sfreq                       = 1/mean(diff(foi));
psd.(datafield)                 = rmfield(data.(datafield),dataparam);
[psd.(datafield)(1:ncond).time] = deal(foi);
tmp                             = zeros(nchan,nfrq,ntrl);
[psd.(datafield)(1:ncond).data] = deal(tmp);
clear tmp

% calculate periodogram
for c = 1:ncond
  for ch = 1:nchan
    for trl = 1:ntrl
      tmp   = double(data.(datafield)(c).(dataparam{1})(ch,:,trl));
      psd.(datafield)(c).data(ch,:,trl) = periodogram(tmp,w,foi,round(Fs));
      clear tmp
    end
  end
end
