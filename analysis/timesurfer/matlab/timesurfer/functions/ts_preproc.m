function data = ts_preproc(data,varargin)
% input:    epoch_data
% output:   epoch_data or avg_data
% 
%  'bandpass_flag' - [1|0] Toggle bandpass fft filter before averaging
%     {default: 0}
%  'bandpass_low_cf'  - low cutoff frequency (high-pass filter) (Hz)
%     { default: 0 }
%  'bandpass_low_tb'  - low cutoff transition band (Hz)
%     { default: 0 }
%  'bandpass_high_cf' - high cutoff frequency (low-pass filter) (Hz)
%     { default: 100 }
%  'bandpass_high_tb' - high cutoff transition band (Hz)
%     { default: 0 }
%  'bandpass_detrend_flag' - whether to subtract a linear fit before
%  filtering {default: 1}
%  'dsfact'  - downsampling factor -- must be an integer
%     data is downsampled to lower sampling frequency
%     e.g. original sampling freq = 1000, dsfact = 4,
%         resulting sampling freq = 250
%     { default: 1 (no downsampling) }
%  'detrend_flag' - [1|0] Toggle detrending of single trials
%     before averaging
%     {default: 1}
%  'baseline_flag' - [1|0] Toggle baseline subtraction of single trials
%     before averaging
%     {default: 1}
%  'baseline_start' - start time of baseline period (msec)
%     relative to trigger onset; negative times occur before trigger
%     { default: -Inf } (start at beginning of prestimulus period)
%  'baseline_end'   - end time of baseline period (msec)
%     relative to trigger onset; negative times occur before trigger
%     { default: 0 } (end at trigger onset)
%  'ncov_ex_evnts' - vector of event codes that should not
%     be used in calculating the noise covariance matrix
%     use this, for example, to exclude events with pre-stimulus activity
%     that is not noise -- for example, a button press
%     { default: [] }
%  'badchans'    - vector of bad channel indices -- will be set to zero
%     { default: [] }
%  'badchanfile' - name of text file containing bad channel labels
%    {default: []}
%   'combinations'  - this is the list of combinations to produce 
%                provided as a cell array.  Each combination is
%                supplied as a string.
%   'calc' - 'weighted', 'avg', or 'sum' - specify whether to perform
%     a weighted average, straight average, or addition when combining
%     conditions with '+' operation
%     {default = 'weighted'}
%   'neweventcodes' - a list of new event codes to assign each of the new
%     combinations, if none is supplied it will default creating new event
%     codes starting with the largest current event code in the data
%     structure conditions.  It is recommended to specify your own.
%   'fieldtrip_flag' - [0|1] Toggle fieldtrip preprocessing { default: 0 }
%     this flag determines what functions will be used to perform the
%     preprocessing (fieldtrip or timesurfer).
%   'saveepochs_flag' - [0|1]
%   'saveaverages_flag' - [0|1]
%   'returnepochs_flag' - [0|1]
%   'returnaverages_flag' - [0|1]
%
% FIELDTRIP options (if fieldtrip_flag is set to 1; you can set the following parameters):
% The configuration can contain
%   lpfilter      = 'no' or 'yes'  lowpass filter
%   hpfilter      = 'no' or 'yes'  highpass filter
%   bpfilter      = 'no' or 'yes'  bandpass filter
%   bsfilter      = 'no' or 'yes'  bandstop filter
%   lnfilter      = 'no' or 'yes'  line noise removal using notch filter
%   dftfilter     = 'no' or 'yes'  line noise removal using discrete fourier transform
%   medianfilter  = 'no' or 'yes'  jump preserving median filter
%   lpfreq        = lowpass  frequency in Hz
%   hpfreq        = highpass frequency in Hz
%   bpfreq        = bandpass frequency range, specified as [low high] in Hz
%   bsfreq        = bandstop frequency range, specified as [low high] in Hz
%   lnfreq        = line noise frequency in Hz, default 50Hz
%   dftfreq       = line noise frequencies for DFT filter, default [50 100 150] Hz
%   lpfiltord     = lowpass  filter order
%   hpfiltord     = highpass filter order
%   bpfiltord     = bandpass filter order
%   bsfiltord     = bandstop filter order
%   lnfiltord     = line noise notch filter order
%   lpfilttype    = digital filter type, 'but' (default) or 'fir'
%   hpfilttype    = digital filter type, 'but' (default) or 'fir'
%   bpfilttype    = digital filter type, 'but' (default) or 'fir'
%   bsfilttype    = digital filter type, 'but' (default) or 'fir'
%   lpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   hpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   bpfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   bsfiltdir     = filter direction, 'twopass' (default), 'onepass' or 'onepass-reverse'
%   medianfiltord = length of median filter
%   blc           = 'no' or 'yes'
%   blcwindow     = [begin end] in seconds, the default is the complete trial
%   detrend       = 'no' or 'yes', this is done on the complete trial
%   polyremoval   = 'no' or 'yes', this is done on the complete trial
%   polyorder     = polynome order (default = 2)
%   hilbert       = 'no', 'abs', 'complex', 'real', 'imag', 'absreal', 'absimag' or 'angle' (default = 'no')
%   rectify       = 'no' or 'yes'
%   precision     = 'single' or 'double' (default = 'double')
%
% Preprocessing options that you should only use for EEG data are
%   reref         = 'no' or 'yes' (default = 'no')
%   refchannel    = cell-array with new EEG reference channel(s)
%   implicitref   = 'label' or empty, add the implicit EEG reference as zeros (default = [])
%   montage       = 'no' or a montage structure (default = 'no')
%
%
% Created by Jason Sherfey on 10-Apr-2009
% Modified last by JSS on 12-Apr-2009

parms = mmil_args2parms(varargin,...
						{'timesurfer_flag',1,{0,1},...
             'events',[],[],...
             'fieldtrip_flag',0,{0,1},...
             'bandpass_flag',false,[false true],...
             'bandpass_detrend_flag',true,[false true],...
             'bandpass_baseline_flag',false,[false true],...
             'bandpass_low_cf',0.2,[],...
             'bandpass_low_tb',0,[],...
             'bandpass_high_cf',30,[],...
             'bandpass_high_tb',0,[],...
             'notch_flag',false,[false true],...
             'notch_cf',[],[],...
             'lowpass_flag',false,[false true],...
             'lowpass_cf',[],[],...
             'lowpass_tb',0,[],...
             'highpass_cf',[],[],...
             'highpass_tb',0,[],...
             'dsfact',1,[],...
             'detrend_flag',false,[false true],...
             'baseline_flag',false,[false true],...
             'baseline_start',-Inf,[-Inf,Inf],...
             'baseline_end',Inf,[-Inf,Inf],...     
             'combinations',[],[],...
             'comboeventcodes',[],[],...
             'calc','weighted',{'weighted','avg','sum'},...
             'cfg',[],[],...
             'feedback','none',{'non','gui','dial','textbar','text','textcr','textnl','none','no','yes'},...
             'lpfilter','no',{'yes','no'},...
             'hpfilter','no',{'yes','no'},...
             'bpfilter','no',{'yes','no'},...
             'bsfilter','no',{'yes','no'},...
             'lnfilter','no',{'yes','no'},...
             'dftfilter','no',{'yes','no'},...
             'medianfilter','no',{'yes','no'},...
             'lpfreq',30,[],...
             'hpfreq',[],[],...
             'bpfreq',[],[],...
             'lnfreq',60,[],...
             'dftfreq',[60 120 180],[],...              
             'lpfiltord',6,[],...
             'hpfiltord',6,[],...
             'bpfiltord',4,[],...
             'bsfiltord',4,[],...
             'lnfiltord',4,[],...
             'lpfilttype','but',{'but','fir'},...
             'hpfilttype','but',{'but','fir'},...
             'bpfilttype','but',{'but','fir'},...
             'bsfilttype','but',{'but','fir'},...
             'lpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'hpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'bpfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'bsfiltdir','twopass',{'onepass','onepass-reverse','twopass'},...
             'medianfiltord',9,[],...
             'blc','no',{'yes', 'no'},...
             'blcwindow',[],[],...
             'detrend','no',{'yes','no'},...     
             'polyremoval','no',{'yes','no'},...
             'polyorder',2,[],...
             'hilbert','no',{'no','abs','complex','real','imag','absreal','absimag','angle'},...
             'rectify','no',{'yes','no'},...
             'precision',[],{'single','double'},...
             'reref','no',{'yes','no'},...
             'refchannel',[],[],...
             'implicitref',[],[],...
             'montage','no',[],...
             'saveepochs_flag',0,{0,1},...
             'saveaverages_flag',0,{0,1},...
             'returnepochs_flag',1,{0,1},...
             'returnaverages_flag',0,{0,1},...
             'saveavgs','no',{'yes','no'},...
             'saveepochs','no',{'yes','no'},...
             'returnavgs','no',{'yes','no'},...
             'returnepochs','no',{'yes','no'},...
             'rootoutdir',pwd,[],...
             'prefix','preproc',[],...
             'filename',[],[],...
             'overwrite',0,{0,1},...
             'verbose',1,{0,1},...
             'logfile',      [],[],...
             'logfid',       [1],[], ...                
						 },false);

data  = ts_checkdata_header(data,'events',parms.events);

% store original data structure parms
if isfield(data,'parms'), parms.previous = data.parms; end

[datatype,datafield,dataparam] = ts_object_info(data,varargin{:});
% data = ts_data_selection(data,varargin{:});

if isempty(parms.precision),       parms.precision = class(data.(datafield)(1).data); end
if strcmp(parms.saveavgs,'yes'),   parms.saveaverages_flag = 1; end
if strcmp(parms.saveepochs,'yes'), parms.saveepochs_flag = 1;   end
if strcmp(parms.returnavgs,'yes'),   parms.returnaverages_flag = 1; end
if strcmp(parms.returnepochs,'yes'), parms.returnepochs_flag = 1;   end

avgflag = parms.saveaverages_flag || parms.returnaverages_flag;
epoflag = parms.saveepochs_flag || parms.returnepochs_flag;

mmil_logstr(parms,'%s: preprocessing %g conditions, %g channels\n',mfilename,length(data.(datafield)),data.num_sensors);
% fprintf('%s: preprocessing %g conditions, %g channels\n',mfilename,length(data.(datafield)),data.num_sensors);

%% TimeSurfer preprocessing    
if parms.timesurfer_flag && ~parms.fieldtrip_flag
  % convert fieldtrip parameters to timesurfer parameters
  if parms.baseline_flag || strcmp(parms.blc,'yes')
    if strcmp(parms.blc,'yes'), 
      parms.baseline_flag = 1;
    end
    if isnumeric(parms.blcwindow) && ~isempty(parms.blcwindow)
      parms.baseline_start = parms.blcwindow(1);
      parms.baseline_end   = parms.blcwindow(end);
    elseif ischar(parms.blcwindow) && strcmp(parms.blcwindow,'all')
      parms.baseline_start = -Inf;
      parms.baseline_end   =  Inf;
    end
  end
  if parms.lowpass_flag || strcmp(parms.lpfilter,'yes')
    if strcmp(parms.lpfilter,'yes')
      parms.lowpass_flag = 1;
    end
    if isempty(parms.lowpass_cf)
      parms.lowpass_cf = parms.lpfreq;
    end
  end
  if parms.bandpass_flag || strcmp(parms.bpfilter,'yes')
    if strcmp(parms.bpfilter,'yes')
      parms.bandpass_flag = 1;
    end
    if ~isempty(parms.bpfreq)
      parms.bandpass_low_cf  = parms.bpfreq(1);
      parms.bandpass_high_cf = parms.bpfreq(2);
    end    
  end
  if parms.notch_flag || strcmp(parms.lnfilter,'yes')
    if strcmp(parms.lnfilter,'yes')
      parms.notch_flag = 1;
    end
    if isempty(parms.notch_cf)
      parms.notch_cf = parms.lnfreq;
    end
  end  
  if strcmp(parms.detrend,'yes')
    parms.detrend_flag = 1;
  end
  % loop over conditions
  for c = 1:length(data.(datafield))
    if ~isempty(parms.events) && ~ismember(data.(datafield)(c).event_code,parms.events)
      continue;
    end
    mmil_logstr(parms,'%s: preprocessing %g of %g\n',mfilename,c,length(data.(datafield)));
%     fprintf('%s: preprocessing %g of %g\n',mfilename,c,length(data.(datafield)));
%     if ischar(parms.blcwindow) && strcmp(parms.blcwindow,'all')
%       parms.baseline_start = data.(datafield)(c).time(1);
%       parms.baseline_end   = data.(datafield)(c).time(end);
%     end
    % get dimensions for preallocation
    sx = size(data.(datafield)(c).data,1);
    sy = size(data.(datafield)(c).data,2);
    sz = size(data.(datafield)(c).data,3);
    if parms.dsfact ~= 1
      data.(datafield)(c).time = downsample(data.(datafield)(c).time,parms.dsfact);
      sy = length(data.(datafield)(c).time);
    end
    % preallocate space for the preprocessed data
    tmp = zeros(sx,sy,sz);
    % convert to single precision
    if ~strcmp(parms.precision,'double')
      tmp = cast(tmp,parms.precision);
    end    
    % loop over trials
    for trl = 1:sz
      % extract one trial
      dat = squeeze(data.(datafield)(c).data(:,:,trl));
      if trl == sz
        data.(datafield)(c).data = [];
      end
      % force to double precision
      if ~isa(dat,'double')
        dat = double(dat);
      end
      % lowpass filter
      if parms.lowpass_flag
        dat = ts_freq_filt(dat',data.sfreq,parms.lowpass_cf,parms.lowpass_tb,'lowpass')';
      end      
      % bandpass filter
      if parms.bandpass_flag
        if parms.bandpass_baseline_flag
          dat = blc(dat);
        end        
        szthresh = 100; % arbitrary size threshold
        if numel(dat)*8*2/10E6 > szthresh
          % loop over channels to handle large matrices
          if parms.bandpass_detrend_flag
            for ch = 1:size(dat,1)
              dat(ch,:) = detrend(dat(ch,:)')';
            end
          end
          for ch = 1:size(dat,1)
            dat(ch,:) = ts_freq_filt( dat(ch,:)',data.sfreq,[parms.bandpass_low_cf,parms.bandpass_high_cf],...
                                [parms.bandpass_low_tb,parms.bandpass_high_tb],'bandpass' )';
          end
        else
          if parms.bandpass_detrend_flag
            dat = detrend(dat')';
          end
          dat = ts_freq_filt( dat',data.sfreq,[parms.bandpass_low_cf,parms.bandpass_high_cf],...
                              [parms.bandpass_low_tb,parms.bandpass_high_tb],'bandpass' )';        
        end
      end
      % highpass filter
      if strcmpi(parms.hpfilter,'yes')
        dat = ts_freq_filt(dat',data.sfreq,parms.hpfreq,parms.highpass_tb,'highpass')';
      end
      % notch filter
      if parms.notch_flag
        for fc = 1:length(parms.notch_cf)
          dat = ts_freq_filt(dat',data.sfreq,parms.notch_cf(fc),5,'notch')';
        end
      end
      % convert to single precision
      if ~strcmp(parms.precision,'double')
        dat = cast(dat,parms.precision);
      end      
      % downsample
      if parms.dsfact ~= 1
        dat = downsample(dat',parms.dsfact)'; 
      end
      % linear trend removal
      if parms.detrend_flag
        dat = detrend(dat')';
      end
      % baseline correction
      if parms.baseline_flag
        bid = [nearest(data.(datafield)(c).time,parms.baseline_start),...
               nearest(data.(datafield)(c).time,parms.baseline_end)];
        bmu = mean(dat(:,bid(1):bid(2)),2);
        dat = dat - bmu*ones(1,size(dat,2));
        clear bmu bid
      end
      tmp(:,:,trl) = dat;
      clear dat
    end
%     % convert to single precision
%     if ~strcmp(parms.precision,'double')
%       tmp = cast(tmp,parms.precision);
%     end
    % remove data from structure prior to copy to prevent out of memory error
    data.(datafield)(c).data = [];
    data.(datafield)(c).data = tmp;
    clear tmp
  end
  data.sfreq = data.sfreq / parms.dsfact;
  if avgflag
    data = ts_trials2avg(data);
  end
end

%% FieldTrip preprocessing
if parms.fieldtrip_flag
  n = 1;
  for c = 1:length(data.(datafield))
    if ~isempty(parms.events) && ~ismember(data.(datafield)(c).event_code,parms.events)
      continue;
    end
    mmil_logstr(parms,'%s: preprocessing %g of %g\n',mfilename,c,length(data.(datafield)));
%     fprintf('%s: preprocessing %g of %g\n',mfilename,c,length(data.(datafield)));
    if ~isempty(parms.cfg)
      cfg = parms.cfg;
    else
      cfg = parms;
    end
    % force to double precision
    if ~isa(data.(datafield)(c).data,'double')
      data.(datafield)(c).data = double(data.(datafield)(c).data);
    end
    tmp = ts_data2fieldtrip(data,'condition',c,'dimord','chan_time');
    if avgflag
      cfg.keeptrials = 'no';
    else
      cfg.keeptrials = 'yes';
    end
    ftout{n} = timelockanalysis(cfg,tmp);
    clear tmp;
    if ~strcmp(parms.precision,'double')
      tmp = cast(tmp,parms.precision);
    end
    n = n + 1;
  end
  data = rmsubfield(data,'epochs.data');
  if avgflag
    data = ts_fieldtrip2data(ftout,'averages',data);
  else
    data = ts_fieldtrip2data(ftout,'epochs',data);
  end
  clear ftout
end

% calculate composite conditions
if ~isempty(parms.combinations)
  data = ts_combine_conditions(data,'combinations',parms.combinations,'neweventcodes',parms.comboeventcodes,'calc',parms.calc);
end

if isfield(parms,'previous'), data.parms = parms; end

if parms.saveepochs_flag && ~avgflag
  for c = 1:length(data.epochs)
    epoch_data = rmfield(data,'epochs');
    epoch_data.epochs(1) = data.epochs(c);
    if isempty(parms.filename) || ~ischar(parms.filename)
      parms.filename = sprintf('%s/%s_epochs_cond%g.mat',parms.rootoutdir,parms.prefix,c);
    end
    if exist(parms.filename) && ~parms.overwrite
      mmil_logstr(parms,'%s: not overwriting %s\n',mfilename,parms.filename);
%       fprintf('%s: not overwriting %s\n',mfilename,parms.filename);
    else
      mmil_logstr(parms,'%s: saving epoch_data to %s\n',mfilename,parms.filename);
%       fprintf('%s: saving epoch_data to %s\n',mfilename,parms.filename);
      save(parms.filename,'epoch_data');
    end
    clear epoch_data;
  end
end
if parms.saveaverages_flag
%   data = ts_trials2avg(data);
  for c = 1:length(data.averages)
    avg_data = rmfield(data,'averages');
    avg_data.averages(1) = data.averages(c);
    if isempty(parms.filename) || ~ischar(parms.filename)
      parms.filename = sprintf('%s/%s_averages_cond%g.mat',parms.rootoutdir,parms.prefix,c);
    end
    if exist(parms.filename) && ~parms.overwrite
      mmil_logstr(parms,'%s: not overwriting %s\n',mfilename,parms.filename);
%       fprintf('%s: not overwriting %s\n',mfilename,parms.filename);
    else
      mmil_logstr(parms,'%s: saving avg_data to %s\n',mfilename,parms.filename);
%       fprintf('%s: saving avg_data to %s\n',mfilename,parms.filename);
      save(parms.filename,'avg_data');
    end
    clear avg_data;
  end  
end
