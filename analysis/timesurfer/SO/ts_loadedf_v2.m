function data = ts_loadedf(varargin)
% Reading continuous data:
% datafile = '/somepath/cl_screening_v2_032010.edf';
% data     = ts_loadedf('datafile',datafile);
% save(outfile,'data','-v7.3');
% 
% Reading epoched data:
% datafile = '/somepath/cl_screening_v2_032010.edf';
% evntfile = '/somepath/cl_screening_v2_032010.evt';
% prestim  = .2;
% poststim = .6;
% data     = ts_loadedf('datafile',datafile,'evntfile',evntfile,'prestim',prestim,'poststim',poststim);
%
% Example evntfile:
% ----------------------------------------
% 15-Apr-2010
% filename= cl_screening_v2_032010.edf
% SR= 1998.8242
% sample point	stim_num	stim_category
% 10682	0	fixation
% 20712	1	face
% 22744	2	object
% 24710	2	object
% 26710	1	face
% ----------------------------------------
% Note: SR = sampling rate (Hz).
% This would produce epoch_data with three conditions (event codes 0,1,2).
% The first condition would have one trial and the other two would have two
% each.
% 
% Troubleshooting:
% If you have a problem with sopen.m, you need a newer version of biosig
% which is part of EEGLAB.  Solution: download the latest version of EEGLAB.
% 
% If you have a conflict between biosig and matlab versions of str2double.m:
% make a copy of this function, uncomment the commands below (calls to
% addpath, unix, and the lines setting mv1 & mv2). Also, in those commands,
% change the paths to str2double.m.
% 
% Created by JSS on April-2010

addpath(genpath('/home/jsherfey/svn/dev/packages/eeglab7_2_9_18b'));
mv1      = 'mv /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double_rmbyJSS.m /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double.m';
mv2      = 'mv /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double.m /home/jsherfey/svn/dev/packages/eeglab7_2_9_18b/external/biosig-20090130/t200/str2double_rmbyJSS.m';

parms = mmil_args2parms( varargin, ...
                         { 'datafile',[],[],...
                           'evntfile',[],[],...
                           'channels',[],[],...
                           'prestim' ,.2,[],...
                           'poststim',.6,[],...
                         }, ...
                         false );

if isempty(parms.datafile)
	% ask user
  [filename,pathname] = uigetfile({'*.mat;*.fif;*.eeg;*.avg;*.cnt;*.vhdr;*.set'},'Pick a file.','MultiSelect','on');
  if isequal(filename,0) || isequal(pathname,0), return; end
  if iscell(filename)
    parms.datafile = cellfun(@(x)fullfile(pathname,x),filename,'uniformoutput',false);
  else
    parms.datafile = [pathname filename];
  end
end
if iscell(parms.datafile)
  if length(parms.datafile) == 1
    parms.datafile = parms.datafile{1};
  else
    error('%s: support is limited to importing one file only.',mfilename);
  end
end
% read continuous data
cfg = [];
cfg.dataset    = parms.datafile;
cfg.continuous = 'yes';
unix(mv1);
dat  = preprocessing(cfg);
unix(mv2);

if isempty(parms.evntfile)
  % convert to timesurfer format
  data = ts_fieldtrip2data(dat,'epoch');
elseif ~exist(parms.evntfile,'file')
  error('Event file does not exist.');
else
  % read edf header
  unix(mv1); 
  hdr = read_header(parms.datafile);
  unix(mv2);
  nsamp = hdr.nSamples * hdr.nTrials;
  nchan = hdr.nChans;

  if isempty(parms.channels)
    chans = 1:nchan;
  else
    chans = parms.channels;
  end

  % read event file
  [content,res] = readtext(parms.evntfile,'\t');

  Fs   = content{3,1};
  Fs   = Fs(regexp(Fs,'[\d.]'));
  Fs   = str2num(Fs);
  T    = [-parms.prestim:1/Fs:parms.poststim];

  evnt = content(res.numberMask(:,1),:);
  samp = [evnt{:,1}]';
  keep = ~((samp-parms.prestim*Fs)<1 | (samp+parms.poststim*Fs)>nsamp); % not outbounds
  samp = samp(keep);
  cond = [evnt{keep,2}]';
  name = evnt(keep,3);

  [conds,ii] = unique(cond);
  names      = name(ii);
  nevnt      = length(samp);
  ncond      = length(conds);

  begsample = round(samp - parms.prestim*Fs );
  endsample = round(samp + parms.poststim*Fs);

  % initialize timesurfer structure
  tmp          = rmfield(dat,{'trial','time'});
  tmp.trial{1} = dat.trial{1}(:,1);
  tmp.time{1}  = dat.time{1}(1);
  data         = ts_fieldtrip2data(tmp,'epoch');
  data.epochs.time = [];
  data.epochs.data = [];
  data.epochs(1:ncond) = data.epochs;
  clear tmp
  
  % epoch the continuous data
  for c = 1:ncond
    ii  = cond==conds(c);
    s0  = num2cell(begsample(ii));
    sf  = num2cell(endsample(ii));
    tmp = cellfun(@(x,y)(dat.trial{1}(chans,x:y)),s0,sf,'UniformOutput',false);
    data.epochs(c).event_code = conds(c);
    data.epochs(c).time       = T;
    data.epochs(c).num_trials = length(s0);
    data.epochs(c).data       = cat(3,tmp{:}); % channel x time x trial
    clear tmp
  end
end
