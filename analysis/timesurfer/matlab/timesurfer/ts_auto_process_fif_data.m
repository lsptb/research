function ts_auto_process_fif_data(indata, varargin)
%function ts_auto_process_fif_data(indata, varargin)
%
% Purpose: creates and runs a script for processing fif data
%   with default settings, either combining across multiple raw fif
%   files in a single input directory, or across fif files given as 
%   an input character cell array.
%
% Usage:
%   ts_auto_process_fif_data('/home/halgdev/incoming/MEG_UCSD/STU_SUBJ_070523')
%     -- This will process fif files in the fictional directional above
%     -- and place output files in the current working directory
%
%   ts_auto_process_fif_data('/home/halgdev/incoming/MEG_UCSD/STU_SUBJ_070523'...
%     /home/halgdev/analysis/MEG_UCSD/STU_SUBJ_070523)
%     -- This will process fif files in the fictional directional above
%     -- and place output files in the analysis directory (will create if necessary)
%
%   ts_auto_process_fif_data( {'/home/ktravis/Subjects/MattH/rawdata/matt_h_test_ses3_raw.fif'} )
%     -- This will process the fif files specified (only one here)
%     -- and place output files in the current working directory
%     -- NOTICE the {} around the FILENAMES (not used for above 
%     -- examples with directories)
%
% Required Input:
%   indata: either:
%          1. full path of input directory containing one or more raw fif
%             files
%             OR
%          2. a cell array of input fif files.
%
% Optional Input:
%  'outdir' - full path of output directory which will contain
%     automatically generated processing script and matfiles
%     subdirectory
%     {default = current working directory}
%  'saveepochs_flag' - [1|0] instead of calculating averages, extract
%     epochs (single trial time courses)
%     {default: 0}
%  'event_recode_rules' - cell array of strings defining how to 
%    change one event's condition to a new value.  See
%    ts_recode_events for more details.
%     {default: []}
%  'trig_minduration' - minimum # of samples for a trigger code to be ON
%    before it is considered ACTUALLY ON and not just noise.
%    {default: []} (if empty, treat all event codes as valid)
%
% Created:  05/23/07 by Don Hagler
% Last Mod: 08/27/07 by Ben Cipollini
%
% See also: ts_process_fif_data
%

if (~mmil_check_nargs(nargin,1)) return; end;

% Accept OLD style to call: only 2 args, second is "outdir"
if (nargin==2)
  varargin = {'outdir', varargin{:}};
end;

% input args are ONLY for specifying WHAT to process
% or HOW to process it, 
% but NOT for numerical values input into the processing stream.
parms = mmil_args2parms(varargin, { ...
  'event_recode_rules', [], [], ...
  'trig_minduration', 5, [1 inf], ...
  'saveepochs_flag', 0, [0 1], ...
  'outdir', pwd, [] ...
},true);

if isempty(parms.outdir) | ismember(parms.outdir,{'.','./'}), parms.outdir = pwd; end;
if isstr(parms.event_recode_rules), parms.event_recode_rules = { parms.event_recode_rules; }; end;

% we got an input list of fif files
if (iscell(indata)) 
  flist   = indata;
  if isempty(flist)
    error('%s: no fif files specified in input (e.g. input was empty!)',mfilename);
  else
    for i=1:length(flist)
      if ~exist(flist{i},'file')
        error('%s: fif file doesn''t exist: %s',mfilename, flist{i});
      end;
    end;
  end;

% we got nothin'
elseif (~exist(indata, 'file'))
  error('%s: input source not found: %s', mfilename, indata);

% we got a single input fif file
elseif (exist(indata, 'file') & ~exist(indata, 'dir'))
  flist = {indata};

% we got a directory to process
else
  indir   = indata;
  if ismember(indir,{'.','./'}), indir = pwd; end;

  fif_files = dir(sprintf('%s/*.fif',indir));
  flist = {};
  for i=1:length(fif_files)
    flist{end+1} = fullfile(indir, fif_files(i).name);
  end;

  if isempty(flist)
    error('%s: no fif files found in input directory %s',mfilename,indir);
  end;
end;

if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('%s: failed to create output directory %s',mfilename,parms.outdir);
  end;
end;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

script_name = 'run_process_data';
% this filename must correspond to the hard-coded value run at bottom of this script.
fname_script = sprintf('%s/%s.m',parms.outdir, script_name);
if exist(fname_script,'file')
  fprintf('%s: running existing script %s...\n',mfilename,fname_script);
else
  fprintf('%s: automatically generating processing script %s...\n',...
    mfilename,fname_script);

  fid = fopen(fname_script,'wt');
  if fid<0
    error('%s: failed to create output script %s\n',mfilename,fname_script);
    return;
  end;

  fprintf(fid,'%% auto_run_process_data\n\n');
  fprintf(fid,'%% automatically generated on %s\n\n',date);
  fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'%% set variables\n');

  fprintf(fid,'datafile={};\n');
  fprintf(fid,'datafile{end+1} = ''%s'';\n',flist{:});

  fprintf(fid,'\n');
  fprintf(fid,'badchanfile = [];\n');
  fprintf(fid,'%% specify a text file listing the bad channels\n');
  fprintf(fid,'valid_event_codes = [];\n');
  fprintf(fid,'%% specify event codes to be averaged\n');
  fprintf(fid,'%%   (empty means include all events found in files)\n');
  if (isempty(parms.event_recode_rules))
    fprintf(fid,'event_recode_rules = [];\n');
  else
%% todo: test that this actually works
    fprintf(fid,'event_recode_rules = { %s };\n', ...
            sprintf('''%s'' ', parms.event_recode_rules{:}));
  end;
  fprintf(fid,'%% specify rules to change event condition values\n');
  fprintf(fid,'%%   (empty means include all events found in files)\n');
  fprintf(fid,'code_excl = [];\n');
  fprintf(fid,'%% event code(s) around which other events are excluded\n');
  fprintf(fid,'%%   (empty means do not exclude events like this)\n');
  fprintf(fid,'time_excl_pre = 500;\n');
  fprintf(fid,'%% duration before code_excl that other events are excluded\n');
  fprintf(fid,'time_excl_post = 500;\n');
  fprintf(fid,'%% duration after code_excl that other events are excluded\n');
  fprintf(fid,'trig_minduration = %d;\n', parms.trig_minduration);
  fprintf(fid,'%% minimum duration for data to appear on a trigger channel\n');
  fprintf(fid,'%% before the trigger is considered to be "on"\n');
  fprintf(fid,'stim_delay = 32.5;\n');
  fprintf(fid,'%% delay between trigger and stimulus onset\n');
  fprintf(fid,'reject_mag = 10000;\n');
  fprintf(fid,'%% auto-rejection threshold for magnetometer channels (fT)\n');
  fprintf(fid,'reject_grad = 6000;\n');
  fprintf(fid,'%% auto-rejection threshold for gradiometer channels (fT/cm)\n');
  fprintf(fid,'reject_eeg = 0;\n');
  fprintf(fid,'%% auto-rejection threshold for EEG channels (uV)\n');
  fprintf(fid,'reject_eog = 200;\n');
  fprintf(fid,'%% auto-rejection threshold for EOG channels (uV)\n');
  fprintf(fid,'prestim_dur = 100;\n');
  fprintf(fid,'%% pre-stimulus duration (msec)\n');
  fprintf(fid,'poststim_dur = 400;\n');
  fprintf(fid,'%% post-stimulus duration (msec)\n');
  fprintf(fid,'bandpass_low_cf = 0.2;\n');
  fprintf(fid,'%% bandpass filter low cut-off frequency (Hz)\n');
  fprintf(fid,'bandpass_low_tb = 0.4;\n');
  fprintf(fid,'%% bandpass filter low transition band (Hz)\n');
  fprintf(fid,'bandpass_high_cf = 50;\n');
  fprintf(fid,'%% bandpass filter high cut-off frequency (Hz)\n');
  fprintf(fid,'bandpass_high_tb = 10;\n');
  fprintf(fid,'%% bandpass filter high transition band (Hz)\n');
  fprintf(fid,'dsfact = 4;\n');
  fprintf(fid,'%% downsampling factor\n');
  fprintf(fid,'baseline_start = -80;\n');
  fprintf(fid,'%% start of baseline period (msec) relative to stim onset\n');
  fprintf(fid,'baseline_end = -5;\n');
  fprintf(fid,'%% end of baseline period (msec) relative to stim onset\n');
  fprintf(fid,'null_event = [];\n');
  fprintf(fid,'%% event code used for "null" subtraction\n');
  fprintf(fid,'post_dsfact = 1;\n');
  fprintf(fid,'%% downsampling factor applied after averaging\n');
  fprintf(fid,'\n');
  fprintf(fid,'bandpass_flag = 1;\n');
  fprintf(fid,'%% whether (0 or 1) to apply bandpass filtering before averaging\n');
  fprintf(fid,'detrend_flag = 1;\n');
  fprintf(fid,'%% whether (0 or 1) to detrend before averaging\n');
  fprintf(fid,'baseline_flag = 1;\n');
  fprintf(fid,'%% whether (0 or 1) to baseline subtract before averaging\n');
  fprintf(fid,'post_subnull_flag = 0;\n');
  fprintf(fid,'%% whether (0 or 1) to do null subtraction after averaging\n');
  fprintf(fid,'post_bandpass_flag = 0;\n');
  fprintf(fid,'%% whether (0 or 1) to apply bandpass filtering after averaging\n');
  fprintf(fid,'post_detrend_flag = 0;\n');
  fprintf(fid,'%% whether (0 or 1) to detrend after averaging\n');
  fprintf(fid,'post_baseline_flag = 0;\n');
  fprintf(fid,'%% whether (0 or 1) to baseline subtract after averaging\n');
  fprintf(fid,'\n');
  fprintf(fid,'saveepochs_flag = %d;\n', parms.saveepochs_flag);
  fprintf(fid,'%% whether (0 or 1) to save trial-by-trial data\n');
  fprintf(fid,'\n');
  fprintf(fid,'browseraw = 0;\n');
  fprintf(fid,'%% launch browser for one of the raw data files\n');
  fprintf(fid,'%%   If 0, no browsing, if 1, browse datafile{1}, if 2, browse datafile{2},...\n');
  fprintf(fid,'\n');
  fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
  fprintf(fid,'%% run ts_process_fif_data\n');
  fprintf(fid,'ts_process_fif_data(...\n');
  fprintf(fid,'  datafile,...\n');
  fprintf(fid,'  ''badchanfile'',badchanfile,...\n');
  fprintf(fid,'  ''valid_event_codes'',valid_event_codes,...\n');
  fprintf(fid,'  ''event_recode_rules'',event_recode_rules,...\n');
  fprintf(fid,'  ''trig_minduration'',trig_minduration,...\n');
  fprintf(fid,'  ''stim_delay'',stim_delay,...\n');
  fprintf(fid,'  ''reject_mag'',reject_mag,...\n');
  fprintf(fid,'  ''reject_grad'',reject_grad,...\n');
  fprintf(fid,'  ''reject_eeg'',reject_eeg,...\n');
  fprintf(fid,'  ''reject_eog'',reject_eog,...\n');
  fprintf(fid,'  ''prestim_dur'',prestim_dur,...\n');
  fprintf(fid,'  ''poststim_dur'',poststim_dur,...\n');
  fprintf(fid,'  ''bandpass_low_cf'',bandpass_low_cf,...\n');
  fprintf(fid,'  ''bandpass_low_tb'',bandpass_low_tb,...\n');
  fprintf(fid,'  ''bandpass_high_cf'',bandpass_high_cf,...\n');
  fprintf(fid,'  ''bandpass_high_tb'',bandpass_high_tb,...\n');
  fprintf(fid,'  ''dsfact'',dsfact,...\n');
  fprintf(fid,'  ''baseline_start'',baseline_start,...\n');
  fprintf(fid,'  ''baseline_end'',baseline_end,...\n');
  fprintf(fid,'  ''noise_start'',baseline_start,...\n');
  fprintf(fid,'  ''noise_end'',baseline_end,...\n');
  fprintf(fid,'  ''null_event'',null_event,...\n');
  fprintf(fid,'  ''post_dsfact'',post_dsfact,...\n');
  fprintf(fid,'  ''code_excl'',code_excl,...\n');
  fprintf(fid,'  ''time_excl_pre'',time_excl_pre,...\n');
  fprintf(fid,'  ''time_excl_post'',time_excl_post,...\n');
  fprintf(fid,'  ''bandpass_flag'',bandpass_flag,...\n');
  fprintf(fid,'  ''detrend_flag'',detrend_flag,...\n');
  fprintf(fid,'  ''baseline_flag'',baseline_flag,...\n');
  fprintf(fid,'  ''post_subnull_flag'',post_subnull_flag,...\n');
  fprintf(fid,'  ''post_bandpass_flag'',post_bandpass_flag,...\n');
  fprintf(fid,'  ''post_detrend_flag'',post_detrend_flag,...\n');
  fprintf(fid,'  ''post_baseline_flag'',post_baseline_flag,...\n');
  fprintf(fid,'  ''saveepochs_flag'',saveepochs_flag,...\n');
  fprintf(fid,'  ''browseraw'',browseraw);\n');
  fclose(fid);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  fprintf('%s: running automatically generated processing script on %d fif file(s)\n',...
    mfilename, length(flist));
end;

% 'run' uses eval, which makes code un-debuggable through the debugger.
% we can avoid that, so let's do so.
if (strcmp(script_name,'run_process_data')~=1)
  error('%s: (PROGRAMMING ERROR): processing script name was changed, but name of called script was not changed.');
else
  old_wd=pwd;
  cd(parms.outdir);
  % must correspond to file created above:
  run_process_data
  cd(old_wd);
end;
