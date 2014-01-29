function ts_freqanalysis( indata, varargin )
%function ts_timefreq( indata, varargin )
%
% Purpose: 
%   Do time-frequency analysis on epoch data.
%
%
% Required Params (ONE OF THE FOLLOWING):
%
%   indata:    epoch_data structure
%
% Optional Params
%
%  tool: fieldtrip (which calculation / program to use)
%              default: fieldtrip
%  outdir: output directory to write timefreq_data file
%  
%  prefix:  prefix to prepend to output file
%
%  See ts_timefreq_fieldtrip and ts_timefreq_mao
%  for information on analysis input arguments
%
%
%  Created by:        Ben Cipollini  07/10/2007
%  Last Modified by:  Rajan Patel    12/20/2007
%
% See also: ts_freqanalysis_fieldtrip
%
% TO DO: Allow different inputs: epoch file, fif files, combination of
%        stuff, ect.

if (~mmil_check_nargs(nargin, 1))
    return;
end;
  
  
  % First, scrub and validate the params.
  parms           = mmil_args2parms( varargin, ...
                                     { 'tool', 'fieldtrip', {'fieldtrip'}, ...
                                       ... %I/O args
                                       'outdir',      '',       [],...
                                       'prefix', '', [], ...
                                       ... %user prefs
                                       'verbose',      true,     sort([false true]), ...
                                       'logfile',      [],       [],...
                                       'logfid',       1,        [] ...
                                       
                                     }, ...
                                     false );
 
  matfile = fullfile(parms.outdir,'matfiles', [parms.prefix 'timefreq_data.mat']);
                                   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %   Run the analysis
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (exist(matfile, 'file'))
    mmil_logstr(parms, 'NOTICE: not re-analyzing existing timefreq data.');
    
  else
    % remove known params
    args = mmil_parms2args(rmfield(parms, {'tool', 'outdir', 'prefix'}));
    switch (parms.tool)
      case 'fieldtrip'
        timefreq_data = ts_freqanalysis_fieldtrip(indata, varargin{:});
      otherwise
        mmil_error(parms, 'Analysis tool NYI: %s', parms.tool);
    end;
  end;    
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %  Write the output files.
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (exist('timefreq_data', 'var'))
    mmil_logstr(parms, 'Saving output mat files.');

    % First, save the full file
    save(matfile, 'timefreq_data');

  end;
  
  mmil_logstr(parms, 'Done.');
  
