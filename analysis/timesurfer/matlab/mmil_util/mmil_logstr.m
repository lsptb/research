function mmil_logstr(varargin)
 %function mmil_logstr(varargin)
%
% Purpose:
%   allow powerful string logging with consitent formatting:
%     - always output mfilename and newline
%     - allow optional "parms" input to specify:
%         logfile : file where you want to log to 
%         OR       
%         logfid  : file ID where to log to (default: 1, or stdout)
%         verbose : whether to output or not.
%   
%   varargin can contain both a string and parameters to print to the
%   string, a-la fprintf
%
%
% Thus, your parameter parsing would simply need to define these
%   fields, and you get powerful logging!
%
% Example parameter parsing:
%
%  parms   = mmil_args2parms( varargin, ...
%                            { 'verbose',        true,     sort([false true]), ...
%                              'logfile',        [],       [],...
%                              'logfid',         1,        [] 
%                            });
%
% Input Parameters:
%   parms    : (optional) input parameter object with above fields defined
%   str      : string to "log"
%   varargin : arguments to "print" into the format string str
%
% See also: fprintf
%
%
% Created By:       Ben Cipollini on 07/15/2007
% Last Modified By: Ben Cipollini on 08/15/2007]

  if (~mmil_check_nargs(nargin, 1))
    return;
  end;
  
  logfid = 1;
  
  % no parms specified
  if (~isstruct(varargin{1}))
    logfid = 1;
    
  % parms object tells us whether to print or not.
  else
    mmil_check_nargs(nargin, 2, 'Not enough input args; must specify a string to log.');
    
    parms = varargin{1};
    
    % Shouldn't be logging anything.
    if (isfield(parms,'verbose') && ~parms.verbose)
      return;
    end;
    
	  if ((isfield(parms,'logfile') && ~isempty(parms.logfile)) ...
       && (~isfield(parms,'logfid') || parms.logfid==1))
  	  logfid = fopen(parms.logfile, 'a');
	    we_opened_logfid = true;
      
    elseif (isfield(parms,'logfid'))
      logfid = parms.logfid;
      
    else
	    logfid = 1;
    end;
    
    % shift args so output string is 1
    varargin = varargin(2:end);
  end;

  % if we're here, we should print
  M = dbstack;
    
  if (length(M)<2)
    caller = '';
  else
    caller = M(2).name;
    if (strcmp(caller, 'logstr')==1 ||strcmp(caller,'mmil_error')==1)
      if (length(M) < 3)
        caller = '';
      else
        caller = M(3).name;
      end;
    end;
  end;
  
  if (strcmp(caller, '')==1)
    prefix = '';
  else
    prefix = sprintf('%s: ', caller);
  end;
  
  % prep the args
  formatstr = varargin{1};
  if (length(varargin) > 1)
    args = {varargin{2:end}};
  else
    args = {};
  end;

  % output the info
  outstr = sprintf(formatstr, args{:});
  fprintf(logfid, [prefix outstr '\n']);
  
  % if necessary, close the logfid (e.g. if we opened it)
  if (exist('we_opened_logfid', 'var') && we_opened_logfid)
    fclose(logfid);
  end;
