function varargout = smartunix(cmd, outtype, verbose)
% varargout = smartunix(cmd, [outtype], [verbose])
%
% Purpose:
%   Execute a unix command while addressing the following issues:
%   - easily specify whether output should be shown to screen or not
%   - easily choose which return values to return
%   - return a true/false value if function returned a "succeded" SC (0)
%   - move MATLAB paths to bottom of LD_LIBRARY_PATH stack, to simulate
%     better non-MATLAB shell environment.
%
% Input Arguments:
%   cmd     : UNIX command-line to execute.
%   outtype : string specifying which output argument to return ('result', 'output')
%             type 'result': true/false value indicating if command succeeded
%             type 'output': stdout/stderr output from running cmd
%             {default: result AND output are returned, in that order)
%  verbose  : whether to print cmd output to screen or not
%             {default: false}
%
% Output Arguments:
%   'result' and/or 'output'; see 'outtype' argument above for details.
%
% See also: unix
%
% Created By:       Ben Cipollini on 08/01/2007
% Last Modified By: Ben Cipollini on 08/09/2007

  varargout={};

  % push MATLAB to the back of the path stack.
  path = getenv('LD_LIBRARY_PATH');
  paths = splitstr(path,':');
  mat_paths = find(~isempty2(regexp(paths,'matlab'),true));
  non_mat_paths = setdiff(1:length(paths),mat_paths);

  newpaths = { paths{non_mat_paths} };
  if (isempty(newpaths))
    newpath  = '';
  else
    newpath  = [sprintf('%s:',newpaths{1:end-1}) newpaths{end}];
  end;
  
  % 
  setenv('LD_LIBRARY_PATH',newpath);
  if (exist('verbose','var') && verbose)
    [r,o] = unix(cmd)
  else
    [r,o] = unix(cmd);
  end;    
  setenv('LD_LIBRARY_PATH',path);

  % return the appropriate type
  if (~exist('outtype','var'))
    varargout{end+1}=(r==0);
    varargout{end+1}=o;
  elseif (strcmp(outtype, 'result'))
    varargout{end+1}=(r==0);
  else
    varargout{end+1} = o;
  end;


