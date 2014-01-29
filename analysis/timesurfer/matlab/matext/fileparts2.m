function s = fileparts2(f,fp)
%
% Purpose: 
%   Works like FILEPARTS, but allows you to 
%   specify WHICH filepart(s) to return.
%
%   Requested 'parts' are concatenated together
%   'appropriately' (if no 'appropriate' concatenation is
%   found, then the parts are just concatenated in the order
%   that they were requested by param fp)
%
% Input Arguments:
%   f   : path to parse (string)
%   fp  : path parts to return (cell array)
%     'path', 'name', 'ext', or 'versn'
%
% Output Arguments:
%   s   : output string of concatenated parts
%
% See also: fileparts
%
% Created By:       Ben Cipollini on 07/15/2007
% Last Modified By: Don Hagler    on 08/21/2007

  if (~iscell(fp))
    fp = {fp};
  end;

  [pathstr name ext versn]   = fileparts(f);
  s   = '';

  for i=1:length(fp)
    switch (fp{i})
      case {'path','pathstr'}
                      s   = fullfile(s,pathstr);
      case 'name',    s   = [s name];
      case 'ext',     s   = [s ext];
      case 'versn',   s   = [s versn];
      otherwise
        error('%s: unknown filepart: %s', mfilename, fp);
    end;
  end;


