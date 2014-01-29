function status=fs_reorient_vol(fname_in,fname_out,orientation,overwrite_flag);
%function status=fs_reorient_vol(fname_in,fname_out,[orientation],[overwrite_flag]);
%
% Required input:
%  fname_in: full or relative path of input file name
%  fname_out: full or relative path of output file name
%
% Optional input:
%  orientation: orientation string
%    {default = 'LPI'}
%  'overwrite_flag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% created: 02/08/06 by Don Hagler
% last modified: 02/08/07 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

status=0;

% check parameters
if nargin<2, help(mfilename); return; end;  

if ~exist('orientation','var'), orientation = []; end;
if isempty(orientation), orientation = 'LPI'; end;
if ~exist('overwrite_flag','var'), overwrite_flag = []; end;
if isempty(overwrite_flag), overwrite_flag = 0; end;

if ~exist(fname_in,'file')
  fprintf('%s: ERROR: input file %s not found\n',mfilename,fname_in);
  return;
end;

orientation = upper(orientation);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_out,'file') | overwrite_flag
  % create mri_convert cmd string
  cmd = sprintf('mri_convert --out_orientation %s %s %s',...
    orientation,fname_in,fname_out);
  % run cmd
  fprintf('%s\n',cmd);
  [status,result]=unix(cmd);
  if status
    fprintf('%s: ERROR:\n',mfilename);
    disp(result);
  end;
end;

