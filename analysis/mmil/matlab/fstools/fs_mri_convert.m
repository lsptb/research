function fs_mri_convert(fname_in,fname_out,varargin )
%function fs_mri_convert(fname_in,fname_out,varargin )
%
% Purpose: call freesurfer's mri_convert for volume file format conversions
%
% Usage:
%  fs_mri_convert(fname_in,fname_out,'key1', value1,...); 
%
% Required Parameters:
%   fname_in: full or relative path name of input file
%   fname_out: full or relative path name of output file
%
% Optional parameters:
%  'out_orient' - output orientation (e.g. 'LPS', 'RAS', etc.)
%     if empty or omitted, will keep original orientation
%    {default = []}
%  'options' - other options (see mri_convert --help)
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Created:  07/19/08 Don Hagler
% Last Mod: 06/15/09 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, { ...
  'out_orient',[],[],...
  'options',[],[],...
  'forceflag',false,sort([false true]),...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(fname_in,fname_out)
  error('fname_in and fname_out must be different');
end;
if ~exist(fname_in,'file')
  error('file %s not found',fname_in);
end;

if ~exist(fname_out,'file') || parms.forceflag
  fprintf('%s: converting %s to %s...\n',mfilename,fname_in,fname_out);
  cmd = 'mri_convert';
  if ~isempty(parms.out_orient)
    cmd = [cmd ' --out_orientation ' upper(parms.out_orient)];
  end;
  if ~isempty(parms.options)
    cmd = [cmd ' ' parms.options];
  end;
  cmd = [cmd ' ' fname_in ' ' fname_out];
  [status,msg] = unix(cmd);
  if status || ~exist(fname_out,'file')
    error('failed to convert %s to %s:\n%s\n%s\n',...
      fname_in,fname_out,cmd,msg);
  end;
end;

