function outfile = fs_surf2ico(infile,subject,varargin )
%function outfile = fs_surf2ico(infile,subject,varargin )
%
% Purpose: resample a surface stats file (mgh or w formats) to icosahedral sphere 
%   (wrapper for freesurfer's mri_surf2surf)
%
% Usage:
%  fs_surf2ico(infile,subject,'key1', value1,...); 
%
% Example:
%  fs_surf2ico('/usr/data/surfstats-lh.mgh','bert','lh',...
%              'outfile','/usr/data/surfstats-ico-lh.mgh')
%
% Required Parameters:
%   infile:  full pathname of input file
%   subject: freesurfer recon subject name
%
% Optional parameters:
%  'outfile'  - full path of output file
%     if ommitted, outfile will be constructed from infile
%     {default: []}
%  'outdir'  - directory for output file
%     if ommitted, will use same directory as input file
%     if outfile is supplied, outdir is ignored
%     {default: []}
%  'hemi' - cortical hemisphere (lh or rh)
%     necessary only if input file name does not have hemi tag at end
%       e.g. 'stem-lh.mgh' or 'stem-rh.mgh'
%     if file name does have hemi tag, this parameter will be ignored
%     {default: []}
%  'intype' - input file type ('mgh' or 'w')
%  'outtype' - output file type ('mgh' or 'w')
%  'icolevel' - icosahedron order number:
%               Order  Number of Vertices
%                 0              12
%                 1              42
%                 2             162
%                 3             642
%                 4            2562
%                 5           10242
%                 6           40962
%                 7          163842
%    {default: 7}
%  'smooth_in' - smoothing steps on surface before sampling to ico
%    {default: 0}
%  'smooth_out' - smoothing steps on surface after sampling to ico
%    {default: 0}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%  'verbose' - [0|1] toggle output of mri_surf2surf output
%    {default: 1}
%
% Created:       05/30/2007 by Ben Cipollini
% Last Modified: 12/19/2008 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, ...
        { 'outfile',      '', [],...
          'outdir',       '', [],...            
          'hemi',         '', [{'lh','rh'}],...
          'intype',       'mgh', [{'mgh','w'}],...
          'outtype',      'mgh', [{'mgh','w'}],...
          'subjdir',       getenv('SUBJECTS_DIR'), [], ...
          'smooth_in',  0, [0 Inf],...
          'smooth_out', 0, [0 Inf],...
          'icolevel',      7, [1 7],...
          'forceflag',     false, sort([false true]),...
          'verbose',       true, sort([false true])}...
        );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check the input file name, extract the instem
[indir,inname,inext] = fileparts(infile);
if ~strcmp(['.' parms.intype],inext)
  error('infile %s missing correct file extension for intype %s',infile,parms.intype);
end;
pat = sprintf('(?<instem>.+)-(?<hemi>[lr]h)');
n = regexp(inname,pat,'names');
if isempty(n)
  if isempty(parms.hemi)
    error('input file name %s missing hemi tag (lh or rh) at end of name and hemi parameter not specified',inname);
  else
    fprintf('%s: WARNING: input file name %s missing hemi tag (lh or rh) at end of name',inname);
  end;
  instem = inname;
else
  instem = n.instem;
  parms.hemi = n.hemi;
end;
if (isempty(parms.outdir))
  parms.outdir = indir;
end;
if ~exist(parms.outdir,'dir')
  [success,msg] = mkdir(parms.outdir);
  if ~success
    error('failed to create outdir %s: %s',parms.outdir,msg);
  end;
end;

% construct the output file name
if (isempty(parms.outfile))
  outstem = instem;
  if (parms.smooth_in > 0)
    outstem = sprintf('%s-sm%d',outstem,parms.smooth_in);
  end;
  outstem = sprintf('%s-ico%d',outstem,parms.icolevel);
  if (parms.smooth_out > 0)
    outstem = sprintf('%s-sm%d',outstem,parms.smooth_out);
  end;
  parms.outfile = sprintf('%s-%s.%s',outstem,parms.hemi,parms.outtype);
  parms.outfile = fullfile(parms.outdir,parms.outfile);
end;

if (parms.forceflag || ~exist(parms.outfile,'file'))
  % construct the unix command
  cmd = sprintf('mri_surf2surf --srcsubject %s --trgsubject ico --hemi %s',subject,parms.hemi);
  cmd = sprintf('%s  --trgicoorder %d --sval %s --tval %s',cmd,...
                parms.icolevel,infile,parms.outfile);
  cmd = sprintf('%s --sfmt %s --tfmt %s',cmd,...
                parms.intype,parms.outtype);
  if parms.smooth_in>0
    cmd = sprintf('%s --nsmooth-in %d',cmd,parms.smooth_in);
  end;
  if parms.smooth_out>0
    cmd = sprintf('%s --nsmooth-out %d',cmd,parms.smooth_out);
  end;
  cmd = sprintf('%s --noreshape', cmd);
  if (parms.verbose), display(cmd); end;
  oldsdir     = getenv('SUBJECTS_DIR');
  setenv('SUBJECTS_DIR', parms.subjdir);
  % run the unix command
  [status,result]=unix(cmd);
  % restore the original subjects_dir
  setenv('SUBJECTS_DIR', oldsdir);
  % check for success or failure
  if (status || ~isempty(findstr(result, 'could not open')) || ~exist(parms.outfile,'file'))
    error('%s\n\ncmd:\n%s',result,cmd);
  end;
  if (parms.verbose), display(result); end;
end;

outfile = parms.outfile;


