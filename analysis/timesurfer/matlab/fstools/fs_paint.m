function fnames_out = fs_paint(subj,funcname,varargin);
%function fnames_out = fs_paint(subj,funcname,[options]);
%
% Purpose:
%   a wrapper around freesurfer binaries for sampling
%   volume data to surface, sampling to sphere, and smoothing
%
% Usage:
%  fs_paint(subj,funcname,'key1', value1,...);
%
% Required parameters:
%  subj is a string specifying the subject name
%  funcname is a string specifying the full or relative path
%    of the functional volume (must be mgh format)
%    (can be empty if 'meas' is supplied -- see below)
%
% Optional parameters:
%  'hemi' - should be either 'lh' for left hemisphere or 'rh' for right hemi
%    {default = both}
%  'meas' - surface measure such as 'thickness', or 'avgarea'
%     Note: for 'avgarea', must be an "avg" subject
%    {default = []}
%  'outstem'  - output file stem (omit extension, hemi)
%    {default = funcname or meas}
%  'infix'  - add extra suffix to outstem before hemi and extension
%    {default = []}
%  'outtype' - output file type ('w' or 'mgh')
%    {default: 'mgh'}
%  'intype' - input file type (e.g. 'mgh', 'analyze', 'bfloat')
%    {default: 'mgh'}
%  'regfile' - register.dat file containing 4x4 registration matrix
%    {default: []}
%  'regmat' - 4x4 registration matrix (ignored if 'regfile' is specified)
%    {default: identity matrix}
%  'inplane' - inplane resolution (mm) (ignored if 'regfile' is specified)
%    {default: 1}
%  'slicethick' - slice thickness (mm) (ignored if 'regfile' is specified)
%    {default: 1}
%  'sphere_flag' - [0|1] whether to sample to icosohedral sphere
%    {default: 0}
%  'projfrac_flag' - [0|1] whether to use projdist (0) or projfract (1)
%    {default: 0}
%  'projdist' - distance (mm) to project along surface vertex normal
%    {default: 1}
%  'projfrac' - fractional distance to project along surface vertex normal
%    relative to cortical thickness
%    {default: 0.5}
%  'smoothsteps' - smoothing steps on surface after painting
%    {default: 0}
%  'sphsmoothsteps' - smoothing steps on spherical surface
%    {default: 0}
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'surfname' - surface to paint onto
%    {default: white}
%  'mask_midbrain_flag' - [0|1] toggle mask out thalamus, corpus collosum
%    {default: 1}
%  'overwrite_flag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Output:
%   fnames_out: cell array of output file names (e.g. left and right)
%     with multiple steps (e.g. sphere, smoothing, mask), only the final
%     output files are included in fnames_out
%
% created:  10/19/06 by Don Hagler
% last mod: 04/01/09 by Don Hagler

%% todo: surf data as input for resampling to sphere, smoothing only
%% todo: use mmil_args2parms
%% todo: optionally resample to T1 before painting

%% todo: output as mgz
%  'outtype' - output file type ('w','mgh' or 'mgz')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set defaults
DEFAULT_INTYPE = 'mgh';
DEFAULT_SURFNAME = 'white';
DEFAULT_OUTTYPE = 'mgh';
DEFAULT_PROJFRAC = 0.5;
DEFAULT_PROJDIST = 1;
mask_roinums = [1,5]; % 'unknown' and 'corpuscallosum'

fnames_out = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

if ~isfield(g,'hemi'), g.hemi = []; end;
if ~isfield(g,'hemilist'), g.hemilist = []; end;
if ~isfield(g,'outstem'), g.outstem = []; end;
if ~isfield(g,'infix'), g.infix = []; end;
if ~isfield(g,'outtype'), g.outtype = DEFAULT_OUTTYPE; end;
if ~isfield(g,'regfile'), g.regfile = []; end;
if ~isfield(g,'regmat'), g.regmat = eye(4); end;
if ~isfield(g,'inplane'), g.inplane = 1; end;
if ~isfield(g,'slicethick'), g.slicethick = 1; end;
if ~isfield(g,'sphere_flag'), g.sphere_flag = 0; end;
if ~isfield(g,'subjdir'), g.subjdir = []; end;
if ~isfield(g,'surfname'), g.surfname = DEFAULT_SURFNAME; end;
if ~isfield(g,'projfrac'), g.projfrac = DEFAULT_PROJFRAC; end;
if ~isfield(g,'projdist'), g.projdist = DEFAULT_PROJDIST; end;
if ~isfield(g,'projfrac_flag'), g.projfrac_flag = 0; end;
if ~isfield(g,'smoothsteps'), g.smoothsteps = 0; end;
if ~isfield(g,'sphsmoothsteps'), g.sphsmoothsteps = 0; end;
if ~isfield(g,'intype'), g.intype = DEFAULT_INTYPE; end;
if ~isfield(g,'meas'), g.meas = []; end;
if ~isfield(g,'mask_midbrain_flag'), g.mask_midbrain_flag = 1; end;
if ~isfield(g,'overwrite_flag'), g.overwrite_flag = 0; end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'hemi' 'hemilist'...
         'outstem' 'outtype' 'regmat' 'regfile' 'inplane' 'slicethick'...
         'sphere_flag' 'subjdir' 'surfname' 'projfrac' 'intype'...
         'overwrite_flag' 'infix' 'projdist' 'projfrac_flag'...
         'smoothsteps' 'sphsmoothsteps' 'meas' 'mask_midbrain_flag'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
hemi = g.hemi;
hemilist = g.hemilist;
outstem = g.outstem;
infix = g.infix;
outtype = g.outtype;
regfile = g.regfile;
regmat = g.regmat;
inplane = g.inplane;
slicethick = g.slicethick;
sphere_flag = g.sphere_flag;
subjdir = g.subjdir;
surfname = g.surfname;
projfrac = g.projfrac;
projdist = g.projdist;
projfrac_flag = g.projfrac_flag;
smoothsteps = g.smoothsteps;
sphsmoothsteps = g.sphsmoothsteps;
intype = g.intype;
meas = g.meas;
mask_midbrain_flag = g.mask_midbrain_flag;
overwrite_flag = g.overwrite_flag;
clear g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parameters
if nargin<2, help(mfilename); return; end;  

if ~exist('subjdir','var'), subjdir = []; end;
if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as an environment variable');
  end;
end;

if isempty(hemilist)
  if isempty(hemi)
    hemilist = {'lh','rh'};
  else
    hemilist = {hemi};
  end;
end;
if any(~ismember(hemilist,{'lh','rh'}))
  error('hemi must be lh or rh',hemi);
end;

if isempty(meas)
  if ~exist(funcname,'file')
    error('functional volume file %s not found',funcname);
  end;
else
  if ~ismember(meas,{'thickness','avgarea'})
    error('unsupported meas %s',meas);
  end;
end;

for h=1:length(hemilist)
  hemi = hemilist{h};
  surffile = sprintf('%s/%s/surf/%s.%s',subjdir,subj,hemi,surfname);
  if ~exist(surffile,'file')
    error('surface file %s not found',surffile);
  end
end;
setenv('SUBJECTS_DIR',subjdir);

if ~isempty(regfile) & ~exist(regfile,'file')
  error('regfile %s not found',regfile);
elseif any(size(regmat)~=[4 4])
  error('regmat is wrong size');
end;

if isempty(outstem)
  if isempty(meas)
    [tmp_path,tmp_fstem,tmp_ext,tmp_ver] = fileparts(funcname);
    outstem = [tmp_path,'/',tmp_fstem];
  else
    outstem = sprintf('%s/%s/stats/%s',...
      subjdir,subj,meas);
  end;
end;
if ~strcmp(outstem(1),'/')
  outstem = sprintf('./%s',outstem);
end;
[outdir,tmp_fstem] = fileparts(outstem);


%if ~ismember(outtype,{'w','mgh','mgz'})
%  error('outtype must be w, mgh, or mgz (is %s)\n',outtype);
if ~ismember(outtype,{'w','mgh'})
  error('outtype must be w or mgh (is %s)\n',outtype);
end;

if ~ismember(intype,{'mgh','bfloat' 'analyze'})
  error('intype must be mgh, analyze, or bfloat (is %s)',intype);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(meas) % i.e. painting from volume to surface
%% todo: and not surf file
  if ~isempty(regfile)
   [regmat,tmp_subj,inplane,slicethick] = fs_read_regdat(regfile);
  end;

  % create tmp register.dat file (in case regfile not specified or subj is different)
  tmp_regfile = sprintf('%s/%s-tmpreg.dat',outdir,hemi);

  fs_write_regdat(tmp_regfile,...
    'M',regmat,...
    'subj',subj,...
    'inplane',inplane,...
    'slicethick',slicethick,...
    'forceflag',1);
end;

for h=1:length(hemilist)
  hemi = hemilist{h};

  % sample from volume to surface
  if isempty(infix)
    prefix_out = outstem;
  else
    prefix_out = sprintf('%s%s',outstem,infix);
  end;
  outfile = sprintf('%s-%s.%s',prefix_out,hemi,outtype);
  fnames_out{h} = outfile;
  if ~exist(outfile,'file') | overwrite_flag
    % create mri_vol2surf cmd string
    if ~isempty(meas)
      fprintf('%s: painting cortical %s\n',mfilename,meas);
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_surf2surf --srcsubject %s --hemi %s',...
        subjdir,subj,hemi);
      switch meas
        case 'thickness'
          cmd = sprintf('%s --srcsurfval thickness --src_type curv',cmd);
        case 'avgarea'
          cmd = sprintf('%s --sval %s.white.avg.area.mgh',cmd);
      end;
      cmd = sprintf('%s --tval %s --trgsubject %s',cmd,outfile,subj);
    else
      fprintf('%s: painting func volume from %s for hemi %s\n',mfilename,funcname,hemi);
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_vol2surf --src %s --src_type %s --srcreg %s --hemi %s',...
        subjdir,funcname,intype,tmp_regfile,hemi);
      cmd = sprintf('%s --surf %s --out %s --out_type %s --trgsubject %s',...
        cmd,surfname,outfile,outtype,subj);
      %% todo: can surf data be saved as mgz?
      if ~projfrac_flag
        cmd = sprintf('%s --projdist %0.3f',...
          cmd,projdist);
      else
        cmd = sprintf('%s --projfrac %0.3f',...
          cmd,projfrac);
      end;
    end;
    cmd = sprintf('%s --noreshape',cmd);
    % run cmd
    [status,result]=unix(cmd);
    if status
      error('cmd %s failed:\n%s',cmd,result);
    end;
  end;

  % mask out midbrain
  if mask_midbrain_flag
    prefix_in = prefix_out;
    prefix_out = sprintf('%s-mbmask',prefix_in);
    infile = sprintf('%s-%s.%s',prefix_in,hemi,outtype);
    outfile = sprintf('%s-%s.%s',prefix_out,hemi,outtype);
    fnames_out{h} = outfile;
%% todo: if outtype is w, then need to load file and mask surfstats
    fs_mask_surfmgh_with_aparc(subj,hemi,infile,outfile,subjdir,mask_roinums);
  end;

  % smooth on surface
  if smoothsteps
    prefix_in = prefix_out;
    prefix_out = sprintf('%s-sm%d',prefix_in,smoothsteps);
    infile = sprintf('%s-%s.%s',prefix_in,hemi,outtype);
    outfile = sprintf('%s-%s.%s',prefix_out,hemi,outtype);
    fnames_out{h} = outfile;
    if ~exist(outfile,'file') | overwrite_flag
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_surf2surf --s %s --hemi %s',...
        subjdir,subj,hemi);
      cmd = sprintf('%s --sval %s --tval %s',...
        cmd,infile,outfile);
      cmd = sprintf('%s --sfmt %s --nsmooth-out %d',...
        cmd,outtype,smoothsteps);
      cmd = sprintf('%s --noreshape',cmd);
      % run cmd
      [status,result]=unix(cmd);
      if status
        error('cmd %s failed:\n%s',cmd,result);
      end;
    end;
  end;

  % sample to sphere
  if sphere_flag
    prefix_in = prefix_out;
    if sphsmoothsteps
      prefix_out = sprintf('%s-sphere-sm%d',prefix_in,sphsmoothsteps);
    else
      prefix_out = sprintf('%s-sphere',prefix_in);
    end;    
    infile = sprintf('%s-%s.%s',prefix_in,hemi,outtype);
    outfile = sprintf('%s-%s.%s',prefix_out,hemi,outtype);
    fnames_out{h} = outfile;
    if ~exist(outfile,'file') | overwrite_flag
      cmd = sprintf('setenv SUBJECTS_DIR %s; mri_surf2surf --srcsubject %s --trgsubject ico --hemi %s',...
        subjdir,subj,hemi);
      cmd = sprintf('%s  --trgicoorder 7 --sval %s --tval %s',...
        cmd,infile,outfile);
      cmd = sprintf('%s --sfmt %s --nsmooth-out %d',...
        cmd,outtype,sphsmoothsteps);
      cmd = sprintf('%s --noreshape',cmd);
      % run cmd
      [status,result]=unix(cmd);
      if status
        error('cmd %s failed:\n%s',cmd,result);
      end;
    end;
  end;
end;

if exist('tmp_regfile','var') & exist(tmp_regfile,'file')
  delete(tmp_regfile);
end;
