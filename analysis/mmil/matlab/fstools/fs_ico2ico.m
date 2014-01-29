function outfile = fs_surf2ico(infile,varargin )
%function outfile = fs_surf2ico(infile,varargin )
%
% Purpose: resample a surface stats file (mgh format)
%   from one icoshedral order to another
%
% Usage:
%  fs_ico2ico(infile,subject,'key1', value1,...); 
%
% Required Parameters:
%   infile:  full pathname of input file
%
% Optional parameters:
%  'ico_subj': name of icosahedral freesurfer subject
%     {default: 'fsaverage'}
%  'outfile'  - full path of output file
%     if ommitted, outfile will be constructed from infile
%     {default: []}
%  'sparsesmooth' - number of sparse smoothing steps
%    {default = 0}
%    note: sparse smoothing is a fast way to fill
%          in gaps between sparse vertices
%  'postsmooth' - number of normal smoothing steps
%    {default = 0}
%    note: postsmoothing is additional nearest-neighbor average
%          smoothing applied after sparse smoothing
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the avgerage subject directory
%    {default = $FREESURFER_HOME/subjects}
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
%
% NOTES:
%   Input file must have number of values corresponding to one of the
%     following icosahedral orders:
%
%   Order  Number of Vertices
%     1              42
%     2             162
%     3             642
%     4            2562
%     5           10242
%     6           40962
%     7          163842
%
%  Output file will have ico order corresponding to ico_subj  
%
% Created:       09/20/2008 by Don Hagler
% Last Modified: 02/23/2009 by A.I. 
%
% 02/23/2009 - A.I.: changed vol to invals, sparsesmooth to parms.sparsesmooth

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ns = [42 162 642 2562 10242 40962 163842];
outfile = [];

if (~mmil_check_nargs(nargin, 1)), return; end;
parms = mmil_args2parms(varargin, {...
  'ico_subj','fsaverage',[],...
  'outfile',[],[],...
  'sparsesmooth',0,[0 1000],...
  'postsmooth',0,[0 1000],...
  'subjdir',[],[],...
  'forceflag',false,[false true],...
... % hidden
  'hemi','lh',{'lh','rh'},...
});

if isempty(parms.subjdir)
  tmp = getenv('FREESURFER_HOME');
  if isempty(tmp)
    error('FREESURFER_HOME environment variable not set');
  end;
  parms.subjdir = [tmp '/subjects'];
end;

fullsubj = [parms.subjdir '/' parms.ico_subj];
if ~exist(fullsubj,'dir')
  error('ico subj %s not found',fullsubj);
end;

surf = fs_load_subj(parms.ico_subj,parms.hemi,[],1,parms.subjdir); % read nverts only
ico_out = find(surf.nverts==Ns);
if isempty(ico_out)
  error('ico subj %s has non-ico number of vertices (%d)',fullsubj,surf.nverts);
end;

% set outfile
if isempty(parms.outfile)
  [indir,inname,inext] = fileparts(infile);
  pat = sprintf('(?<instem>.+)-(?<hemi>[lr]h)');
  n = regexp(inname,pat,'names');
  if isempty(n)
    fprintf('%s: WARNING: input file name %s missing hemi tag (lh or rh) at end of name',inname);
    instem = inname;
    hemi = [];
  else
    instem = n.instem;
    hemi = n.hemi;
    parms.hemi = hemi;
  end;
  parms.outfile = sprintf('%s/%s-ico%d',indir,instem,ico_out);
  if parms.sparsesmooth>0
    parms.outfile = sprintf('%s-spsm%d',parms.outfile,parms.sparsesmooth);
  end;
  if parms.postsmooth>0
    parms.outfile = sprintf('%s-sm%d',parms.outfile,parms.postsmooth);
  end;
  if ~isempty(hemi)
    parms.outfile = sprintf('%s-%s',parms.outfile,hemi);
  end;
  parms.outfile = [parms.outfile '.mgh'];
end;

if exist(parms.outfile) && ~parms.forceflag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load input file
invals = fs_load_mgh(infile);
if size(invals,2)>1 || size(invals,3)>1
  error('input file %s appears to contain volume data',infile);
end;
nverts = size(invals,1);
nframes = size(invals,4);
ico_in = find(nverts==Ns);
if isempty(ico_in)
  error('input file has non-ico number of vertices (%d)',nverts);
end;
invals = reshape(invals,[nverts,nframes]);

% load ico_subj
if parms.sparsesmooth || parms.postsmooth
  surf = fs_load_subj(parms.ico_subj,parms.hemi,[],[],parms.subjdir);
end;

% create output
fprintf('%s: resampling from ico%d to ico%d...\n',mfilename,ico_in,ico_out);
outvals = zeros(surf.nverts,nframes);
if ico_in<ico_out
  outvals(1:Ns(ico_in),:) = invals;
else
  outvals = invals(1:Ns(ico_out),:);
end;

% apply smoothing
if parms.sparsesmooth>0 || parms.postsmooth>0
  fprintf('%s: smoothing %d frames...\n',mfilename,nframes);
  for f=1:nframes
    tmpvals = outvals(:,f);
    if parms.sparsesmooth>0
      tmpvals = fs_smooth_sparse(surf,tmpvals,parms.sparsesmooth);
    end;
    if parms.postsmooth>0
      tmpvals = fs_smooth(surf,tmpvals,parms.postsmooth);
    end;
    outvals(:,f) = tmpvals;
  end;
end

% save output
outvals = reshape(outvals,[surf.nverts,1,1,nframes]);
fs_save_mgh(outvals,parms.outfile);
outfile = parms.outfile;

