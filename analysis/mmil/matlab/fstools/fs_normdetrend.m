function fname_out = fs_normdetrend(fname,varargin)
%function fname_out = fs_normdetrend(fname,varargin)
%
% Purpose: normalize and detrend timeseries volume
%   input mgh file, get back mgh file as output
%
% Usage:
%  fs_normdetrend(fname,'key1', value1,...); 
%
% Required Parameters:
%   fname: full or relative path name of input mgh file
%
% Optional parameters:
%  'fname_out' - full or relative path name of output mgh file
%    If empty, will create output name based on input fname
%    If specified, 'outdir' is ignored
%    {default: []}
%  'outdir' - output directory
%    If empty, will write output file to directory containing input fname
%    {default: []}
%  'skipTRs' - number of initial repetitions to omit from analysis
%    {default: 0}
%  'norm_flag' - [0|1] whether to normalize input timeseries
%    by mean for each voxel (new mean = 100) before doing Fourier calculations
%    {default: 1}
%  'detrend' - [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    {default: 2}
%  'out_ext' - output file extension ('.mgh' or '.mgz')
%    if empty, will use input file extension
%    {default: []}
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Output:
%   fname_out: output file name
%
% Created:  09/01/08 Don Hagler (from fs_fourier)
% Last Mod: 09/01/08 Don Hagler
%

%% todo: fname_out as input parameter?

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fnames_fstats = [];

if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_out',[],[],...
  'outdir',[],[],...
  'skipTRs',0,[0 Inf],...
  'norm_flag',true,[false true],...
  'detrend',2,[0:2],...
  'thresh',10,[0 Inf],...
  'out_ext','.mgh',{'.mgh','.mgz'},...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct output file name

if isempty(parms.fname_out)
  [fpath,fstem,fext] = fileparts(fname);
  if ~ismember(fext,{'.mgh','.mgz'})
    error('input file must be mgh or mgz file type (has %s extension)',...
      fext);
  end;
  if isempty(parms.out_ext), parms.out_ext = fext; end;
  if isempty(parms.outdir), parms.outdir = fpath; end;
  fname_out = sprintf('%s/%s_normdetrend%s',parms.outdir,fstem,parms.out_ext);
else
  fname_out = parms.fname_out;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(fname_out) || parms.forceflag
  % load input mgh volume
  [vol,M] = fs_load_mgh(fname);
  [nx,ny,nz,nframes] = size(vol);

  % strip dummy TRs
  if parms.skipTRs>0 && parms.skipTRs<nframes
    vol = vol(:,:,:,parms.skipTRs+1:end);
    [nx,ny,nz,nframes] = size(vol);
  end;

  if parms.norm_flag || parms.detrend
    % reshape into nvox x ntime
    vec = reshape(vol,[nx*ny*nz,nframes]);

    if parms.norm_flag
      % normalize to mean
      vec_mean = mean(vec,2);
      vec = 100*vec./(vec_mean*ones(1,nframes)+eps);
      vec(find(vec_mean<parms.thresh),:) = 0;
    end;

    if parms.detrend
      t = [0:nframes-1]';
      switch parms.detrend
        case 1
          % linear detrend
          A= [t ones(size(t))];
        case 2
          % quadratic detrend
          A = [t t.^2 ones(size(t))];
      end;
      beta = (pinv(A)*vec')';
      trend = beta*A';
      vec = vec - trend;
    end;

    % reshape back to volume
    vol = reshape(vec,size(vol));
  end;

  fs_save_mgh(vol,fname_out,M);
end;

