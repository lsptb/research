function results = fs_fourier_avg(fnamelist,varargin)
%function results = fs_fourier_avg(fnamelist,[options])
%
% Usage:
%  results = fs_fourier_avg(fnamelist,'key1', value1,...);
%
% Required Input:
%   fnamelist: cell matrix of input full path file names
%     Each row should correspond to a single scan / subject
%     Must have two columns corresponding to real and imaginary components
%
%     All input files should be mgh format surface or volume files
%       for a single subject or registered to common space
%       (e.g. spherical ico surface)
%
%  Optional Input:
%   'frames' - vector of frame numbers (e.g. multiple time points)
%     If empty, will loop over all frames
%     {default = []}
%   'revflags' - vector of 0's and 1's for each subject (each row of fnamelist)
%     determining whether to set the imaginary components negative
%     (e.g. reverse phase)
%     {default = [0 0 0....]}
%   'phase_offset' - phase angle (in cycles) subtracted from
%       complex valued input before applying revflags
%       (e.g. to correct for hemodynamic delay)
%      {default = 0}
%   'phase_offset_postrev' - phase angle (in cycles) subtracted from
%       complex valued input after applying revflags
%       (e.g. to correct for stimulus delay)
%     {default = 0}
%   'verbose_flag' - [0|1] whether to print status messages
%     {default = 1}
%
%
% Output:
%   results: structure containing cross-scan average
%
% Created:  08/18/08 by Don Hagler
% Last Mod: 08/20/08 by Don Hagler
%

%% todo: optionally calculate complex f-stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'frames',[],[],...
  'phase_offset',0,[-1,1],...
  'phase_offset_postrev',0,[-1,1],...
  'revflags',[],[],...
  'verbose_flag',true,sort([false true]),...
});
results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input

if isempty(fnamelist)
  error('input file name list is empty');
elseif ~iscell(fnamelist)
  fnamelist = {fnamelist};
end;
[nscans,nconds] = size(fnamelist);
if nconds~=2
  error('input file name list must have two columns corresponding to real and imaginary components');
end;

if ~isempty(parms.frames) && min(parms.frames)<1
  error('frame numbers must be > 1');
end;

if parms.verbose_flag
  tic
  fprintf('%s: checking input files...\n',mfilename);
end;
% check files
results.orig_volsz = [];
results.M = [];
for s=1:nscans
  for c=1:nconds
    if isempty(fnamelist{s,c})
      error('fnamelist{%d,%d} is empty',s,c);
    end;
    [vol,M,mrparms,volsz] = fs_load_mgh(fnamelist{s,c},[],[],1); % header-only
    if isempty(parms.frames), parms.frames = [1:volsz(4)]; end;
    if max(parms.frames)>volsz(4)
      error('only %d frames in file %s',volsz(4),fnamelist{s,c});
    end;
    if isempty(results.orig_volsz)
      results.orig_volsz = volsz;
      results.M = M;
    end;
    if any(results.orig_volsz~=volsz)
      error('size of input volume for\n%s (%d,%d,%d,%d)\ndoes not match that of\n%s (%d,%d,%d,%d)',...
        fnamelist{s,c},volsz(1),volsz(2),volsz(3),volsz(4),...
        fnamelist{1,1},results.orig_volsz(1),results.orig_volsz(2),results.orig_volsz(3),results.orig_volsz(4));
    end;
  end;
end;
if parms.verbose_flag, toc; end;

results.nvals = prod(results.orig_volsz(1:3));
results.nframes = length(parms.frames);
results.volsz = [results.orig_volsz(1:3) results.nframes];
results.frames = parms.frames;

if isempty(parms.revflags)
  parms.revflags = zeros(nscans,1);
else
  if length(parms.revflags)~=nscans
    error('number of revflags must match number of rows of fnamelist');
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data and calculate sums
results.n = 0;
results.mean_r = zeros(results.nvals,results.nframes,'single');
results.mean_i = zeros(results.nvals,results.nframes,'single');
for f=1:results.nframes
  frame = parms.frames(f);
  if parms.verbose_flag
    tic
    fprintf('%s: summing data for frame %d of %d...\n',mfilename,frame,results.orig_volsz(4));
  end;
  for s=1:nscans
    vec_r = fs_load_mgh(fnamelist{s,1},[],frame);
    vec_r = reshape(vec_r,[prod(size(vec_r)) 1]);
    vec_i = fs_load_mgh(fnamelist{s,2},[],frame);
    vec_i = reshape(vec_i,[prod(size(vec_i)) 1]);
    if parms.phase_offset~=0
      vec_ampl = hypot(vec_r,vec_i);
      vec_phase = atan2(vec_i,vec_r) - parms.phase_offset*2.0*pi;
      vec_r = vec_ampl.*cos(vec_phase);
      vec_i = vec_ampl.*sin(vec_phase);
    end;
    if parms.revflags(s)
      vec_i = -vec_i;
    end;
    if parms.phase_offset_postrev~=0
      vec_ampl = hypot(vec_r,vec_i);
      vec_phase = atan2(vec_i,vec_r) - parms.phase_offset_postrev*2.0*pi;
      vec_r = vec_ampl.*cos(vec_phase);
      vec_i = vec_ampl.*sin(vec_phase);
    end;
    if f==1
      results.n = results.n + 1;
    end;
    results.mean_r(:,frame) = results.mean_r(:,frame) + vec_r;
    results.mean_i(:,frame) = results.mean_i(:,frame) + vec_i;
  end;
  if parms.verbose_flag, toc; end;
end;

% calculate mean
if parms.verbose_flag
  tic
  fprintf('%s: calculating means...\n',mfilename);
end;
results.mean_r = reshape(results.mean_r,results.volsz);
results.mean_r = results.mean_r / (eps+results.n);
results.mean_i = reshape(results.mean_i,results.volsz);
results.mean_i = results.mean_i / (eps+results.n);
if parms.verbose_flag, toc; end;

%keyboard

%% todo: calculate complex f-stats
