function fs_erode_mask(fname_in,fname_out,sigma,forceflag)
%function fs_erode_mask(fname_in,fname_out,sigma,forceflag)
%
% Purpose: Creating eroded ROIs to avoid boundaries
%    and prevent partial voluming in ROI analysis
% Description: input file is binarized, smoothed
%    and then thresholded, creating an ROI with an outer boundary some
%    distance inside the "true" boundary
%
% Required Input:
%   fname_in: full path to input mask volume
%   fname_out: full path to output eroded mask volume
%
% Optional Input:
%   sigma: width (in voxels) of 3D isotropic Gaussian blurring kernel
%     {default = 1}
%   forceflag: [0|1] toggle overwrite of existing output file
%     {default = 0}
%
% created:  04/01/09 by Don Hagler
% last mod: 04/01/09 by Don Hagler
%

if (~mmil_check_nargs(nargin, 2)), return; end;

smf = 10^-2;
results = [];

if ~exist('sigma','var') | isempty(sigma), sigma=1; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;

if ~exist(fname_in,'file')
  error('input file %s not found',fname_in);
end;
if exist(fname_out,'file') & ~forceflag, return; end;

% load input segmentation volume
fprintf('%s: loading input mask volume %s...\n',mfilename,fname_in);
[vol,M,mr_parms,volsz] = fs_load_mgh(fname_in);

% binarize
vol(vol>0) = 1;
vol(vol<=smf) = 0;

fprintf('%s: eroding mask...\n',mfilename);
vol = double(real(smooth3d(single(vol),sigma,sigma,sigma)));
vol(vol>=1-smf) = 1;
vol(vol<1-smf) = 0;

fprintf('%s: saving output mask volume %s...\n',mfilename,fname_out);
fs_save_mgh(vol,fname_out,M,mr_parms);

