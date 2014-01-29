function fs_erode_aseg(fname_in,fname_out,sigma,forceflag)
%function fs_erode_aseg(fname_in,fname_out,sigma,forceflag)
%
% Purpose: Creating eroded segmentation ROIs to avoid boundaries
%    and prevent partial voluming in ROI analysis
% Description: Each ROI in the input segmentation file is smoothed
%    and then thresholded, creating an ROI with an outer boundary some
%    distance inside the "true" boundary
%
% Required Input:
%   fname_in: full path to input segmentation volume (e.g. aseg.mgz)
%   fname_out: full path to output segmentation volume (e.g. aseg_eroded.mgz)
%
% Optional Input:
%   sigma: width (in voxels) of 3D isotropic Gaussian blurring kernel
%     {default = 1}
%   forceflag: [0|1] toggle overwrite of existing output file
%     {default = 0}
%
% created:  05/14/07 by Don Hagler
% last mod: 05/14/07 by Don Hagler
%

if nargin<2
  help(mfilename);
  return;
end;  

smf = 10^-2;
results = [];

if ~exist('sigma','var') | isempty(sigma), sigma=1; end;
if ~exist('forceflag','var') | isempty(forceflag), forceflag=0; end;

if ~exist(fname_in,'file')
  fprintf('%s: ERROR: input egmentation volume file %s not found\n',...
    mfilename,fname_in);
  return;
end;
if exist(fname_out,'file') & ~forceflag, return; end;

% load input segmentation volume
fprintf('%s: loading input segmentation volume %s...\n',mfilename,fname_in);
[segvol,M,mr_parms,volsz] = fs_load_mgh(fname_in);

% find unique roi numbers in segvol
segroicodes = unique(segvol(:));
segroicodes = segroicodes(find(segroicodes>0));
nsegrois = length(segroicodes);

fprintf('%s: eroding %d ROIs...\n',mfilename,nsegrois);
segvol2 = zeros(volsz);
for i=1:nsegrois
  fprintf('%s: eroding ROI %d...\n',mfilename,i);
  tmp_vol = zeros(volsz);
  tmp_vol(segvol==segroicodes(i)) = 1.0;
  tmp_vol = double(real(smooth3d(single(tmp_vol),sigma,sigma,sigma)));
  segvol2(tmp_vol>=1-smf) = segroicodes(i);
end;

fprintf('%s: saving output segmentation volume %s...\n',mfilename,fname_out);
fs_save_mgh(segvol2,fname_out,M,mr_parms);

