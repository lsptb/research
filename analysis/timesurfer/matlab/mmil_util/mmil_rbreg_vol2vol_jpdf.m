function [M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_jpdf(vol1,vol2,type1,type2,volmask,...
  maskflag,initflag,interpm,mask_smoothing,Mreg0,scales)
%function [M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_jpdf(vol1,vol2,[type1],[type2],...
%  [volmask],[maskflag],[initflag],[interpm],[mask_smoothing],[Mreg0],[scales])
%
% Required Input:
%   vol1: 3D volume of reference scan (e.g. structural scan)
%   vol2: 3D volume of scan to be registered (e.g. first frame of functional scan)
%
% Optional Input:
%   type1: image type for vol1 (reference image)
%     (e.g. 'MPR','FLASHHI','FLASHLO', or 'FSNU' (i.e. Freesurer nu.mgz))
%     {default = 'MPR'} 
%   type2: image type for vol2
%     (e.g.'BOLD','DTI','GECAL','PET','FLASHLO', or 'FLASHHI')
%     {default = 'LOFA'}
%   volmask: mask for vol1 (e.g. brain mask)
%     {default = []}
%   maskflag: [0|1] toggle whether to use brain mask for vol1
%     if maskflag=1 and volmask=[], will
%       automatically generate one by registering to atlas
%       (will only work if vol1 is T1-weighted)
%     (this option is ignored -- mask is used no matter what)
%     {default = 1}
%   initflag: [0|1|2] indicate desired behavior for initalizing jpdf
%     0: if jpdf does not already exist, exit with error message
%     1: if jpdf does not already exist, initialize it
%        (this requires that vol1 and vol2 are already well-aligned)
%     2: use existing jpdf if it exists (otherwise initialize)
%        and reinitialize with values from vol1 and vol2 after registration
%     {default = 0}
%   interpm:  0:Nearest 1:Linear (default) 2:Cubic
%             3:Key's spline 4:Cubic spline 5:Hamming Sinc
%   mask_smoothing: vector of x,y,z smoothing sigma (voxels) for mask
%     {default = [3,3,3]}
%   Mreg0: 4x4 matrix defining initial estimate of registration matrix
%     {default = identity}
%   scales: vector of scales (i.e. displacement sizes) for multi-scale search
%     {default:[0 83 49 27 16 9 5 3 2 1]}
%
% Output:
%   M_v1_to_v2: registration matrix from vol1 to vol2
%   volmask: brain mask for vol1
%
% Note: vol2, vol1, and volmask should be in ctx format
%   use vol_ctx=ctx_mgh2ctx(vol,M);
%
% Based on code by Anders Dale
% Created : 12/16/06 by Don Hagler
% Last Mod: 02/16/09 by Don Hagler
%

if nargin<2
  help(mfilename);
  return;
end;

M_v1_to_v2 = [];

if ~exist('type1','var') | isempty(type1), type1 = 'MPR'; end;
if ~exist('type2','var') | isempty(type2), type1 = 'LOFA'; end;
if ~exist('volmask','var'), volmask = []; end;
if ~exist('maskflag','var') | isempty(maskflag), maskflag = 1; end;
if ~exist('initflag','var') | isempty(initflag), initflag = 0; end;
if ~exist('interpm','var') | isempty(interpm), interpm = 1; end;
if ~exist('mask_smoothing','var') | isempty(mask_smoothing)
  mask_smoothing = [3 3 3];
end;
if ~exist('Mreg0','var') | isempty(Mreg0), Mreg0 = eye(4); end;
if ~exist('scales','var') || isempty(scales)
  scales=[0 83 49 27 16 9 5 3 2 1];
end;
if (scales(1)~=0) scales=[0,scales]; end;

type1 = upper(type1);
type2 = upper(type2);

jpdffname = sprintf('/home/mmildev/matlab/atlases/%s_%s_jpdf.mat',type1,type2);

% load jpdf file if it exists
if exist(jpdffname,'file')
  fprintf('%s: using jpdf file %s\n',mfilename,jpdffname);
  load(jpdffname);
  if initflag==1, initflag=0; end;
elseif ~initflag
  fprintf('%s: ERROR: jpdf file %s not found\n',mfilename,jpdffname);
  return;
else
  fprintf('%s: initializing jpdf file %s\n',mfilename,jpdffname);
end;

% compute scalefactor for vol1
[hc,bv] = hist(vol1.imgs(:),1000);
if ~initflag
  sf1 = ComputeHistScalefactor(hc,bv,v1_hc,v1_bv,v1_bn0);
  vol1.imgs = sf1*vol1.imgs;
  fprintf('%s: sf1 = %f\n',mfilename,sf1);
else
  sf1 = 1;
  v1_hc = hc;
  v1_bv = bv;
  v1_bn0 = 50;
end

% compute scalefactor for vol2
[hc,bv] = hist(vol2.imgs(:),1000);
if ~initflag
  sf2 = ComputeHistScalefactor(hc,bv,v2_hc,v2_bv,v2_bn0);
  vol2.imgs = sf2*vol2.imgs;
  fprintf('%s: sf2 = %f\n',mfilename,sf2);
else
  sf2 = 1;
  v2_hc = hc;
  v2_bv = bv;
  v2_bn0 = 50;
end

% make brain mask for v1
if isempty(volmask)
  brainmesh = dctMorph_SkullStrip_amd(vol1, [4 4 4], [5 5 5], false);
  [maskedvol, volmask] = getmaskvol(vol1, brainmesh, eye(4));
  clear maskedvol
  volmask.imgs = double(real(smooth3d(single(volmask.imgs),...
    mask_smoothing(1),mask_smoothing(2),mask_smoothing(3))));
  volmask.imgs = 1.0*(volmask.imgs>0.01);
end;
fprintf('\n');

% register and resample vol2 to vol1
if initflag
  vol2_res = vol_resample_pad(vol2, vol1, Mreg0, 1);
  vol2_res.imgs(find(isnan(vol2_res.imgs))) = 0;
  [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=vols_jhist_mask_amd(vol1, vol2_res, volmask, 1, eye(4), eye(4));
  jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
end
[M_v1_to_v2, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol1, vol2, volmask,...
  Mreg0, 0, 4, scales, [], [], [],...
  jpdf_mask, bins1_mask, bins2_mask, interpm);
vol2_res = vol_resample_pad(vol2, vol1, M_v1_to_v2, 1);
vol2_res.imgs(find(isnan(vol2_res.imgs))) = 0;

if initflag
  % calculate jpdf from registered volumes and reregister
  [sum_log10_val, jpdf_mask, bins1_mask, bins2_mask, jentropy]=vols_jhist_mask_amd(vol1, vol2_res, volmask, 1, eye(4), eye(4));
  jpdf_mask = max(0,smooth2(jpdf_mask,11,11)); % Smooth jpdf function
  Mreg0 = M_v1_to_v2;
  [M_v1_to_v2, min_cost] = rbreg_vol2vol_jpdf_mask_amd(vol1, vol2, volmask,...
    Mreg0, 0, 4, scales, [], [], [],...
    jpdf_mask, bins1_mask, bins2_mask,interpm);
  vol2_res = vol_resample_pad(vol2, vol1, M_v1_to_v2, 1);
  vol2_res.imgs(find(isnan(vol2_res.imgs))) = 0;

  fprintf('%s: saving jpdf file %s\n',mfilename,jpdffname);
  save(jpdffname,'jpdf_mask','bins1_mask','bins2_mask','v2_hc','v2_bv','v2_bn0','v1_hc','v1_bv','v1_bn0');
end

fprintf('%s: min_cost = %f\n',mfilename,min_cost);

