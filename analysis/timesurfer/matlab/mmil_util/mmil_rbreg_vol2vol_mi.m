function [M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_mi(vol1,vol2,...
  volmask,maskflag,mask_smoothing,scales)
%function [M_v1_to_v2,volmask]=mmil_rbreg_vol2vol_mi(vol1,vol2,...
%  [volmask],[maskflag],[mask_smoothing],[scales])
%
% Required Input:
%   vol1: 3D volume of reference scan (e.g. structural scan)
%   vol2: 3D volume of scan to be registered (e.g. first frame of functional scan)
%
% Optional Input:
%   volmask: mask for vol1 (e.g. brain mask)
%     {default = []}
%   maskflag: [0|1] toggle whether to use brain mask for vol1
%     if maskflag=1 and volmask=[], will
%       automatically generate one by registering to atlas
%       (will only work if vol1 is T1-weighted)
%     {default = 1}
%   mask_smoothing: vector of x,y,z smoothing sigma (voxels) for mask
%     {default = [3,3,3]}
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
% Created : 05/05/08 by Don Hagler
% Last Mod: 05/05/08 by Don Hagler
% Last Mod: 07/02/08 by Don Hagler
%

if (~mmil_check_nargs(nargin,2)) return; end;

M_v1_to_v2 = [];

if ~exist('volmask','var'), volmask = []; end;
if ~exist('maskflag','var') || isempty(maskflag), maskflag = 1; end;
if ~exist('mask_smoothing','var') || isempty(mask_smoothing)
  mask_smoothing = [3 3 3];
end;
if ~exist('scales','var') || isempty(scales)
  scales=[0 83 49 27 16 9 5 3 2 1];
end;
if (scales(1)~=0) scales=[0,scales]; end;

if maskflag
  % make brain mask for v1
  if isempty(volmask)
    brainmesh = dctMorph_SkullStrip_amd(vol1, [4 4 4], [5 5 5], true);
    [maskedvol, volmask] = getmaskvol(vol1, brainmesh, eye(4));
    clear maskedvol
    volmask.imgs = double(real(smooth3d(single(volmask.imgs),...
      mask_smoothing(1),mask_smoothing(2),mask_smoothing(3))));
    volmask.imgs = 1.0*(volmask.imgs>0.01);
    fprintf('\n');
  end;
  % register vol2 to vol1
  mstep = 1;
  [M_v1_to_v2,min_cost] = rbreg_vol2vol_mi_mask_djh(vol1,vol2,volmask,true,...
    mstep,scales);
else
  % register vol2 to vol1
  range=[-70 -70 -70; 70 70 70; 2. 2. 2.];
  [M_v1_to_v2,min_cost] = rbreg_vol2vol_mi(vol1,vol2,true,...
    range,scales);
end;

fprintf('%s: min_cost = %f\n',mfilename,min_cost);

