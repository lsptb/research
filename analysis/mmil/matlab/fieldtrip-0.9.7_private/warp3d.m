function [warped] = warp3d(M, input, method);

% WARP3D is deprecated, please use WARP_APPLY

% warning('this function is deprecated, please use WARP_APPLY');
[warped] = warp_apply(M, input, method);

