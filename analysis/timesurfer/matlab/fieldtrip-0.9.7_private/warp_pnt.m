function [result, M] = warp_pnt(input, target, method);

% WARP_PNT is deprecated, please use WARP_OPTIM

% warning('this function is deprecated, please use WARP_OPTIM');
[result, M] = warp_optim(input, target, method);


