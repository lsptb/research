function [dist] = warpfun(M, input, target, varargin);

% WARPFUN is deprecated, please use WARP_ERROR

% warning('this function is deprecated, please use WARP_ERROR');
[dist] = warp_error(M, input, target, varargin{:});

