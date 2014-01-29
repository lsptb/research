function [smoothedvals] = fs_smooth(surf,vals,niter)
%function [smoothedvals] = fs_smooth(surf,vals,niter)
%
% surf is a structure containing:
%   coords: x,y,z coords for each surface vertex
%   nbrs:   vertex numbers of neighbors for each vertex
%   nverts: number of vertices
%
% vals is vector with nverts members
%
% niter is number of iterations (smoothing steps)
%
% Created: <06/11/06 Don Hagler
% Last Mod: 03/13/08  by Don Hagler
%

smoothedvals = [];

if nargin<3
  help(mfilename);
  return;
end;

if isempty(vals)
  error('vals is empty');
end;

if size(vals,2)~=1
  error('vals must have a single column (it has %d)',size(vals,2));
end;

if size(vals,1)~=surf.nverts
  error('number of vals (%d) does not match number of verts (%d)',...
    size(vals,1),surf.nverts);
end;

%fprintf('%s(%d): smoothing',mfilename,niter);
vals = [0;vals]; % create dummy vertex with index 1, value 0
maxnbrs = size(surf.nbrs,2);
nbrs = [ones(maxnbrs,1)';surf.nbrs + 1]; % nbrs contains zeros, change to 1
num_nbrs = sum(nbrs>1,2)+1; % including vertex and its neighbors
for iter=1:niter
%  fprintf('.');
  Y=[vals, vals(nbrs(:,:))]; % sum values from vertex and its neighbors
  vals=sum(Y,2)./num_nbrs;
end;
smoothedvals = vals(2:end);

%fprintf('\n');
