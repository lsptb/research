function [H] = traditional(f);

% TRADITIONAL creates the homogenous spatial transformation matrix
%   for a 9 parameter traditional "Talairach-model" transformation
%
% H = traditional(f)
%
% The transformation vector f should contain the 
%	x-shift
%	y-shift
%	z-shift
% followed by the
%	pitch (rotation around x-axis)
%	roll  (rotation around y-axis)
%	yaw   (rotation around z-axis)
% followed by the 
%	x-rescaling factor
%	y-rescaling factor
%	z-rescaling factor
%
% See also WARP3D

% Copyright (C) 2000-2004, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% $Log: traditional.m,v $
% Revision 1.2  2004/05/19 09:57:07  roberto
% added GPL copyright statement, added CVS log item
%

% compute the homogenous transformation matrix for the translation
T = eye(4,4);
if isa( f, 'sym')
    T = sym(T); 
end;
T(1,4) = f(1);
T(2,4) = f(2);
T(3,4) = f(3);

% precompute the sin/cos values of the angles
cX = cos(f(4));
cY = cos(f(5));
cZ = cos(f(6));
sX = sin(f(4));
sY = sin(f(5));
sZ = sin(f(6));

% compute the homogenous transformation matrix for the rotation
R = eye(4,4);
if isa( f, 'sym')
    R = sym(R); 
end;
R(1,1) = cZ*cY + sZ*sX*sY;
R(1,2) = sZ*cY + cZ*sX*sY;
R(1,3) =            cX*sY;
R(2,1) = -sZ*cX;
R(2,2) =  cZ*cX;
R(2,3) =     sX;
R(3,1) =  sZ*sX*cY - cZ*sY;
R(3,2) = -cZ*sX*cY - sZ*sY;
R(3,3) =             cX*cY;

% compute the homogenous transformation matrix for the scaling
S = eye(4,4);
if isa( f, 'sym')
    S = sym(S); 
end;
S(1,1) = f(7);
S(2,2) = f(8);
S(3,3) = f(9);

% compute the homogenous coordinate transformation matrix for use by WARP3D
H = T*R*S;

