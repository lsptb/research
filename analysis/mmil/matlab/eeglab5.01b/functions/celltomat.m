% celltomat() - convert cell array to matrix
%
% Usage: >> M = celltomat( C );
%
% Author: Arnaud Delorme, CNL / Salk Institute, Jan 25 2002
%
% Note: This function overloads the neuralnet toolbox function CELLTOMAT,
% but does not have all its capacities. Delete this version if you have 
% the toolbox.

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Jan 25 2002 Arnaud Delorme, CNL / Salk Institute  
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

% $Log: celltomat.m,v $
% Revision 1.3  2005/05/24 17:52:29  arno
% cell2mat -> celltomat
%
% Revision 1.2  2003/01/01 17:55:53  scott
% header edits
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%

function M = celltomat( C, varargin );

if nargin < 1
	help celltomat;
	return;
end;

for i=1:size(C,1)
    for j=1:size(C,2)
        M(i,j) = C{i,j};
    end;
end;   
