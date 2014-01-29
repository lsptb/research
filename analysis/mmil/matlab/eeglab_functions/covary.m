%  covary() - For vectors, covary(X) returns the variance of X.
%             For matrices, covary(X)is a row vector containing the
%             variance of each column of X.
%
% Notes:          
%   covary(X) normalizes by N-1 where N is the sequence length.  
%   This makes covary(X) the best unbiased estimate of the 
%   covariance if X are from a normal distribution.
%   Does not require the Matlab Signal Processing Toolbox
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 2000 

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2000 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: covary.m,v $
% Revision 1.2  2002/11/15 01:50:16  arno
% header for web
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function covout = covary(data)

data = data - mean(mean(data));
if size(data,1) == 1
    data = data';   % make column vector
end
covout = sum(data.*data)/(size(data,1)-1);


