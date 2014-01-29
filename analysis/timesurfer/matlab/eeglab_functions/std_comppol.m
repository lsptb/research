% std_comppol() - inverse component polarity in a component cluster
%
% Usage: [compout pol] = std_comppol(compin);
%
% Inputs:
%    compin  - component scalp maps, one per column.
%
% Outputs:
%    compout - component scalp maps some of them with inverted
%              polarities, one per column.
%    pol     - logical vector of component with inverted 
%              polarities (same length as the number of rows in 
%              compin)
%
% Author: Arnaud Delorme & Hilit Serby, SCCN, INC, UCSD, 2004

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2004 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: std_comppol.m,v $
% Revision 1.2  2006/03/12 04:22:29  arno
% comppol -> std_comppol
%
% Revision 1.1  2004/08/26 23:06:27  arno
% Initial revision
%

function [compin, pol] = std_comppol(compin);

if nargin < 1
    help std_comppol;
    return;
end;

% run several iterations
% ----------------------
compini = compin;
count = 1;
while count < size(compin,2)
    count = count+1;
    
    % remove diagonal and put 0 and 1
    % -------------------------------
    r = corrcoef( compin );
    r = r - eye(size(r));
    
    % invert component polarities
    % ---------------------------
    for index = 1:size(compin,2)
        if sign(sum(r(:,index))) < 0
            r(:,index) = -r(:,index);
            r(index,:) = -r(index,:);
            r(index,index) = r(index,index);
            compin(:,index) = -compin(:,index);
        end;
    end;
end;

% get polarities
% --------------
for index = 1:size(compin,2)
    if compin(1,index) == compini(1,index), pol(index) = 1;
    else                                    pol(index) = -1;
    end;
end;

    