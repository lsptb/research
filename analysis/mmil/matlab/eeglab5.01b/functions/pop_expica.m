% pop_expica() - export ICA weights or inverse matrix
%
% Usage:
%   >> pop_expica( EEG, whichica);             % a window pops up
%   >> pop_expica( EEG, whichica, filename );
%
% Inputs:
%   EEG         - EEGLAB dataset
%   whichica    - ['weights'|'inv'] export ica 'weights' or ica inverse
%                 matrix ('inv'). Note: for 'weights', the function 
%                 export the product of the sphere and weights matrix.
%   filename    - text file name
% 
% Author: Arnaud Delorme, CNL / Salk Institute, Mai 14, 2003
%
% See also: pop_export()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Mai 14, 2003, Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: pop_expica.m,v $
% Revision 1.4  2004/07/29 19:00:17  arno
% saving as double
% .,
%
% Revision 1.3  2003/05/22 15:04:30  arno
% header typo
%
% Revision 1.2  2003/05/14 15:28:19  arno
% sphere matrixc
%
% Revision 1.1  2003/05/14 15:10:20  arno
% Initial revision
%

function com = pop_expica(EEG, whichica, filename); 
    
com = '';
if nargin < 1 
    help pop_expica;
    return;
end;

if nargin < 2
    whichica = 'weights';
end;
switch lower(whichica)
 case {'weights' 'inv'}, ;
 otherwise error('Unrecognized option for ''whichica'' parameter');
end;

if nargin < 3
	% ask user
	[filename, filepath] = uiputfile('*.*', [ 'File name for ' ...
                        fastif(strcmpi(whichica, 'inv'), 'inverse', 'weight') ' matrix -- pop_expica()']); 
    drawnow;
	if filename == 0 return; end;
	filename = [filepath filename];
end;

% save datas
% ----------
if strcmpi(whichica, 'inv')
    tmpmat = double(EEG.icawinv);
else
    tmpmat = double(EEG.icaweights*EEG.icasphere);
end;
save(filename, '-ascii', 'tmpmat');

com = sprintf('EEG = pop_expica(%s, ''%s'', ''%s'');', inputname(1), whichica, filename); 

return;
