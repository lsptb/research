% std_readtopo() - returns the scalp map of a specified ICA component, assumed
%                  to have been saved in a Matlab file, [dataset_name].icatopo, 
%                  in the same directory as the dataset file. If this file does 
%                  not exist, use std_topo() to create it, else a pre-clustering 
%                  function that calls it: pop_preclust() or eeg_preclust().  
% Usage:    
%   >> [grid, y, x ] = std_readtopo(ALLEEG, setindx, component);  
%
% Inputs:
%   ALLEEG     - vector of EEG datasets (can also be one EEG set). 
%                must contain the dataset of interest (see 'setindx' below).
%   setind     - [integer] an index of an EEG dataset in the ALLEEG
%                structure, for which to get the component ERP.
%   component  - [integer] index of the component for which the scalp map 
%                grid should be returned. 
% Outputs:
%   grid      - square scalp-map color-value grid for the requested ICA component 
%               in the specified dataset, an interpolated Cartesian grid as output 
%               by topoplot(). 
%   y         - y-axis values for the interpolated grid
%   x         - x-axis values of the interpolated grid
%
%  See also  std_topo(), std_preclust()
%
% Authors: Arnaud Delorme, Hilit Serby, SCCN, INC, UCSD, February, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, October 11, 2004, hilit@sccn.ucsd.edu
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

% $Log: std_readtopo.m,v $
% Revision 1.8  2006/03/14 02:39:40  scott
% help msg
%
% Revision 1.7  2006/03/11 07:23:51  arno
% header
%
% Revision 1.6  2006/03/10 00:37:17  arno
% error msg
%
% Revision 1.5  2006/03/09 19:00:42  arno
% reading Matlab file
%
% Revision 1.4  2006/03/09 00:00:54  arno
%  now saving Matlab file
%
% Revision 1.3  2006/03/08 20:32:48  arno
% rename func
%
% Revision 1.2  2006/03/07 22:14:40  arno
% use fullfile
%

function [grid, yi, xi ] = std_readtopo(ALLEEG, abset, comp)

grid = [];
yi = [];
xi = [];
filename = fullfile( ALLEEG(abset).filepath,[ ALLEEG(abset).filename(1:end-3) 'icatopo']);
try
    topo = load( '-mat', filename, ...
                 [ 'comp' int2str(comp) '_grid'], ...
                 [ 'comp' int2str(comp) '_x'], ...
                 [ 'comp' int2str(comp) '_y'] );
catch
    error( [ 'Cannot read file ''' filename '''' ]);
end;
    
grid = getfield(topo, [ 'comp' int2str(comp) '_grid']);
yi   = getfield(topo, [ 'comp' int2str(comp) '_y']);
xi   = getfield(topo, [ 'comp' int2str(comp) '_x']);

return;
