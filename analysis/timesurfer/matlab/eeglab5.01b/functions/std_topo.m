% std_topo() - uses topoplot() to get the interpolated Cartesian grid of the 
%               specified component topo maps. The topo map grids are saved
%               into a (.icatopo) file and a pointer to the file is stored 
%               in the EEG structure. If such a file already exists, 
%               loads the information from it. 
%
%               Returns the topo map grids of all the requested components. Also
%               returns the EEG sub-structure etc (i.e EEG.etc), which is modified 
%               with a pointer to the float file and some information about the file. 
% Usage:
%               >> X = std_topo(EEG, components, option);  
%
%                  % Returns the ICA topo map grid for a dataset. 
%                  % Updates the EEG structure in the Matlab environment and re-saves
% Inputs:
%   EEG        - an EEG dataset structure. 
%   components - [numeric vector] components in the EEG structure to compute topo maps
%                      {default|[] -> all}      
%   option     - ['gradient'|'laplacian'|'none'] compute gradient or laplacian of
%                the scale topography. {default is 'none' = the interpolated topo map}
% Outputs:
%   X          - the topo map grid of the requested ICA components, each grid is 
%                     one ROW of X. 
%
% File output: [dataset_name].icatopo
%  
% Authors:  Hilit Serby, Arnaud Delorme, SCCN, INC, UCSD, January, 2005
%
% See also  topoplot(), std_erp(), std_ersp(), std_spec(), std_preclust()

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

% $Log: std_topo.m,v $
% Revision 1.15  2006/03/10 00:20:39  arno
% same
%
% Revision 1.14  2006/03/10 00:19:49  arno
% removing reference to etc field
%
% Revision 1.13  2006/03/09 19:27:39  arno
% fix function crash
%
% Revision 1.12  2006/03/09 19:21:29  arno
% header
%
% Revision 1.11  2006/03/09 19:21:08  arno
% header
%
% Revision 1.10  2006/03/09 19:20:45  arno
% saving topographies
%
% Revision 1.9  2006/03/09 00:00:52  arno
% now saving Matlab file
%
% Revision 1.8  2006/03/08 20:27:38  arno
% rename func
%
% Revision 1.7  2006/03/07 03:46:43  scott
% reworked help msg; made function accept specified component list or [] -sm
%
% Revision 1.6  2006/03/06 23:45:28  arno
% adding gradient and laplacian
%
% Revision 1.5  2006/03/06 23:17:34  arno
% fix resave
%
% Revision 1.4  2006/03/03 21:22:07  arno
% remve change folder
%
% Revision 1.3  2006/03/03 21:19:16  arno
% change topo map file name
%
% Revision 1.2  2006/03/03 00:41:38  arno
% now correctly saving data
%

function [X] = std_topo(EEG, comps, option)

if nargin < 1
    help std_topo;
    return;
end;
if isfield(EEG,'icaweights')
   numc = size(EEG.icaweights,1);
else
   error('EEG.icaweights not found');
end
if nargin < 2
   comps = 1:numc;
elseif isempty(comps)
   comps = 1:numc;
end

if nargin < 3
    option = 'none';
end;

% figure; toporeplot(grid,'style', 'both','plotrad', 0.5, 'intrad', 0.5, 'xsurface' ,Xi, 'ysurface',Yi );

% Topo information found in dataset
% ---------------------------------
if exist(fullfile(EEG.filepath, [ EEG.filename(1:end-3) 'icatopo' ]))
    for k = 1:length(comps)
        tmp = std_readtopo( EEG, 1, comps(k));
        if strcmpi(option, 'gradient')
            [tmpx, tmpy]  = gradient(tmp); %Gradient
            tmp = [tmpx(:); tmpy(:)]';
        elseif strcmpi(option, 'laplacian')
            tmp = del2(tmp); %Laplacian
            tmp = tmp(:)';
        else
            tmp = tmp(:)';
        end;
        
        tmp = tmp(find(~isnan(tmp)));
        if k == 1
            X = zeros(length(comps),length(tmp)) ;
        end
        X(k,:) =  tmp;
    end
    return
end
 
all_topos = [];
for k = 1:numc

    % compute topo map grid (topoimage)
    % ---------------------------------
    [hfig grid plotrad Xi Yi] = topoplot( EEG.icawinv(:,k), EEG.chanlocs, ...
                                          'verbose', 'off',...
                                           'electrodes', 'on' ,'style','both',...
                                           'plotrad',0.5,'intrad',0.5,...
                                           'noplot', 'on', 'chaninfo', EEG.chaninfo);
    all_topos = setfield(all_topos, [ 'comp' int2str(k) '_grid' ], grid);
    all_topos = setfield(all_topos, [ 'comp' int2str(k) '_x' ]   , Xi(:,1));
    all_topos = setfield(all_topos, [ 'comp' int2str(k) '_y' ]   , Yi(:,1));
    
end

% Save topos in file
% ------------------
all_topos.datatype = 'TOPO';
tmpfile = fullfile( EEG.filepath, [ EEG.filename(1:end-3) 'icatopo' ]); 
std_savedat(tmpfile, all_topos);

for k = 1:length(comps)
    tmp =  getfield(all_topos, [ 'comp' int2str(comps(k)) '_grid' ]);
    
    if strcmpi(option, 'gradient')
        [tmpx, tmpy]  = gradient(tmp); % Gradient
        tmp = [tmpx(:); tmpy(:)]';
    elseif strcmpi(option, 'laplacian')
        tmp = del2(tmp); % Laplacian
        tmp = tmp(:)';
    else
        tmp = tmp(:)';
    end;

    tmp = tmp(find(~isnan(tmp)));
    if k == 1
        X = zeros(length(comps),length(tmp)) ;
    end
    X(k,:) =  tmp;
end
