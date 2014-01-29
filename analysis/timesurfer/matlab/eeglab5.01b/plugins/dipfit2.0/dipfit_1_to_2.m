% dipfit_1_to_2() - convert dipfit 1 structure to dipfit 2 structure.
%
% Usage:
%  >> EEG.dipfit = dipfit_1_to_2(EEG.dipfit);
%
% Note:
% For non-standard BESA models (where the radii or the conductances
% have been modified, users must create a new model in Dipfit2 from
% the default BESA model. 
%
% Author: Arnaud Delorme, SCCN, La Jolla 2005

% Copyright (C) Arnaud Delorme, SCCN, La Jolla 2005
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

% $Log: dipfit_1_to_2.m,v $
% Revision 1.3  2005/03/11 16:05:52  arno
% implement .vol
%
% Revision 1.2  2005/03/10 22:11:18  arno
% *** empty log message ***
%

function newdipfit = dipfit_1_to_2( dipfit );

    if isfield( dipfit, 'model')
        newdipfit.model = dipfit.model;
    end;
    if isfield( dipfit, 'chansel')
        newdipfit.chansel = dipfit.chansel;
    end;
    
    ind = 1; % use first template (BESA)
    newdipfit.coordformat = template_models{ind}{2};
    newdipfit.mrifile     = template_models{ind}{3};
    newdipfit.chanfile'   = template_models{ind}{4};
    
    if ~isfield(dipfit, 'vol')
        newdipfit.hdmfile = template_models{ind}{1};
    else
        newdipfit.vol     = dipfit.vol;

        %if length(dipfit.vol) == 4
            %if ~all(dipfit.vol == [85-6-7-1 85-6-7 85-6 85]) | ...
            %   ~all(dipfit.c   == [0.33 1.00 0.0042 0.33]) | ...
            %   ~all(dipfit.o = [0 0 0])
            %    disp('Warning: Conversion from dipfit 1 to dipfit 2 can only deal');
            %    disp('         with standard (not modified) BESA model');
            %    disp('         See "help dipfit_1_to_2" to convert this model');
            %    newdipfit = [];
            %end;
        %end;
    end;
                
                