% eeg_chantype() - Returns the channel indices of the desired channel type(s).
%
% Usage:
%   >> indices = eeg_chantype(struct, types )
%
% Inputs:
%   struct     - EEG.chanlocs data structure returned by readlocs() containing
%                channel location, type and gain information.
%
% Optional input
%   types      - [cell array] cell array containing types ...
%
% Output:
%   indices    -
%
% Author: Toby Fernsler, Arnaud Delorme, Scott Makeig
%
% See also: topoplot()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2005 Toby Fernsler, SCCN, INC, UCSD, toby@sccn.ucsd.edu
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

% $Log: eeg_chantype.m,v $
% Revision 1.5  2006/04/12 02:37:19  toby
% Dealing with no chantype match.
%
% Revision 1.4  2005/11/30 19:38:18  arno
% same
%
% Revision 1.3  2005/11/30 19:36:52  arno
% handling empty channel types
%
% Revision 1.2  2005/10/27 22:22:44  arno
% header
%

function indices = eeg_chantype(data,chantype)

    if nargin < 1
        help eeg_chantype;
    end;
   
    if ischar(chantype), chantype = cellstr(chantype); end
    if ~iscell(chantype), 
        error( 'chantype must be cell array, e.g. {''EEG'', ''EOG''}, or single character string, e.g.''EEG''.'); 
    end
    
    % Define 'datatype' variable, listing the type of each channel.
    % ------------------------------------------------------------
    if isfield(data,'type')
        datatype = {data.type};
    elseif isfield(data,'chanlocs') && isfield(data.chanlocs,'type')
        datatype = {data.chanlocs.type};
    else error('Incorrect ''data'' input. Should be ''EEG'' or ''loc_file'' structure variable in the format associated with EEGLAB.');
    end
    
    % seach for types
    % ---------------
    k = 1;
    plotchans = [];
    for i = 1:length(chantype)
        for j = 1:length(datatype)
            if strcmpi(chantype{i},char(datatype{j}))
                plotchans(k) = j;
                k = k + 1;
            end;
        end
    end
    indices = sort(plotchans);
