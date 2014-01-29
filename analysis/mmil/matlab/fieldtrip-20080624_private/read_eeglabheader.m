% read_eeglabheader() - import EEGLAB dataset files
%
% Usage:
%   >> header = read_eeglabheader(filename);
%
% Inputs:
%   filename - [string] file name
%
% Outputs:
%   header   - FILEIO toolbox type structure
%
% Author: Arnaud Delorme, SCCN, INC, UCSD, 2008-

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2008 Arnaud Delorme, SCCN, INC, UCSD, arno@salk.edu
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

% $Log: read_eeglabheader.m,v $
% Revision 1.2  2008/04/21 18:45:59  roboos
% fixed bug, ori should be orig
%
% Revision 1.1  2008/04/18 14:04:48  roboos
% new implementation by Arno, shoudl be tested
%

function header = read_eeglabheader(filename)

if nargin < 1
  help read_eeglabheader;
  return;
end;

if ~isstruct(filename)
  load('-mat', filename);
else
  EEG = filename;
end;

header.Fs          = EEG.srate;
header.nChans      = EEG.nbchan;
header.nSamples    = EEG.pnts;
header.nSamplesPre = -EEG.xmin*EEG.srate;
header.nTrials     = EEG.trials;
header.label       = { EEG.chanlocs.labels }';
for ind = 1:length( EEG.chanlocs )
  header.elec.label{ind} = EEG.chanlocs(ind).labels;
  if ~isempty(EEG.chanlocs(ind).X)
    % this channel has a position
    header.elec.pnt(ind,1) = EEG.chanlocs(ind).X;
    header.elec.pnt(ind,2) = EEG.chanlocs(ind).Y;
    header.elec.pnt(ind,3) = EEG.chanlocs(ind).Z;
  end;
end;

% remove data
% -----------
%if isfield(EEG, 'datfile')
%    if ~isempty(EEG.datfile)
%        EEG.data = EEG.datfile;
%    end;
%else
%    EEG.data = 'in set file';
%end;
EEG.icaact = [];

header.orig = EEG;
