% pop_averef() - Convert an EEG dataset to average reference.
%                This function is obsolete. See pop_reref() instead.
%
%
% Usage:
%       >> EEGOUT = pop_averef( EEG, confirm);
%
% Inputs:
%   EEG         - input dataset
%   confirm     - [0|1] ask for confirmation
%
% Inputs:
%   EEGOUT      - output dataset (with variable EEG.rmave added)
%
% Author: Arnaud Delorme, CNL / Salk Institute, 22 March 2002
%
% See also: eeglab(), averef()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 22 March 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: pop_averef.m,v $
% Revision 1.13  2003/12/10 02:45:14  arno
% same
%
% Revision 1.12  2003/12/10 02:44:03  arno
% msg
%
% Revision 1.11  2002/11/12 19:02:37  arno
% debugging command line call
%
% Revision 1.10  2002/09/05 00:31:12  scott
% added EEG.rmave to output dataset - to allow re-referencing -sm
%
% Revision 1.9  2002/08/26 22:04:08  arno
% debug
%
% Revision 1.8  2002/08/21 02:14:35  arno
% debug
%
% Revision 1.7  2002/08/21 02:13:39  arno
% more messages
%
% Revision 1.6  2002/08/19 21:56:28  arno
% debug for MAC
%
% Revision 1.5  2002/08/12 18:25:49  arno
% questdlg2
%
% Revision 1.4  2002/08/12 16:22:53  arno
% questdlg2
%
% Revision 1.3  2002/04/18 16:10:10  scott
% changed 'yes' to 'Yes' for uniformity with other flags -sm
%
% Revision 1.2  2002/04/11 17:58:20  arno
% computing average reference of components
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

function [EEG, com] = pop_averef( EEG, confirm);

com = '';
if nargin < 1
   help pop_averef;
   return;
end;   
if isempty(EEG.data)
    error('Pop_averef: cannot process empty data');
end;

if nargin < 2 | confirm == 1
    % which set to save
	% -----------------
	 ButtonName=questdlg2( strvcat('Convert the data to average reference?', ...
								   'Note: ICA activations will also be converted if they exist...'), ...
	        'Average reference confirmation -- pop_averef()', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', return;
	 end;
	 confirm = 0;
end;

EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
if ~isempty(EEG.icaweights)
	disp('pop_averef(): converting ICA weight matrix to average reference (see >> help averef)');
	[EEG.data EEG.icaweights EEG.icasphere EEG.rmave] = averef(EEG.data,EEG.icaweights,EEG.icasphere);
	EEG.icawinv = [];
	if size(EEG.icaweights,1) > EEG.nbchan
		disp('Warning: one or more channels may have been removed; component weight re-referencing may be inaccurate'); 
	end;
	if size(EEG.icasphere,1) <  EEG.nbchan
		disp('Warning: one or more components may have been removed; component weight re-referencing could be inaccurate'); 
	end;
else
	EEG.data = averef(EEG.data);
end;	
EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
EEG.averef = 'Yes';
EEG.icaact = [];
EEG = eeg_checkset(EEG);

com = sprintf('%s = pop_averef( %s, %d);', inputname(1), inputname(1), confirm);
return;
