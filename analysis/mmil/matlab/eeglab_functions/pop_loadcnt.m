% pop_loadcnt() - load a neuroscan CNT file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_loadcnt; % pop-up window mode
%   >> EEG = pop_loadcnt( filename, 'key', 'val', ...);
%
% Graphic interface:
%   "Data fomat" - [checkbox] 16-bits or 32-bits. We couldn't find in the
%                   data file where this information was stored. Command
%                   line equivalent in loadcnt() 'dataformat'.
%   "Time interval in seconds" - [edit box] specify time interval [min max]
%                   to import portion of data. Command line equivalent
%                   in loadcnt: 't1' and 'lddur'
%   "Import keystrokes" - [checkbox] set this option to import keystroke
%                   event types in dataset. Command line equivalent
%                   'keystroke'.
%   "loadcnt() 'key', 'val' params" - [edit box] Enter optional loadcnt()
%                   parameters.
%
% Inputs:
%   filename       - file name
%
% Optional inputs:
%   'keystroke'    - ['on'|'off'] set the option to 'on' to import 
%                    keystroke event types. Default is off.
%   Same as loadcnt() function.
% 
% Outputs:
%   EEG            - EEGLAB data structure
%
% Note: 
% 1) This function extract all non-null event from the CNT data structure.
% Null events are usually associated with internal signals (recalibrations...).
% 2) The "Average reference" edit box had been remove since the re-referencing
% menu of EEGLAB offers more options to re-reference data.
% 3) The 'blockread' has been disabled since we found where this information
% was stored in the file.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: loadcnt(), eeglab()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2001 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: pop_loadcnt.m,v $
% Revision 1.30  2005/10/24 23:33:15  arno
% fix filename
%
% Revision 1.29  2005/05/24 17:28:06  arno
% remove cell2mat
%
% Revision 1.28  2005/05/12 15:51:23  arno
% deal with keypad etc...
%
% Revision 1.27  2004/12/08 18:01:55  arno
% debug offset
%
% Revision 1.26  2004/11/17 00:13:49  arno
% adding 1 to offset
%
% Revision 1.25  2004/09/14 23:31:11  arno
% nothing
%
% Revision 1.24  2004/04/29 17:13:17  arno
% ignoring keystroke events
%
% Revision 1.23  2004/03/17 22:52:29  arno
% import keystroke event type
%
% Revision 1.22  2003/12/17 23:20:49  arno
% debug last
%
% Revision 1.21  2003/12/17 23:19:51  arno
% importing channel labels
%
% Revision 1.20  2003/10/30 19:38:39  arno
% removing warnings
%
% Revision 1.19  2003/08/21 17:38:40  arno
% nothing
%
% Revision 1.18  2003/08/01 21:16:30  arno
% header info
%
% Revision 1.17  2003/07/29 21:26:34  arno
% new loadcnt version
%
% Revision 1.16  2003/07/28 22:07:44  arno
% no more block factor
%
% Revision 1.15  2003/07/25 00:59:26  arno
% removing blockread option
%
% Revision 1.14  2003/06/19 16:14:59  arno
% makeur
%
% Revision 1.13  2003/05/20 01:49:01  arno
% allowing to import keystrokes
%
% Revision 1.12  2003/05/20 00:46:11  arno
% debug if no events
%
% Revision 1.11  2003/05/14 17:16:27  arno
% putting time range in gui
%
% Revision 1.10  2003/04/23 21:29:49  arno
% removing filepath
%
% Revision 1.9  2003/04/10 17:56:52  arno
% debuging function and history
%
% Revision 1.8  2003/04/10 17:35:29  arno
% header and history
%
% Revision 1.7  2003/03/05 19:48:42  arno
% removing matlab warning
%
% Revision 1.6  2002/11/23 21:10:02  arno
% importing type of event
%
% Revision 1.5  2002/10/22 23:57:21  arno
% change default blockread
%
% Revision 1.4  2002/10/15 17:01:13  arno
% drawnow
%
% Revision 1.3  2002/08/12 02:40:59  arno
% inputdlg2
%
% Revision 1.2  2002/08/06 21:33:30  arno
% spelling
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, command] = pop_loadcnt(filename, varargin); 
command = '';
EEG = [];

if nargin < 1 

	% ask user
	[filename, filepath] = uigetfile('*.CNT;*.cnt', 'Choose a CNT file -- pop_loadcnt()'); 
    drawnow;
	if filename == 0 return; end;

	% popup window parameters
	% -----------------------
    callback16 = 'set(findobj(gcbf, ''tag'', ''32''), ''value'', ~get(gcbo, ''value''));';
    callback32 = 'set(findobj(gcbf, ''tag'', ''16''), ''value'', ~get(gcbo, ''value''));';
    uigeom       = { [1.3 0.5 0.5] [1 0.5] [1.09 0.13 0.4] [1 0.5] };
    uilist       = { { 'style' 'text' 'string' 'Data format' } ...
                     { 'style' 'checkbox' 'tag' '16' 'string' '16-bits' 'value' 1 'callback' callback16 } ...
                     { 'style' 'checkbox' 'tag' '32' 'string' '32-bits' 'value' 0 'callback' callback32 } ...
                     { 'style' 'text' 'string' 'Time interval in seconds (i.e. [0 100]; default all):' } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' 'Check to Import keystrokes:' } ...
                     { 'style' 'checkbox' 'string' '' } { } ...
                     { 'style' 'text' 'string' 'loadcnt() ''key'', ''val'' params' } ...
                     { 'style' 'edit' 'string' '' } };
	result = inputgui( uigeom, uilist, 'pophelp(''pop_loadcnt'')', 'Load a CNT dataset');    
	if length( result ) == 0 return; end;

	% decode parameters
	% -----------------
    options = [];
    if result{1}, options = [ options ', ''dataformat'', ''int16''' ];
    else          options = [ options ', ''dataformat'', ''int32''' ];
    end;
    if ~isempty(result{3}), 
        timer =  eval( [ '[' result{3} ']' ]);
        options = [ options ', ''t1'', ' num2str(timer(1)) ', ''lddur'', '  num2str(timer(2)-timer(1)) ]; 
    end;   
    if result{4}, options = [ options ', ''keystroke'', ''on''' ]; end;
    if ~isempty(result{5}), options = [ options ',' result{5} ]; end;
else
	options = vararg2str(varargin);
end;

% load datas
% ----------
EEG = eeg_emptyset;
if exist('filepath')
	fullFileName = sprintf('%s%s', filepath, filename);
else
	fullFileName = filename;
end;	
if nargin > 0
	if ~isempty(varargin)
		r = loadcnt( fullFileName, varargin{:});
	else
		r = loadcnt( fullFileName);
	end;	
else
	eval( [ 'r = loadcnt( fullFileName ' options ');' ]);
end;

if isfield(r, 'dat')
    error('pop_loadcnt is not compatible with current loadcnt version, please use latest loadcnt() version');
end;
EEG.data            = r.data;
EEG.comments        = [ 'Original file: ' fullFileName ];
EEG.setname 		= 'CNT file';
EEG.nbchan          = r.header.nchannels; 

% inport events
% -------------
I = 1:length(r.event);
if ~isempty(I)
    EEG.event(1:length(I),1) = [ r.event(I).stimtype ];
    EEG.event(1:length(I),2) = [ r.event(I).offset ]+1;
    EEG.event = eeg_eventformat (EEG.event, 'struct', { 'type' 'latency' });
end;

% modified by Andreas Widmann  2005/05/12  14:15:00
temp = find([r.event.accept_ev1] == 14 | [r.event.accept_ev1] == 11); % 14: Discontinuity, 11: DC reset
if ~isempty(temp)
    disp('pop_loadcnt note: event field ''type'' set to ''boundary'' for data discontinuities');
    for index = 1:length(temp)
        EEG.event(temp(index)).type = 'boundary';
    end;
end
% end modification

% process keyboard entries
% ------------------------
if ~isempty(findstr('keystroke', lower(options)))
    tmpkbd  = [ r.event(I).keyboard ];
    tmpkbd2 = [ r.event(I).keypad_accept ];
    for index = 1:length(EEG.event)
        if EEG.event(index).type == 0
            if r.event(index).keypad_accept,
                EEG.event(index).type = [ 'keypad' num2str(r.event(index).keypad_accept) ];
            else
                EEG.event(index).type = [ 'keyboard' num2str(r.event(index).keyboard) ];
            end;
        end;
    end;
else
    % removeing keystroke events
    % --------------------------
    rmind = [];
    for index = 1:length(EEG.event)
        if EEG.event(index).type == 0
            rmind = [rmind index];
        end;
    end;
    if ~isempty(rmind)
        fprintf('Ignoring %d keystroke events\n', length(rmind));
        EEG.event(rmind) = [];
    end;
end;

% import channel locations (Neuroscan coordinates are not wrong)
% ------------------------
%x            = celltomat( { r.electloc.x_coord } );
%y            = celltomat( { r.electloc.y_coord } );
for index = 1:length(r.electloc)
    names{index} = deblank(char(r.electloc(index).lab'));
    if size(names{index},1) > size(names{index},2), names{index} = names{index}'; end;
end;
EEG.chanlocs  = struct('labels', names);
%EEG.chanlocs = readneurolocs( { names x y } );
%disp('WARNING: Electrode locations imported from CNT files may not reflect true locations');

EEG.srate    = r.header.rate;
EEG.nbchan   = size(EEG.data,1);
EEG.trials   = 1;
EEG.pnts     = size(EEG.data,2);
EEG          = eeg_checkset(EEG, 'eventconsistency');
EEG          = eeg_checkset(EEG, 'makeur');

if length(options) > 2
    command = sprintf('EEG = pop_loadcnt(''%s'' %s);',fullFileName, options); 
else
    command = sprintf('EEG = pop_loadcnt(''%s'');',fullFileName); 
end;
return;
