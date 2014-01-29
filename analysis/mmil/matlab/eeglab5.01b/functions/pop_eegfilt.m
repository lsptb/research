% pop_eegfilt() - interactively filter EEG dataset data using eegfilt()
%
% Usage:
%   >> EEGOUT = pop_eegfilt( EEG, locutoff, hicutoff, filtorder);
%
% Graphical interface:
%   "Lower edge ..." - [edit box] Lower edge of the frequency pass band (Hz) 
%                 Same as the 'locutoff' command line input.
%   "Higher edge ..." - [edit box] Higher edge of the frequency pass band (Hz) 
%                 Same as the 'hicutoff' command line input.
%   "Notch filter" - [edit box] provide the notch range, i.e. [45 55] 
%                 for 50 Hz). This option overwrites the low and high edge limits
%                 given above. Set the 'locutoff' and 'hicutoff' values to the
%                 values entered as parameters, and set 'revfilt to 1, to swap
%                 from bandpass to notch filtering.
%   "Filter length" - [edit box] Filter lenghth in point (default: see 
%                 >> help eegfilt). Same as 'filtorder' optional input.
%
% Inputs:
%   EEG       - input dataset
%   locutoff  - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz) {0 -> highpass}
%   filtorder - length of the filter in points {default 3*fix(srate/locutoff)}
%   revfilt   - [0|1] Reverse filter polarity (from bandpass to notch filter). 
%                     Default is 0 (bandpass).
% Outputs:
%   EEGOUT   - output dataset
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegfilt(), eegfiltfft(), eeglab()

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

% $Log: pop_eegfilt.m,v $
% Revision 1.31  2006/03/10 19:47:45  arno
% remove option
%
% Revision 1.30  2005/10/11 16:20:45  arno
% warning for bandpass
%
% Revision 1.29  2005/09/10 00:44:35  arno
% adding multiple dataset capabilities
%
% Revision 1.28  2005/06/29 22:36:45  arno
% typo
%
% Revision 1.27  2005/05/24 17:31:18  arno
% remove cell2mat
%
% Revision 1.26  2005/05/13 15:13:19  arno
% floor instead of round
%
% Revision 1.25  2004/12/10 22:24:01  arno
% erasing ICA matrix
%
% Revision 1.24  2004/08/20 23:32:01  arno
% debug filtering between boundaries
%
% Revision 1.23  2003/12/03 18:31:32  arno
% implementing eegfiltfft
%
% Revision 1.22  2003/09/01 18:14:17  arno
% fixing nargin problem -thanks Petr Janata
%
% Revision 1.21  2003/08/06 00:22:05  arno
% removing debug message
%
% Revision 1.20  2003/08/02 21:29:50  arno
% text
%
% Revision 1.19  2003/07/28 17:39:07  arno
% nothing
%
% Revision 1.18  2003/07/22 17:36:44  arno
% subtract DC for small portion of data
%
% Revision 1.17  2003/07/20 19:19:10  scott
% clarify processing message
%
% Revision 1.16  2003/04/24 22:13:51  arno
% typo
%
% Revision 1.15  2003/04/24 22:13:17  arno
% removing 0 warning
%
% Revision 1.14  2003/04/24 22:08:55  arno
% updating error message
%
% Revision 1.13  2003/04/24 22:06:31  arno
% show message when processing data portion
%
% Revision 1.12  2003/04/24 22:00:16  arno
% debuging boundaries
%
% Revision 1.11  2003/02/17 02:43:43  arno
% reformating text for new functionality in help2html
%
% Revision 1.10  2003/02/16 23:10:43  arno
% adding GUI info
%
% Revision 1.9  2003/01/24 04:03:37  scott
% edits msgs -sm
%
% Revision 1.8  2003/01/24 00:23:35  arno
% debugged revfilt parameter
%
% Revision 1.7  2002/11/15 01:45:53  scott
% can not -> cannot
%
% Revision 1.6  2002/10/16 21:42:25  arno
% default for highcuroff
%
% Revision 1.5  2002/08/12 16:25:02  arno
% inputdlg2
%
% Revision 1.4  2002/08/09 01:43:18  arno
% [6~[6~same
%
% Revision 1.3  2002/08/09 01:42:40  arno
% debugging filter over smal time period
%
% Revision 1.2  2002/08/09 00:41:22  arno
% updating for boundaries
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

function [EEG, com] = pop_eegfilt( EEG, locutoff, hicutoff, filtorder, revfilt);

com = '';
if nargin < 1
	help pop_eegfilt;
	return;
end;	
if isempty(EEG(1).data)
    disp('Pop_eegfilt() error: cannot filter an empty dataset'); return;
end;    
if nargin < 2
	% which set to save
	% -----------------
   	promptstr = { 'Highpass: lower edge of the frequency pass band (Hz) (0 -> lowpass)', ...
   				  'Lowpass: higher edge of the frequency pass band (Hz) (0 -> highpass)', ...
   				  strvcat('Notch filter the data. Give the notch range, i.e. [45 55] for 50 Hz)', ...
                  '(this option overwrites the low and high edge limits given above)'), ...
                  'Filter length in points (default: see >> help pop_eegfilt)' };
	inistr       = { '0', '0', '', '' };
   	result       = inputdlg2( promptstr, 'Filter the data -- pop_eegfilt()', 1,  inistr, 'pop_eegfilt');
	if size(result, 1) == 0 return; end;
	locutoff   	 = eval( result{1} );
	hicutoff 	 = eval( result{2} );
	if isempty( result{3} )
		 revfilt = 0;
	else 
        revfilt    = eval( [ '[' result{3} ']' ] );
        locutoff = revfilt(1);
        hicutoff = revfilt(2);
        revfilt = 1;
	end;
	if locutoff == 0 & hicutoff == 0 return; end;
	if isempty( result{4} )
		 filtorder = [];
	else filtorder    = eval( result{4} );
	end;
else
    if nargin < 3
        hicutoff = 0;
    end;
    if nargin < 4
        filtorder = [];
    end;
    if nargin < 5
        revfilt = 0;
    end;
end;

if locutoff & hicutoff
    disp('WARNING: BANDPASS FILTERS SOMETIMES DO NOT WORK')
    disp('WARNING: PLOT SPECTRUM AFTER FILTERING TO ASSESS FILTER EFFICIENCY')
    disp('WARNING: IF FILTER FAILS, LOWPASS DATA THEN HIGHPASS DATA')
end;
    
% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG com ] = eeg_eval( 'pop_eegfilt', EEG, 'warning', 'on', 'params', ...
                            { locutoff, hicutoff, filtorder, revfilt } );
    return;
end;

options = { EEG.srate, locutoff, hicutoff };
if ~isempty( filtorder )
	options = { options{:} filtorder };
else 
	options = { options{:} 0 };
end;
if revfilt ~= 0
	options = { options{:} revfilt };
end;

% warning
% -------
if exist('filtfilt') ~= 2
    disp('Warning: cannot find the signal processing toolbox');
    disp('         a simple fft/inverse fft filter will be used');
end;

if EEG.trials == 1 
	if ~isempty(EEG.event) & isfield(EEG.event, 'type') & isstr(EEG.event(1).type)
		boundaries = strmatch('boundary', {EEG.event.type});
		if isempty(boundaries)
            if exist('filtfilt') == 2
                EEG.data = eegfilt( EEG.data, options{:}); 
            else
                EEG.data = eegfiltfft( EEG.data, options{:}); 
            end;
		else
			options{4} = 0;
			disp('Pop_eegfilt:finding continuous data boundaries');
			tmplat = [ EEG.event.latency ];
            boundaries = tmplat(boundaries);
            boundaries = [0 floor(boundaries-0.49) EEG.pnts];
            try, warning off MATLAB:divideByZero
            catch, end;
			for n=1:length(boundaries)-1
				try
                    fprintf('Processing continuous data (%d:%d)\n',boundaries(n),boundaries(n+1)); 
                    if exist('filtfilt') == 2
                        EEG.data(:,boundaries(n)+1:boundaries(n+1)) = ...
                            eegfilt(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{:});
                    else
                        EEG.data(:,boundaries(n)+1:boundaries(n+1)) = ...
                            eegfiltfft(EEG.data(:,boundaries(n)+1:boundaries(n+1)), options{:});
                    end;
				catch
					fprintf('\nFilter error: continuous data portion too narrow (DC removed if highpass only)\n');
                    if locutoff ~= 0 & hicutoff == 0
                        tmprange = [boundaries(n)+1:boundaries(n+1)];
                        EEG.data(:,tmprange) = ...
                            EEG.data(:,tmprange) - repmat(mean(EEG.data(:,tmprange),2), [1 length(tmprange)]);
                    end;
				end;
			end
            try, warning on MATLAB:divideByZero
            catch, end;
		end
	else
        if exist('filtfilt') == 2
            EEG.data = eegfilt( EEG.data, options{:});
        else
            EEG.data = eegfiltfft( EEG.data, options{:});
        end;
	end;
else
    EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);
    if exist('filtfilt') == 2
        EEG.data = eegfilt( EEG.data, options{:});
    else
        EEG.data = eegfiltfft( EEG.data, options{:});
    end;
	% Note: reshape does not reserve new memory while EEG.data(:,:) does
end;	
EEG.icaact = [];

com = sprintf( '%s = pop_eegfilt( %s, %s, %s, [%s], [%s]);', inputname(1), inputname(1), ...
			num2str( locutoff), num2str( hicutoff), num2str( filtorder ), num2str( revfilt ) );
return;			

	 
