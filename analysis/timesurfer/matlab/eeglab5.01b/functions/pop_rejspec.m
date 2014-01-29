% pop_rejspec() - rejection of artifact in a dataset using 
%                 thresholding of frequencies in the data.
% Usage:
%   >>  pop_rejspec(INEEG, typerej); % pop-up interactive windo mode
%   >> [OUTEEG, Indexes] = pop_rejspec( INEEG, typerej, elec_comp, ...
%         lowthresh, upthresh, startfreq, endfreq, superpose, reject);
%
% Pop-up window options:
%   "Electrode|Component" - [edit box] electrode or component number(s) to 
%                 take into consideration for rejection. Sets the 'elec_comp'
%                 parameter in the command line call (see below).
%   "Lower limits(s)" - [edit box] lower threshold limits(s) (in dB). 
%                 Sets the command line parameter 'lowthresh'. If more than
%                 one, apply to each electrode|component individually. If
%                 fewer than number of electrodes|components, apply the
%                 last values to all remaining electrodes|components.
%   "Upper limits(s)" - [edit box] upper threshold limit(s) in dB. 
%                 Sets the command line parameter 'upthresh'.
%   "Low frequency(s)" - [edit box] low-frequency limit(s) in Hz. 
%                 Sets the command line parameter 'startfreq'.
%   "High frequency(s)" - [edit box] high-frequency limit(s) in Hz. 
%                 Sets the command line parameter 'endfreq'.
%   "Display previous rejection marks?" - [edit box] either YES or NO. 
%                  Sets the command line input option 'superpose'.
%   "Reject marked trials?" - [edit box] either YES or NO. Sets the
%                 command line input option 'reject'.
%
% Command line inputs:
%   INEEG      - input dataset
%   typerej    - [1|0] data to reject on (0 = component activations; 1 = 
%              electrode data). {Default is 1}. 
%   elec_comp  - [e1 e2 ...] electrode|component number(s) to take into 
%              consideration during rejection
%   lowthresh  - lower threshold limit(s) in dB. Can be an array if 
%              several electrodes|components. If fewer values than number 
%              of electrodes|components, the last value is used for the 
%              remaining electrodes|components.
%   upthresh  - upper threshold limit(s) in dB (same syntax as lowthresh)
%   startfreq  - low frequency limit(s) in Hz (same syntax  as lowthresh)
%   endfreq    - high frequency limit(s) in Hz (same syntax  as lowthresh).
%              Options 'startfreq' and 'endfreq' define the frequncy range 
%              used during rejection.
%   superpose  - [0|1] 0 = Do not superpose rejection marks on previous
%              marks stored in the dataset. 1 = Show both previous and
%              current marks using different colors. {Default: 0}.
%   reject     - [0|1] 0 = Do not reject marked trials (but store the 
%              marks. 1 = Reject marked trials. {Default: 1}.
%
% Outputs:
%   OUTEEG     - output dataset with updated spectrograms
%   Indexes    - index of rejected sweeps
%   Note: When eegplot() is called, modifications are applied to the current 
%   dataset at the end of the call to eegplot() (e.g., when the user presses 
%   the 'Reject' button).
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: eegthresh(), eeglab(), eegplot(), pop_rejepoch()

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

% $Log: pop_rejspec.m,v $
% Revision 1.25  2006/01/31 20:16:02  arno
% options
%
% Revision 1.24  2005/08/04 15:37:52  arno
% remove memory saving option
%
% Revision 1.23  2003/12/24 18:56:16  scott
% help msg and pop-window text edits
%
% Revision 1.22  2003/12/17 18:16:35  arno
% reject
%
% Revision 1.21  2003/12/04 23:27:46  arno
% reject trials
%
% Revision 1.20  2003/02/18 22:50:06  arno
% updating header for GUI
%
% Revision 1.19  2003/01/10 01:10:45  arno
% change default position
%
% Revision 1.18  2002/11/15 19:07:24  arno
% updating default position (eegplot)
%
% Revision 1.17  2002/11/15 01:37:17  arno
% header typo
%
% Revision 1.16  2002/11/13 02:06:52  arno
% debugging command line call
%
% Revision 1.15  2002/11/13 01:21:24  arno
% typo in header
%
% Revision 1.14  2002/11/12 23:48:59  luca
% noew saving output
%
% Revision 1.13  2002/08/14 15:26:03  arno
% debug for few channels
%
% Revision 1.12  2002/08/12 21:53:00  arno
% text
%
% Revision 1.11  2002/08/12 02:32:24  arno
% inputdlg2
%
% Revision 1.10  2002/08/07 22:38:42  arno
% editing header
%
% Revision 1.9  2002/07/31 01:02:13  arno
% debugging
%
% Revision 1.8  2002/07/30 23:37:51  arno
% debugging
%
% Revision 1.7  2002/07/30 23:33:47  arno
% new rejection type
%
% Revision 1.6  2002/07/30 15:17:52  arno
% debug display
%
% Revision 1.5  2002/07/30 15:17:16  arno
% debugging
%
% Revision 1.4  2002/07/26 18:09:31  arno
% debugging
%
% Revision 1.3  2002/07/26 18:04:48  arno
% debugging
%
% Revision 1.2  2002/07/26 16:44:35  arno
% switching icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-08-02 reworked spectrum to save space & add eeglab options -ad

function [EEG, I1, com] = pop_rejspec( EEG, icacomp, elecrange, negthresh, posthresh, ...
   						startfreq, endfreq, superpose, reject, topcommand);

I1 = [];
com = '';
if nargin < 1
   help pop_rejspec;
   return;
end;  
if nargin < 2
   icacomp = 1;
end;  

if icacomp == 0
	if isempty( EEG.icasphere )
	    ButtonName=questdlg( 'Do you want to run ICA now ?', ...
                         'Confirmation', 'NO', 'YES', 'YES');
    	switch ButtonName,
        	case 'NO', disp('Operation cancelled'); return;   
        	case 'YES', [ EEG com ] = pop_runica(EEG);
    	end % switch
	end;
end;	
if exist('reject') ~= 1
    reject = 1;
end;
if nargin < 3

	% which set to save
	% -----------------
	promptstr   = { fastif(icacomp==0, 'Component number(s) (Ex: 2 4 5):', ...
                                           'Electrode number(s) (Ex: 2 4 5):'), ...
					'Lower limit(s) (dB):', ...
					'Upper limit(s) (dB):', ...
					'Low frequency(s) (Hz):', ...
					'High frequency(s) (Hz):', ...
               		'Display previous rejection marks? (YES or NO)', ...
         			'Reject marked trial(s)? (YES or NO)' };
	inistr      = { ['1:' int2str(EEG.nbchan)], ...
					'-10', ...
					'10', ...
					'15', ...
					'30', ...
               		'NO', ...
            		'NO' };

	result       = inputdlg2( promptstr, fastif(~icacomp, 'Reject by component spectra -- pop_rejspec()', ...
											   'Reject by data spectra -- pop_rejspec()'), 1,  inistr, 'pop_rejspec');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	negthresh    = result{2};
	posthresh    = result{3};
	startfreq    = result{4};
	endfreq      = result{5};
	switch lower(result{6}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{7}), case 'yes', reject=1; otherwise, reject=0; end;
end;

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	negthresh = eval( [ '[' negthresh ']' ]  );
	posthresh = eval( [ '[' posthresh ']' ]  );
	startfreq = eval( [ '[' startfreq ']' ]  );
	endfreq   = eval( [ '[' endfreq ']' ]  );
else
	calldisp = 0;
end;

sizewin = 2^nextpow2(EEG.pnts);
if icacomp == 1
    [allspec, Irej, tmprejE, freqs ] = spectrumthresh( EEG.data, EEG.specdata, ...
							elecrange, EEG.srate, negthresh, posthresh, startfreq, endfreq);
    rejE = zeros(EEG.nbchan, EEG.trials);
    rejE(elecrange,Irej) = tmprejE;
else
    % test if ICA was computed
    % ------------------------
    eeglab_options; % changed from eeglaboptions 3/30/02 -sm
 	if option_computeica  
    	icaacttmp = EEG.icaact(elecrange, :, :);
	else
        icaacttmp = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, ...
												  EEG.nbchan, EEG.trials*EEG.pnts);
        icaacttmp = reshape( icaacttmp, size(icaacttmp,1), EEG.pnts, EEG.trials);
    end;
    [allspec, Irej, tmprejE, freqs ] = spectrumthresh( icaacttmp, EEG.specicaact, elecrange, ...
                                                      EEG.srate, negthresh, posthresh, startfreq, endfreq);
    rejE = zeros(EEG.nbchan, size(icaacttmp,1));
    rejE(elecrange,Irej) = tmprejE;
end;

fprintf('%d channel selected\n', size(elecrange(:), 1));
fprintf('%d/%d trials marked for rejection\n', length(Irej), EEG.trials);
rej = zeros( 1, EEG.trials);
rej(Irej) = 1;

if calldisp
	nbpnts = size(allspec,2);
    if icacomp == 1 macrorej  = 'EEG.reject.rejfreq';
        			macrorejE = 'EEG.reject.rejfreqE';
    else			macrorej  = 'EEG.reject.icarejfreq';
        			macrorejE = 'EEG.reject.icarejfreqE';
    end;
	colrej = EEG.reject.rejfreqcol;
	eeg_rejmacro; % script macro for generating command and old rejection arrays
	if icacomp == 1
		eegplot(EEG.data(elecrange,:,:), 'winlength', 5, 'position', [100 550 800 500], ...
				'limits', [EEG.xmin EEG.xmax]*1000, 'xgrid', 'off', 'tag', 'childEEG' );
	else
		eegplot(icaacttmp(elecrange,:,:), 'winlength', 5, 'position', [100 550 800 500], 'limits', ...
				[EEG.xmin EEG.xmax]*1000 , 'xgrid', 'off', 'tag', 'childEEG' );
	end;	
	eegplot( allspec(elecrange,:,:), 'srate', EEG.srate, 'freqlimits', [1 EEG.srate/2], 'command', ...
			 command, 'children', gcf, 'position', [100 50 800 500], eegplotoptions{:}); 
end;
if ~isempty(rej)
	if icacomp	== 1
		EEG.reject.rejfreq = rej;
		EEG.reject.rejfreqE = rejE;
	else
		EEG.reject.icarejfreq = rej;
		EEG.reject.icarejfreqE = rejE;
	end;
    if reject
        EEG = pop_rejepoch(EEG, rej, 0);
    end;
end;

% store variables
% ---------------
if icacomp == 1, EEG.specdata = allspec;
else,            EEG.specicaact = allspec;
end;
    
com = [com sprintf('%s = pop_rejspec( %s, %s);', inputname(1), ...
   inputname(1), vararg2str({icacomp, elecrange,  negthresh, posthresh, startfreq, endfreq, superpose, reject })) ]; 
if nargin < 3 & nargout == 2
	I1 = com;
end;

return;

% compute spectrum and reject artifacts
% -------------------------------------
function [specdata, Irej, Erej, freqs ] = spectrumthresh( data, specdata, elecrange, srate, negthresh, posthresh, startfreq, endfreq);
	% compute the fft if necessary - old version
	%freqs = EEG.srate*[1, sizewin]/sizewin/2;
	%EEG.specdata = fft( EEG.data-repmat(mean(EEG.data,2), [1 EEG.pnts 1]), sizewin, 2);
	%EEG.specdata = EEG.specdata( :, 2:sizewin/2+1, :);
	%EEG.specdata = 20*log2(abs( EEG.specdata ).^2);

	[tmp freqs] = pmtm( data(1,:,1) ); % just to get the frequencies 	
    freqs = [freqs(1) freqs(end)]*srate;	
    
	fprintf('Computing spectrum (using slepian tapers; done only once):\n');
	if isempty(specdata)
		for index = 1:size(data,1)
			fprintf('%d ', index);    
			for indextrials = 1:size(data,3)
				[ tmpspec(index,:,indextrials) freqs] = pmtm( data(index,:,indextrials) );
			end;
			freqs = [freqs(1) freqs(end)]*srate;	
		end;
		tmpspec  = 10*log(tmpspec);
		tmpspec  = tmpspec - repmat( mean(tmpspec,3), [1 1 size(data,3)]);
		specdata = tmpspec;
	end;
	
	% perform the rejection
	% ---------------------	
	[I1 Irej NS Erej] = eegthresh( specdata(elecrange, :, :), size(specdata,2), 1:length(elecrange), negthresh, posthresh, ...
										 freqs, startfreq, max(max(freqs), endfreq));
	fprintf('\n');    
	
