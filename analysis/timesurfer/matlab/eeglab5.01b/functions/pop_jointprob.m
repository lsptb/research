% pop_jointprob() - reject artifacts in an EEG dataset using joint 
%                   probability of the recorded electrode or component 
%                   activities observed at each time point.  e.g., Observing 
%                   large absoluate values at most electrodes or components 
%                   is improbable and may well mark the presence of artifact.
% Usage:
%   >> pop_jointprob( INEEG, typerej) % pop-up interative window mode
%   >> [OUTEEG, locthresh, globthresh, nrej] = ...
%		= pop_jointprob( INEEG, typerej, elec_comp, ...
%                   locthresh, globthresh, superpose, reject, vistype);
%
% Graphic interface:
%   "Electrode" - [edit box] electrode|component number(s) to take into
%                 consideration for rejection. Sets the 'elec_comp'
%                 parameter in the command line call (see below).
%   "Single-channel limit(s)" - [edit box] activity probability limit(s) (in 
%                 std. dev.) Sets the 'locthresh' command line parameter.
%                 If more than one, defined individual electrode|channel
%                 limits. If fewer values than the number of electrodes | 
%                 components specified above, the last input value is used 
%                 for all remaining electrodes|components.
%   "All-channel limit(s)" - [edit box] activity probability limit(s) (in std.
%                 dev.) for all channels (grouped). Sets the 'globthresh' 
%                 command line parameter.
%   "Display with previously marked rejections?" - [edit box] either YES or
%                 NO. Sets the command line option 'superpose'.
%   "Reject marked trial(s)?" - [edit box] either YES or NO. Sets the
%                 command line option 'reject'.
%   "visualization mode" - [edit box] either REJECTRIALS or EEGPLOT.
%                 Sets the command line option 'vistype'.
% Inputs:
%   INEEG      - input dataset
%   typerej    - [1|0] data to reject on (0 = component activations; 
%              1 = electrode data). {Default: 1 = electrode data}. 
%   elec_comp  - [n1 n2 ...] electrode|component number(s) to take into 
%              consideration for rejection
%   locthresh  - activity probability limit(s) (in std. dev.) See "Single-
%              channel limit(s)" above.
%   globthresh - global limit(s) (all activities grouped) (in std. dev.)
%   superpose  - [0|1] 0 = Do not superpose rejection marks on previously
%              marks stored in the dataset: 1 = Show both current and 
%              previous marks using different colors. {Default: 0}.
%   reject     - 0 = do not reject marked trials (but store the marks: 
%              1 = reject marked trials {Default: 1}.
%   vistype    - visualization mode: 0 = rejstatepoch(); 1 = eegplot()
%              {Default: 0}.  
%
% Outputs:
%   OUTEEG     - output dataset with updated joint probability array
%   locthresh  - electrodes probability of activity thresholds in terms
%              of standard-dev.
%   globthresh - global threshold (where all electrode activity are 
%              regrouped).
%   nrej       - number of rejected sweeps
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: jointprob(), rejstatepoch(), eegplot(), eeglab(), pop_rejepoch()  

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

% $Log: pop_jointprob.m,v $
% Revision 1.23  2006/01/31 20:20:26  arno
% options
%
% Revision 1.22  2004/05/07 22:09:37  arno
% debug rejecttrials
%
% Revision 1.21  2003/12/24 18:12:50  scott
% help msg and text edits
%
% Revision 1.20  2003/12/24 16:52:15  scott
% pop window text edits
%
% Revision 1.19  2003/12/17 18:14:38  arno
% default reject
%
% Revision 1.18  2003/12/04 23:18:26  arno
% adding code to reject trials
%
% Revision 1.17  2003/02/18 22:41:48  arno
% updating header for GUI
%
% Revision 1.16  2002/08/20 22:25:39  arno
% debug rej_statepoch call
%
% Revision 1.15  2002/08/14 00:51:50  arno
% debug multi-elec
%
% Revision 1.14  2002/08/12 21:51:17  arno
% same
%
% Revision 1.13  2002/08/12 21:50:49  arno
% text
%
% Revision 1.12  2002/08/12 02:42:10  arno
% inputdlg2
%
% Revision 1.11  2002/08/07 23:19:43  arno
% text
%
% Revision 1.10  2002/08/07 22:39:54  arno
% same
%
% Revision 1.9  2002/08/07 22:36:04  arno
% editing
%
% Revision 1.8  2002/08/07 22:33:42  arno
% editing header
%
% Revision 1.7  2002/07/31 00:04:10  arno
% debugging
%
% Revision 1.6  2002/07/30 23:21:55  arno
% debugging
%
% Revision 1.5  2002/07/30 22:27:55  arno
% new type of rejection
%
% Revision 1.4  2002/07/26 17:29:54  arno
% output command debug
%
% Revision 1.3  2002/07/26 17:21:03  arno
% return command if 2 outputs
%
% Revision 1.2  2002/07/26 16:40:54  arno
% switching icacomp
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 03-07-02 added srate argument to eegplot call -ad
% 03-08-02 add eeglab options -ad

function [EEG, locthresh, globthresh, nrej, com] = pop_jointprob( EEG, icacomp, elecrange, ...
                       		locthresh, globthresh, superpose, reject, vistype, topcommand);
com = '';
if nargin < 1
   help pop_jointprob;
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
	promptstr   = { [ fastif(icacomp, 'Electrode', 'Component') ' (number(s); Ex: 2 4 5):' ], ...
					[ fastif(icacomp, 'Single-channel', 'Single-component') ' limit(s) (std. dev(s).: Ex: 2 2 2.5):'], ...
					[ fastif(icacomp, 'All-channel', 'All-component') ' limit(s) (std. dev(s).: Ex: 2 2.1 2):'], ...
               		'Display previously marked rejections? (YES or NO)', ...
         			'Reject marked trial(s)? (YES or NO)', ...
         			'Visualization mode (REJECTRIALS|EEGPLOT)' };
	inistr      = { fastif(icacomp, ['1:' int2str(EEG.nbchan)], ['1:' int2str(size(EEG.icaweights,1))])...
					fastif(icacomp, '3', '5'),  ...
					fastif(icacomp, '3', '5'), ...
               		'YES', ...
            		'NO', ...
            		'REJECTTRIALS' };

	result       = inputdlg2( promptstr, fastif( ~icacomp, 'Reject. improbable comp. -- pop_jointprob()', 'Reject improbable data -- pop_jointprob()'), 1,  inistr, 'pop_jointprob');
	size_result  = size( result );
	if size_result(1) == 0 return; end;
	elecrange    = result{1};
	locthresh    = result{2};
	globthresh   = result{3};
	switch lower(result{4}), case 'yes', superpose=1; otherwise, superpose=0; end;
	switch lower(result{5}), case 'yes', reject=1; otherwise, reject=0; end;
	switch lower(result{6}), case 'rejecttrials', vistype=0; otherwise, vistype=1; end;
end;

if ~exist('vistype') vistype = 0; end;
if ~exist('reject') reject = 0; end;
if ~exist('superpose') superpose = 1; end;

if isstr(elecrange) % convert arguments if they are in text format 
	calldisp = 1;
	elecrange = eval( [ '[' elecrange ']' ]  );
	locthresh = eval( [ '[' locthresh ']' ]  );
	globthresh = eval( [ '[' globthresh ']' ]  );
else
	calldisp = 0;
end;

if isempty(elecrange)
	error('No electrode selectionned');
end;	

% compute the joint probability
% -----------------------------
if icacomp == 1
	fprintf('Computing joint probability for channels...\n');
    if isempty(EEG.stats.jpE)
		[ EEG.stats.jpE rejE ] = jointprob( EEG.data, locthresh, EEG.stats.jpE, 1); 
	end;
	[ tmp rejEtmp ] = jointprob( EEG.data(elecrange,:,:), locthresh, EEG.stats.jpE(elecrange,:), 1); 
    rejE    = zeros(EEG.nbchan, size(rejEtmp,2));
	rejE(elecrange,:) = rejEtmp;
	
	fprintf('Computing all-channel probability...\n');
	tmpdata = permute(EEG.data, [3 1 2]);
	tmpdata = reshape(tmpdata, size(tmpdata,1), size(tmpdata,2)*size(tmpdata,3));
	[ EEG.stats.jp rej ] = jointprob( tmpdata, globthresh, EEG.stats.jp, 1); 
else
	eeglab_options; % changed from eeglaboptions 3/30/02 -sm
	if option_computeica  
	    tmpdata = EEG.icaact;
	else
        tmpdata = (EEG.icaweights*EEG.icasphere)*reshape(EEG.data, EEG.nbchan, EEG.trials*EEG.pnts);
        tmpdata = reshape( tmpdata, size(tmpdata,1), EEG.pnts, EEG.trials);
    end;
	fprintf('Computing joint probability for components...\n');
    if isempty(EEG.stats.icajpE)
		[ EEG.stats.icajpE rejE ] = jointprob( tmpdata, locthresh, EEG.stats.icajpE, 1); 
	end;
	[ tmp rejEtmp ] = jointprob( tmpdata(elecrange,:), locthresh, EEG.stats.icajpE(elecrange,:), 1); 
    rejE    = zeros(size(tmpdata,1), size(rejEtmp,2));
	rejE(elecrange,:) = rejEtmp;

	fprintf('Computing global joint probability...\n');
	tmpdata2 = permute(tmpdata, [3 1 2]);
	tmpdata2 = reshape(tmpdata2, size(tmpdata2,1), size(tmpdata2,2)*size(tmpdata2,3));
	[ EEG.stats.icajp  rej] = jointprob( tmpdata2, globthresh, EEG.stats.icajp, 1); 
	clear tmpdata2;
end;
rej = rej' | max(rejE, [], 1);
fprintf('%d/%d trials marked for rejection\n', sum(rej), EEG.trials);

if calldisp
	if vistype == 1 % EEGPLOT -------------------------
	    if icacomp == 1 macrorej  = 'EEG.reject.rejjp';
	        			macrorejE = 'EEG.reject.rejjpE';
	    else			macrorej  = 'EEG.reject.icarejjp';
	        			macrorejE = 'EEG.reject.icarejjpE';
	    end;
		colrej = EEG.reject.rejjpcol;
		eeg_rejmacro; % script macro for generating command and old rejection arrays

	    if icacomp == 1
	        eegplot( EEG.data(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    else
	        eegplot( tmpdata(elecrange,:,:), 'srate', ...
		      EEG.srate, 'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, eegplotoptions{:}); 
	    end;	
    else % REJECTRIALS -------------------------
	  	if icacomp	== 1 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( EEG.data(elecrange,:,:), EEG.stats.jpE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.jp, ...
						'threshold', locthresh, 'thresholdg', globthresh, 'normalize', 'off'  );
		else 
			[ rej, rejE, n, locthresh, globthresh] = ... 
				rejstatepoch( tmpdata(elecrange,:,:), EEG.stats.icajpE(elecrange,:), 'global', 'on', 'rejglob', EEG.stats.icajp, ...
						'threshold', locthresh, 'thresholdg', globthresh, 'normalize', 'off' );
		end;		
		nrej = n;
	end;	
else
	% compute rejection locally
	rejtmp = max(rejE(elecrange,:),[],1);
	rej = rejtmp | rej;
	nrej =  sum(rej);
	fprintf('%d trials marked for rejection\n', nrej);
end;
if ~isempty(rej)
	if icacomp	== 1
		EEG.reject.rejjp = rej;
		EEG.reject.rejjpE = rejE;
	else
		EEG.reject.icarejjp = rej;
		EEG.reject.icarejjpE = rejE;
	end;
    if reject
        EEG = pop_rejepoch(EEG, rej, 0);
    end;
end;
nrej = sum(rej);

com = [ com sprintf('%s = pop_jointprob(%s,%s);', inputname(1), ...
		inputname(1), vararg2str({icacomp,elecrange,locthresh,globthresh,superpose,reject})) ]; 
if nargin < 3 & nargout == 2
	locthresh = com;
end;

return;
