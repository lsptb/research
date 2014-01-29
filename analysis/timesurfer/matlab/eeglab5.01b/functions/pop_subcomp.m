% pop_subcomp() - remove specified components from an EEG dataset.
%                 and subtract their activities from the data. Else,
%                 remove components already marked for rejection.
% Usage:
%   >> OUTEEG = pop_subcomp( INEEG ); % pop-up window mode
%   >> OUTEEG = pop_subcomp( INEEG, components, confirm);
%
% Pop-up window interface:
%   "Component(s) to remove ..." - [edit box] Array of components to 
%                remove from the data. Sets the 'components' parameter 
%                in the command line call (see below).
%   "Component(s) to retain ..." - [edit box] Array of components to
%                to retain in the data. Sets the 'components' parameter in
%                the command line call. Then, comp_to_remove = ...
%                    setdiff([1:size(EEG.icaweights,1)], comp_to_keep)
%                Overwrites "Component(s) to remove" (above).
% Command line inputs:
%   INEEG      - Input EEG dataset.
%   components - Array of components to remove from the data. If empty, 
%                 remove components previously marked for rejection (e.g., 
%                 EEG.reject.gcompreject).
%   confirm    - [0|1] Display the difference between original and processed
%                dataset. 1 = Ask for confirmation. 0 = Do not ask. {Default: 0}
% Outputs:
%   OUTEEG     - output dataset.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: compvar()

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

% $Log: pop_subcomp.m,v $
% Revision 1.17  2006/01/31 20:19:18  arno
% options
%
% Revision 1.16  2006/01/25 21:27:30  arno
% nothing
%
% Revision 1.15  2006/01/25 21:26:15  arno
% add icachansind
%
% Revision 1.14  2006/01/25 21:24:16  arno
% fixing remove components
%
% Revision 1.13  2005/09/27 22:10:02  arno
% change default set name; allow to plot single trials
%
% Revision 1.12  2005/09/05 21:37:50  scott
% clarified trial-ERP and Confirmation window comments. -sm
%
% Revision 1.11  2003/12/24 19:40:58  scott
% edti help msg and text
%
% Revision 1.10  2003/02/19 19:17:25  arno
% update header for GUI
%
% Revision 1.9  2003/01/28 18:34:45  arno
% adding an option to keep components
%
% Revision 1.8  2003/01/16 18:10:56  arno
% debugging for PCA
%
% Revision 1.7  2002/08/12 18:36:35  arno
% questdlg2
%
% Revision 1.6  2002/08/12 02:22:53  arno
% inpudlg2
%
% Revision 1.5  2002/05/03 16:18:50  scott
% icaweight -> icaweights -sm
%
% Revision 1.4  2002/04/11 00:58:12  arno
% updating result size check
%
% Revision 1.3  2002/04/08 02:50:11  scott
% *** empty log message ***
%
% Revision 1.2  2002/04/08 02:48:32  scott
% worked on spelling, made text msg depend on number of pre-labeled comps -sm
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 02-15-02 propagate ica weight matrix -ad sm jorn 

function [EEG, com] = pop_subcomp( EEG, components, plotag )

com='';
if nargin < 1
   help pop_subcomp;
   return;
end;
if nargin < 3
	plotag = 0;
end;	
if nargin < 2
	% popup window parameters
	% -----------------------
	if ~isempty(EEG.reject.gcompreject)
        components = find(EEG.reject.gcompreject == 1);
        components = components(:)';
        promptstr    = { ['Component(s) to remove from the data ([] = marked comps.)'] };
        %promptstr    = { ['Components to subtract from data' 10 '(default: pre-labeled components to reject):'] };
    else
        components = [];
        promptstr    = { ['Component(s) to remove from data:'] };
    end;
    promptstr    = { ['Component(s) to remove from data:'] 'Component(s) to retain (overwrites "Component(s) to remove")' };
	inistr       = { int2str(components) '' };
	result       = inputdlg2( promptstr, 'Remove components from data -- pop_subcomp()', 1,  inistr, 'pop_subcomp');
	if length(result) == 0 return; end;
	components   = eval( [ '[' result{1} ']' ] );
    if ~isempty(result{2}), 
        components   = eval( [ '[' result{2} ']' ] );
        components  = setdiff([1:size(EEG.icaweights,1)], components);
    end;
end;
 
if isempty(components)
	if ~isempty(EEG.reject.gcompreject)
      		components = find(EEG.reject.gcompreject == 0);
   	else
        	fprintf('Warning: no components specified, no rejection performed\n');
         	return;
   	end;
else
    if (max(components) > EEG.nbchan) | min(components) < 1
        error('Component index out of range');
    end;
end;

fprintf('Computing projection ....\n');
eeglab_options; 
component_keep = setdiff(1:size(EEG.icaweights,1), components);
if option_computeica  
    compproj = EEG.icawinv(:, component_keep)*reshape(EEG.icaact(component_keep,:), length(component_keep), EEG.pnts*EEG.trials);
    %[ compproj, varegg ] = compvar( EEG.data, EEG.icaact, EEG.icawinv, setdiff(1:size(EEG.icaweights,1), components));
else
    compproj = EEG.icawinv(:, component_keep)*EEG.icaweights(component_keep,:)*EEG.icasphere ...
                 *reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.pnts*EEG.trials);
    %[ compproj, varegg ] = compvar( EEG.data, { EEG.icasphere EEG.icaweights }, EEG.icawinv, setdiff(1:size(EEG.icaweights,1), components));
end;    
compproj = reshape(compproj, EEG.nbchan, EEG.pnts, EEG.trials);

%fprintf( 'The ICA projection accounts for %2.2f percent of the data\n', 100*varegg);
	
if nargin < 2 | plotag ~= 0

    ButtonName = 'continue';
    while ~strcmpi(ButtonName, 'Cancel') & ~strcmpi(ButtonName, 'Accept')
        ButtonName=questdlg2( [ 'Please confirm. Are you sure you want to remove these components?' ], ...
                             'Confirmation', 'Cancel', 'Plot ERPs', 'Plot single trials', 'Accept', 'Accept');
        if strcmpi(ButtonName, 'Plot ERPs')
            if EEG.trials > 1
                tracing  = [ squeeze(mean(EEG.data,3)) squeeze(mean(compproj,3))];
                figure;   
                plotdata(tracing, EEG.pnts, [EEG.xmin*1000 EEG.xmax*1000 0 0], ...
                    'Trial ERPs (red) with and (blue) without these components');
            else
                warndlg2('Cannot plot ERPs for continuous data');
            end;
        elseif strcmpi(ButtonName, 'Plot single trials')  
        	eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Black = channel before rejection; red = after rejection -- eegplot()', ...
            	 'limits', [EEG.xmin EEG.xmax]*1000, 'data2', compproj); 
        end;
    end;    
    switch ButtonName,
        case 'Cancel', 
        	disp('Operation cancelled');
        	return; 
        case 'Accept',
       		disp('Components removed');
    end % switch
end;
EEG.data  = compproj;
EEG.setname = [ EEG.setname ' pruned with ICA'];
EEG.icaact = [];
EEG.icawinv    = EEG.icawinv(:,setdiff(1:size(EEG.icaweights,1), components));
EEG.icaweights = EEG.icaweights(setdiff(1:size(EEG.icaweights,1), components),:);

com = sprintf('%s = pop_subcomp( %s, [%s], %d);', inputname(1), inputname(1), ...
   int2str(components), plotag);
return;
