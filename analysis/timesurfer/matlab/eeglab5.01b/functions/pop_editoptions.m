% pop_editoptions() - Edit memory-saving eeglab() options. These are stored in 
%                     a file 'eeg_options.m'. With no argument, pop up a window 
%                     to allow the user to set/unset these options. Store
%                     user choices in a new 'eeg_options.m' file in the 
%                     working directory.
%
% Usage: >> pop_editoptions;
%        >> pop_editoptions( 'key1', value1, 'key2', value2, ...);
%
% Graphic interface inputs:
%   "If set, keep at most one dataset in memory ..." - [checkbox] If set, EEGLAB will only retain the current
%                   dataset in memory. All other datasets will be automatically
%                   read and writen to disk. All EEGLAB functionalities are preserved
%                   even for dataset stored on disk. 
%   "If set, write data in same file as dataset ..." - [checkbox] Set -> dataset data (EEG.data) are 
%                   saved in the EEG structure in the standard Matlab dataset (.set) file. 
%                   Unset -> The EEG.data are saved as a transposed stream of 32-bit 
%                   floats in a separate binary file. As of Matlab 4.51, the order 
%                   of the data in the binary file is as in the transpose of EEG.data 
%                   (i.e., as in EEG.data', frames by channels). This allows quick 
%                   reading of single channels from the data, e.g. when comparing 
%                   channels across datasets. The stored files have the extension 
%                   .dat instead of the pre-4.51, non-transposed .fdt. Both file types 
%                   are read by the dataset load function. Command line equivalent is
%                   option_savematlab.
%   "Precompute ICA activations" - [checkbox] If set, all the ICA activation
%                   time courses are precomputed (this requires more RAM). 
%                   Command line equivalent: option_computeica.
%   "If set, remember old folder when reading dataset" - [checkbox] this option
%                   is convinient if the file you are working on are not in the 
%                   current folder.
%
% Commandline keywords:
%   'option_computeica' - [0|1] If 1, compute the ICA component activitations and
%                   store them in a new variable. If 0, compute ICA activations
%                   only when needed (& only partially, if possible) and do not
%                   store the results).
%   NOTE: Turn OFF the options above when working with very large datasets or on 
%                   computers with limited memory.
%   'option_savematlab' - [0|1] If 1, datasets are saved as single Matlab .set files. 
%                   If 0, dataset data are saved in separate 32-bit binary float 
%                   .dat files.  See the corresponding GUI option above for details. 
% Outputs:
%   In the output workspace, variables 'option_computeica', 
%   and 'option_savematlab'  are updated, and a new 'eeg_options.m' file may be
%   written to the working directory. The copy of 'eeg_options.m' placed in your 
%   working directory overwrites system defaults whenever EEGLAB operates in this
%   directory (assuming your working directory is in your MATLABPATH - see path()).
%   To adjust these options system-wide, edit the master "eeg_options.m" file in the
%   EEGLAB directory heirarchy.
%
% Author: Arnaud Delorme, SCCN / INC / UCSD, March 2002
%
% See also: eeg_options()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Arnaud Delorme, CNL / Salk Institute, 09 March 2002, arno@salk.edu
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

% $Log: pop_editoptions.m,v $
% Revision 1.43  2006/03/19 05:29:16  toby
% spelling edit
%
% Revision 1.42  2006/03/11 05:39:03  arno
% header
%
% Revision 1.41  2006/03/11 00:32:08  arno
% header
%
% Revision 1.40  2006/03/06 21:50:28  arno
% remove file after testing
%
% Revision 1.39  2006/02/03 18:53:47  arno
% default option file
%
% Revision 1.38  2006/01/31 20:11:28  arno
% editing option etc...
%
% Revision 1.37  2006/01/31 19:47:31  arno
% change default folder etc...
%
% Revision 1.36  2006/01/31 19:33:39  arno
% adding filepath to new options
%
% Revision 1.35  2006/01/31 19:06:03  arno
% use new features of eeg_readoptions
%
% Revision 1.34  2006/01/31 18:51:19  arno
% now dealing with structure not variables
%
% Revision 1.33  2006/01/31 18:35:05  arno
% change read options
%
% Revision 1.32  2006/01/31 00:15:38  arno
% fixing last changes
%
% Revision 1.31  2006/01/31 00:13:35  arno
% same
%
% Revision 1.30  2006/01/31 00:12:49  arno
% fix checking number of dataset loaded
%
% Revision 1.29  2006/01/31 00:01:52  arno
% adding callback
%
% Revision 1.28  2005/08/04 17:40:40  arno
% debug saving new option file
%
% Revision 1.27  2005/08/04 15:39:22  arno
% edit header
%
% Revision 1.26  2005/08/02 23:24:14  arno
% process comments etc...
%
% Revision 1.25  2005/08/02 01:54:07  arno
% debug header etc...
%
% Revision 1.24  2005/08/02 01:45:39  arno
% dealing with newer option files
%
% Revision 1.23  2005/08/01 14:36:03  arno
% debug wrong alignment etc...
%
% Revision 1.22  2005/07/31 23:01:56  arno
% debuging option writing
%
% Revision 1.21  2004/11/22 17:30:23  scott
% edited help message - now tells the whole story
%
% Revision 1.20  2004/11/21 02:45:57  scott
% help message - new transposed binary file option
%
% Revision 1.19  2003/07/31 22:28:41  arno
% *** empty log message ***
%
% Revision 1.18  2003/04/10 17:31:03  arno
% header edit
%
% Revision 1.17  2002/12/04 19:13:19  arno
% debuging for MacOS6 Matlab 6 directories
%
% Revision 1.16  2002/11/15 01:38:52  scott
% same
%
% Revision 1.15  2002/10/23 15:01:50  arno
% isppc -> computer
%
% Revision 1.14  2002/09/26 14:00:44  scott
% help msg -sm
%
% Revision 1.13  2002/08/20 22:36:47  arno
% debug for windows
%
% Revision 1.12  2002/08/19 21:55:40  arno
% add Mac statement
%
% Revision 1.11  2002/08/14 21:30:48  arno
% debug for windows
%
% Revision 1.10  2002/08/13 23:47:11  arno
% debugging message
%
% Revision 1.9  2002/08/13 18:25:57  scott
% help msg
%
% Revision 1.8  2002/08/13 16:10:43  arno
% debugging
%
% Revision 1.7  2002/08/13 00:31:31  scott
% text
%
% Revision 1.6  2002/08/12 18:40:17  arno
% questdlg2
%
% Revision 1.5  2002/07/25 17:22:31  arno
% adding a clear functions statement
%
% Revision 1.4  2002/04/25 16:56:13  arno
% copying file if read only
%
% Revision 1.3  2002/04/21 01:07:59  scott
% edited help msg -sm
%
% Revision 1.2  2002/04/18 18:09:57  arno
% updating error message
%
% Revision 1.1  2002/04/05 17:46:04  jorn
% Initial revision
%
%02/19/2002 debuging function -ad

function com = pop_editoptions(varargin);

com = '';

datasets_in_memory = 0;
if nargin > 0
    if ~isstr(varargin{1})
        datasets_in_memory = varargin{1};
        varargin = {};
    end;
end;

% parse the eeg_options file
% ----------------------------
eeglab_options;
filename = which('eeg_options.m');
fid = fopen( filename, 'r+');
storelocal = 0;
if	fid == -1
	if exist(filename) == 2 
		if ~popask(['Cannot modify read-only file ''' filename '''' 10 'Do you want to store the option file somewhere else ?']);
			return;
		else 
            warndlg2(strvcat('Warning: you must store the file in a FOLDER always accessible', ...
                             'from Matlab (i.e. a folder in the Matlab path) and not necessarily in', ...
                             'the local folder. Otherwise, every time you restart EEGLAB, the default', ...
                             'EEGLAB options will apply (the path you choose will be added temporarily', ...
                             'for this session). Select a folder in the next pop-up file window.'), 'Warning');
            try
                filepath = uigetdir('', 'Pick a Directory');
            catch,
                [filename filepath] = uigetfile('*.m', 'Pick a folder', 'eeg_options');
            end;
		end;
        if filepath(1) == 0, return; end;
        
        % see if the folder can be written into
        % -------------------------------------
        fid = fopen( fullfile(filepath, 'eeg_options.m'), 'w');
        if fid == -1
            error('Cannot write into this folder');
        end;
        fclose(fid);
        delete(fullfile(filepath, 'eeg_options.m'));
        
        % change default folder option
        % ----------------------------
        W_MAIN = findobj('tag', 'EEGLAB');
        if ~isempty(W_MAIN)
            tmpuserdata    = get(W_MAIN, 'userdata');
            tmpuserdata{3} = filepath;
            set(W_MAIN, 'userdata', tmpuserdata);
        end;
        
        % read variables values and description
        % --------------------------------------
        [ header opt ] = eeg_readoptions( 'eeg_optionsbackup.m' ); 
	else
		error('File not found');
	end;
else 
    [filepath filename ext] = fileparts(filename);
    fprintf('Using option file in directory %s\n', filepath);
    
    % read variables values and description
    % --------------------------------------
    [ header opt ] = eeg_readoptions( 'eeg_optionsbackup.m' ); 
    [ header opt ] = eeg_readoptions( fid, opt  ); % use opt from above as default
end;

if nargin < 2
    geometry = { [6 1] };
    uilist = { ...
         { 'Style', 'text', 'string', '', 'fontweight', 'bold'  }, ...
         { 'Style', 'text', 'string', 'Set/Unset', 'fontweight', 'bold'   } };

    % add all fields to graphic interface
    % -----------------------------------
    for index = 1:length(opt)
        % format the description to fit a help box
        % ----------------------------------------
        descrip = { 'string', opt(index).description }; % strmultiline(description{ index }, 80, 10) };
           
        % create the gui for this variable
        % --------------------------------
        geometry = { geometry{:} [4 0.3 0.1] };
        if strcmpi(opt(index).varname, 'option_storedisk') & datasets_in_memory
            cb_nomodif = [ 'set(gcbo, ''value'', ~get(gcbo, ''value''));' ...
                           'warndlg2(strvcat(''This option may only be modified when at most one dataset is stored in memory.''));' ];
            
        else
            cb_nomodif = '';
        end;
        
        if ~isempty(opt(index).value)
            uilist   = { uilist{:}, { 'Style', 'text', descrip{:}, 'horizontalalignment', 'left' }, ...
                         { 'Style', 'checkbox', 'string', '    ', 'value', opt(index).value 'callback' cb_nomodif } { } }; 
        else
            uilist   = { uilist{:}, { 'Style', 'text', descrip{:}, 'fontweight' 'bold', 'horizontalalignment', 'left' }, ...
                         { } { } }; 
        end;
    end;

    results = inputgui( geometry, uilist, 'pophelp(''pop_editoptions'');', 'Memory options - pop_editoptions()', ...
                        [], 'normal');
    if length(results) == 0, return; end;
   
    % decode inputs
    % -------------
    args = {};
    count = 1;
    for index = 1:length(opt)
        if ~isempty(opt(index).varname)
            args = {  args{:}, opt(index).varname, results{count} }; 
            count = count+1;
        end;
    end;
else 
    % no interactive inputs
    % ---------------------
    args = varargin;
end;

% decode inputs
% -------------
for index = 1:2:length(args)
    ind = strmatch(args{index}, { opt.varname }, 'exact');
    if isempty(ind)
        error(['Variable name ''' args{index} ''' is invalid']);
    else
        opt(ind).value = args{index+1};
    end;
end;        

% write to eeg_options file
% -------------------------
fid = fopen( fullfile(filepath, 'eeg_options.m'), 'w');
addpath(filepath);
if fid == -1
	error('File error, check writing permission');
end;
fprintf(fid, '%s\n', header);
for index = 1:length(opt)
    if isempty(opt(index).varname)
        fprintf( fid, '%% %s\n', opt(index).description);
    else
        fprintf( fid, '%s = %d ; %% %s\n', opt(index).varname, opt(index).value, opt(index).description);
    end;
end;
fclose(fid);    

% generate the output text command
% --------------------------------
com = 'pop_editoptions(';
for index = 1:2:length(args)
    com = sprintf( '%s ''%s'', %d,', com, args{index}, args{index+1});
end;
com = [com(1:end-1) ');'];   
clear functions

% ---------------------------
function  chopedtext = choptext( tmptext )
    chopedtext = '';
    while length(tmptext) > 30
          blanks = findstr( tmptext, ' ');
          [tmp I] = min( abs(blanks - 30) );
          chopedtext = [ chopedtext ''' 10 ''' tmptext(1:blanks(I)) ];
          tmptext  = tmptext(blanks(I)+1:end);
    end;    
    chopedtext = [ chopedtext ''' 10 ''' tmptext];
    chopedtext = chopedtext(7:end);
return;

function num = popask( text )
	 ButtonName=questdlg2( text, ...
	        'Confirmation', 'Cancel', 'Yes','Yes');
	 switch lower(ButtonName),
	      case 'cancel', num = 0;
	      case 'yes',    num = 1;
	 end;
