% pop_ctf_read() - read CTF file as EEGLAB dataset
%
% Usage:
%   >> EEGOUT = pop_ctf_read; % pop up graphic interface
%   >> EEGOUT = pop_ctf_read(folder);
%   >> EEGOUT = pop_ctf_read(folder, chans, time, trials);
%
% Inputs:
%   folder   - [string]  EEGLAB figure
%   chans    - [integer array or string] see ctf_read()
%   time     - [float array or string] see ctf_read()
%   trials   - [integer array or string] see ctf_read()
%
%
% Author: Arnaud Delorme (SCCN, UCSD) and Daren Weber (DNL, UCSF)
%
% See also: ctf_read(), ctf_readmarkerfile()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 2003 Arnaud Delorme, SCCN, UCSD, arno@salk.edu
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

% $Log: pop_ctf_read.m,v $
% Revision 1.10  2005/12/13 18:56:43  psdlw
% Resolved EEG.comment initialization and added comments about EEG.filename and EEG.filepath fields.
%
% Revision 1.9  2005/10/27 05:34:17  arnodelorme
% change file name option in pop_read_ctf.m
%
% Revision 1.8  2004/07/28 15:44:59  arnodelorme
% remove command line debug msg
%
% Revision 1.7  2004/07/28 15:41:25  arnodelorme
% command line history for pop_ctf_read
%
% Revision 1.5  2004/07/06 23:42:21  psdlw
% bug fixing an error when no marker file is present
%
% Revision 1.4  2004/06/15 17:36:08  arnodelorme
% updated reading files in 4 functions, pop_ctf_read.m, ctf_read_markerfile.m, ctf_read_meg4.m, ctf_read_res4.m + lots of minor debuging
%
% Revision 1.3  2004/06/03 00:15:03  psdlw
% updated to handle ctf.data as a 3D matrix rather than a cell array of trials
%
% Revision 1.2  2004/04/02 17:16:41  arnodelorme
% Read this file to set up a menu for this toolbox in EEGLAB
%

% from ctf2eeglab - script to convert and save ctf .ds into eeglab .set data
%
% The script uses a ctf struct in the matlab workspace or a GUI prompt to
% load a CTF .ds folder, then converts the ctf data into an EEGLAB EEG
% struct, saving the resulting dataset into an EEEGLAB .set file, located
% in the same path as the ctf .ds folder.  The GUI prompt for the CTF .ds
% folder also provides access to definition of the channels, time and
% trials to load.
%
% Licence:  GNU GPL, no express or implied warranties
% Modified: 01/2004, Darren.Weber_at_radiology.ucsf.edu
%                    - developed in collaboration with Fredrick Carver of
%                    the NIH, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [EEG, com] = pop_ctf_read(orifolder, varargin)

  EEG = [];
  com = '';
  if nargin < 1
      [tmp orifolder] = uigetfile('*', 'Pick any file in CTF folder - pop_ctf_read()');
      if tmp == 0 return; end;
      orifold = pwd;
      [ folder parentfold ] = formatfolder(orifolder);
      cd(parentfold);
      
      % read info and prompt
      % --------------------
      listchan = { 'all' 'eeg' 'meg' 'ref' 'other' 'megeeg'};
      disp('Reading file info...');
      ctf = ctf_read_res4(folder, 0);
      uigeom       = { [2 1] [2 1] [2 1] [2 1] };
      uilist       = { { 'style' 'text' 'string' 'Channels group subset' } ...
                     { 'style' 'list' 'string' strvcat(listchan{:}) } ...
                     { 'style' 'text' 'string' 'or channel indices (overwrite above)' } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' 'Time range (default all)' } ...
                     { 'style' 'edit' 'string' '' } ...
                     { 'style' 'text' 'string' sprintf('Trial range (default all [1:%d])', ctf.setup.number_trials) } ...
                     { 'style' 'edit' 'string' '' } };
      result = inputgui( uigeom, uilist, 'pophelp(''pop_ctf_read'')', 'Load a CTF dataset');
      if length( result ) == 0 return; end;
      
      % decode inputs
      % -------------
      options{1} = [1:ctf.setup.number_channels];
      if result{1} > 1
          options{1} = listchan{result{1}};
      end;
      if ~isempty(result{2})
          options{1} = eval( [ '[' result{2} ']' ] );
      end;
      if ~isempty(result{3})
          options{2} = eval( [ '[' result{3} ']' ] );
      end;
      if ~isempty(result{4})
          options{3} = eval( [ '[' result{4} ']' ] );
      end;
  else 
      options = varargin;
  end;
  if exist('orifold') ~= 1
      orifold = pwd;
  end;
  [folder parentfold] = formatfolder(orifolder);
  cd(parentfold);
  
  alltrials = [];
  if length(options) == 3
      alltrials = options{3};
  end;
  if length(options) > 1 & isempty(options{2})
      options{2} = 'all';
  end;
  
  % read the data
  % -------------
  ctf = ctf_read(folder,options{:});
  
  % check if the data is averaged
  % -----------------------------
  if ctf.setup.number_trials_averaged > 0,
      warning('this .ds folder is averaged, removing stdev');
  end

  % check if the data is greater than 500 Mb
  data_size = ctf.setup.number_samples * ctf.setup.number_channels * ctf.setup.number_trials;
  data_bytes = data_size * 8;
  if data_bytes > 5e9, warning('data is greater than 500 Mb'); end
  clear data_size data_bytes;

  % ctf.data is a 3D matrix of samples(time) x channels x trials
  % whereas EEGLAB has a 3D data matrix with channels X samples X trials
  data = zeros( [ size(ctf.data,2) size(ctf.data,1) size(ctf.data,3) ] );
  for i = 1:size(ctf.data,3),
      data(:,:,i) = ctf.data(:,:,i)';
  end
  ctf.data = [];
  
  % import the data into the EEGLAB EEG struct
  
  [DSpath,DSfile,DSext] = fileparts(ctf.folder);
  
  EEG = eeg_emptyset;
  EEG.setname = DSfile;
  
  % ---
  % These fields now contain the name of the dataset *once*
  % it has been saved (so reamin empty before the dataset 
  % has been saved).
  %EEG.filename = [DSfile,'.set'];
  %EEG.filepath = DSpath;
  % ---

  EEG.comments = [ 'Original folder: ' ctf.folder ];
  %EEG.comments = ctf.setup.run_description';
  EEG.pnts = ctf.setup.number_samples;
  EEG.nbchan = ctf.setup.number_channels;
  EEG.trials = ctf.setup.number_trials;
  EEG.srate = ctf.setup.sample_rate;
  EEG.xmin = ctf.setup.start_sec;
  EEG.xmax = ctf.setup.end_sec;
  EEG.data = data;
  EEG.ref = 'common';

  for i=1:ctf.setup.number_channels,
      EEG.chanlocs(i).labels = ctf.sensor.label{i};
      EEG.chanlocs(i).X      = ctf.sensor.location(1,i);
      EEG.chanlocs(i).Y      = ctf.sensor.location(2,i);
      EEG.chanlocs(i).Z      = ctf.sensor.location(3,i);
  end
  EEG.chanlocs = convertlocs(EEG.chanlocs, 'cart2all');
  
  % now clear the workspace of the input data
  % -----------------------------------------
  clear data
  clear DSpath DSfile DSext i

  % import event information
  % ------------------------
  %try
      eventarray  = [];
      allfields   = {};
      timefields  = {};
      otherfields = {};
      
      eventstruct = ctf_read_markerfile(ctf.folder, ctf);
      
      if isfield(eventstruct,'markers'),
          if ~isempty(eventstruct.markers),
              
              eventarray = zeros(EEG.trials, length(eventstruct.markers));
              allfields  = { eventstruct.markers.marker_names };
              
              for index = 1:length(eventstruct.markers)
                  if ~isempty( eventstruct.markers(index).trial_times )
                      indval  = eventstruct.markers(index).trial_times(:,1);
                      values  = eventstruct.markers(index).trial_times(:,2);
                      if length(indval) < EEG.trials & length(unique(values)) == 1
                          
                          % non latency field
                          % -----------------
                          otherfields{end+1} =  eventstruct.markers(index).marker_names;
                          if unique(values) == 0
                              eventarray(:,index)      = 1;
                              eventarray(indval,index) = 0;
                          else
                              eventarray(:,index)      = 0;
                              eventarray(indval,index) = values;
                          end;
                      else
                          timefields{end+1} =  eventstruct.markers(index).marker_names;
                          eventarray(indval,index) = values;
                      end;
                  end;
              end;
              
          end
      end
      
      EEG.eventdescription = {};
      if ~isempty(eventarray),
          if isempty(alltrials)
              EEG = pop_importepoch(EEG, eventarray, allfields, 'latencyfields', timefields, 'timeunit', 1);
          else
              EEG = pop_importepoch(EEG, eventarray(alltrials,:), allfields, 'latencyfields', timefields, 'timeunit', 1);
          end;
      end;
  %catch
  %    disp(lasterr);
  %    disp('error (see above) while importing events: events not imported');
  %end;
  
  % command
  cd(orifold);
  com = sprintf('EEG = pop_ctf_read(''%s'', %s)', folder, vararg2str(options));
    
% format folder
% -------------
function [folder,parentfolder] = formatfolder( orifolder )
    delims = [ find(orifolder == '/') find(orifolder == '\') ];
    if ~isempty(delims)
        if delims(end) == length(orifolder)
            folder = orifolder(delims(end-1)+1:end-1);
            parentfolder = orifolder(1:delims(end-1)-1);
        else
            folder = orifolder(delims(end)+1:end-1);
            parentfolder = orifolder(1:delims(end)-1);
        end;
    else
        folder = orifolder;
        parentfolder = '.';
    end;      
