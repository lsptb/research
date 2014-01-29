% pop_writelocs() - load a EGI EEG file (pop out window if no arguments).
%
% Usage:
%   >> EEG = pop_writelocs(chanstruct);             % a window pops up
%   >> EEG = pop_writelocs(chanstruct, filename, 'key', val, ...);
%
% Inputs:
%   chanstruct     - channel structure. See readlocs()
%   filename       - Electrode location file name
%   'key',val      - same as writelocs()
% 
% Author: Arnaud Delorme, CNL / Salk Institute, 17 Dec 2002
%
% See also: writelocs()

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 17 Dec 2002 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: pop_writelocs.m,v $
% Revision 1.12  2005/09/27 22:09:03  arno
% fixing writing locations (removing .ced)
%
% Revision 1.11  2005/03/09 18:50:12  arno
% nothing
%
% Revision 1.10  2005/03/09 18:47:47  arno
% warning message
%
% Revision 1.9  2005/03/09 18:45:57  arno
% make compatible to new format returned by readlocs
%
% Revision 1.8  2003/05/13 23:24:20  arno
% nothing
%
% Revision 1.7  2003/05/13 23:21:27  arno
% debug unicoord
%
% Revision 1.6  2003/05/13 23:09:04  arno
% same
%
% Revision 1.5  2003/05/13 22:32:10  arno
% debuging for pop_chanedit
%
% Revision 1.4  2003/05/13 22:27:10  arno
% debug for pop_chanedit
%
% Revision 1.3  2003/05/13 22:10:43  arno
% debuging
%
% Revision 1.2  2003/05/13 21:10:47  arno
% writing only allow for a subset of file formats
%
% Revision 1.1  2002/12/24 01:25:18  arno
% Initial revision
%
% Revision 1.2  2002/11/14 23:35:36  arno
% header
%
% Revision 1.1  2002/11/13 02:34:22  arno
% Initial revision
%

function com = pop_writelocs(chans, filename, varargin); 
    
com = '';
if nargin < 1
   help pop_writelocs;
   return;
end;

if isfield(chans, 'shrink')
    chans = rmfield(chans, 'shrink');
    disp('Warning: shrink factor ignored');
end;

disp('WARNING: ELECTRODE COORDINATES MUST BE WITH NOSE ALONG THE +X DIMENSION TO BE EXPORTED')
disp('         IF NOT, THE EXPORTED FILE COORDINATES MAY BE INNACURATE')

% get infos from readlocs
% -----------------------
[chanformat listcolformat] = readlocs('getinfos');
chanformat(end)    = [];
listcolformat(end) = []; % remove chanedit
chanformat(end)    = [];
listcolformat(end) = []; % remove chanedit
indformat  = [];
for index = 1:length(chanformat), 
    if ~isstr(chanformat(index).importformat)
        indformat = [ indformat index ];
    end;
    if isempty(chanformat(index).skipline), chanformat(index).skipline = 0; end;
end;
listtype   = { chanformat(indformat).type };
formatinfo = { chanformat(indformat).importformat };
formatskip = [ chanformat(indformat).skipline ];
   
%[listtype formatinfo listcolformat formatskip] = readlocs('getinfoswrite');

listtype{end+1} = 'custom';
formatinfo{end+1} = {};
formatskip = [ formatskip 0];

if nargin < 2
   updatefields = [ 'tmpdata = get(gcf, ''userdata'');' ...
                    'tmpobj = findobj(gcf, ''tag'', ''list2'');' ...
                    'set(tmpobj, ''string'', strvcat(tmpdata{2}));' ...
                    'clear tmpobj tmpdata;' ];
   addfieldcom = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                   'tmpobj = findobj(gcf, ''tag'', ''list1'');' ...
                   'tmpdata{2}{end+1} = tmpdata{1}{get(tmpobj, ''value'')};' ...
                   'set(gcbf, ''userdata'', tmpdata);' ...
                   updatefields ];
   rmfieldcom  = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                   'tmpobj = findobj(gcbf, ''tag'', ''list2'');' ...
                   'try, tmpdata{2}(get(tmpobj, ''value'')) = [];' ...
                   '    set(tmpobj, ''value'', 1);' ...
                   'catch, end;' ...
                   'set(gcbf, ''userdata'', tmpdata);' ...
                   updatefields ];                      
   filetypecom = [ 'tmpdata = get(gcf, ''userdata'');' ...
                   'tmpobj = findobj(gcf, ''tag'', ''formatlist'');' ...
                   'tmpval = get(tmpobj, ''value'');' ...
                   'try, tmpdata{2} = tmpdata{3}{tmpval}; catch, end;' ... %try and catch for custom
                   'set(gcf, ''userdata'', tmpdata);' ...
                   updatefields ...
                   'tmpdata = get(gcf, ''userdata'');' ...
						 'tmpobj1 = findobj(gcf, ''tag'', ''insertcol'');' ... % the lines below
                   'tmpobj2 = findobj(gcf, ''tag'', ''inserttext'');' ... % update the checkbox
                   'try, ' ...                                     % and the edit text box
                   '  if tmpdata{4}(tmpval) == 2,' ...
                   '     set(tmpobj1, ''value'', 1);' ...
                   '  else,' ...
                   '     set(tmpobj1, ''value'', 0);' ...
                   '  end;' ...
                   '  if tmpval == 1,' ... % besa only
                   '     set(tmpobj2, ''string'', ''' int2str(length(chans)) ''');' ...
                   '  else,' ...
                   '     set(tmpobj2, ''string'', '''');' ...
                   '  end;' ...
                   'catch, end;' ... % catch for custom case
                   'tmpobj = findobj(gcf, ''userdata'', ''setfield'');' ...
                   'if tmpval == ' int2str(length(listtype)) ',' ... % disable if non-custom type
                   '   set(tmpobj, ''enable'', ''on'');' ...
                   'else,' ...
                   '   set(tmpobj, ''enable'', ''off'');' ...
                   'end; clear tmpobj tmpobj2 tmpdata tmpval;' ];
                
   geometry = { [1 1 1] [1 1] [1] [1 1 1] [1 1] [1 1 1] [1 0.3 0.7] [1] [1] };
   listui = { ...
         { 'style' 'text' 'string' 'Filename' } ...
         { 'style' 'edit' 'string' '' 'tag' 'filename' 'horizontalalignment' 'left' } ...
         { 'style' 'pushbutton' 'string' 'Browse' 'callback' ...
           [  '[tmpfile tmppath] = uiputfile(''*'', ''Exporting electrode location file -- pop_writelocs()'');' ... 
              'set(findobj(gcbf, ''tag'', ''filename''), ''string'', char([tmppath tmpfile ]));' ...
              'clear tmpfile tmppath;' ] } ...
         { 'style' 'text' 'string' strvcat('Select output file type', ' ', ' ') } ...
         { 'style' 'listbox' 'tag' 'formatlist' 'string' strvcat(listtype) ...
            'value' length(listtype) 'callback' filetypecom } ...
         { 'style' 'text' 'string' 'Select fields to export below' } ...
         { } { 'style' 'pushbutton' 'string' '-> ADD' 'callback' addfieldcom 'userdata' 'setfield' } { } ...
         { 'style' 'listbox' 'tag' 'list1' 'string' strvcat(fieldnames(chans)) 'userdata' 'setfield' } ...
         { 'style' 'listbox' 'tag' 'list2' 'string' '' 'userdata' 'setfield2' } ...
         { } { 'style' 'pushbutton' 'string' 'REMOVE <-' 'callback' rmfieldcom 'userdata' 'setfield' } { } ...
         { 'style' 'text' 'string' 'Insert column names' } ...
         { 'style' 'checkbox' 'tag' 'insertcol' 'value' 1 'userdata' 'setfield' } { } ...
         { 'style' 'text' 'string' 'Enter custom header below' } ...
         { 'style' 'edit' 'userdata' 'setfield' 'tag' 'inserttext' 'horizontalalignment' 'left' 'max' 2 } ...
	};
   
   inputgui(geometry, listui, 'pophelp(''writelocs'');', ...
      'Exporting electrode location file -- pop_writelocs()', { fieldnames(chans) {} formatinfo formatskip }, 'plot', [1 3 1 1 3 1 1 1 3 ]);
   fig = gcf;
   
   % set default format
   tmpobj = findobj(fig, 'tag', 'formatlist'); 
   set(tmpobj, 'value', 6);
   eval(get(tmpobj, 'callback')); 
   
   res = inputgui(geometry, listui, 'pophelp(''writelocs'');', ...
      'Exporting electrode location file -- pop_writelocs()', { listcolformat {} formatinfo formatskip }, fig, [1 3 1 1 3 1 1 1 3 ]);
   if gcf ~= fig, return; end;
   exportfields = get(fig, 'userdata');
   exportfields = exportfields{2};
   close(fig);
   
   % decode the inputs
   filename = res{1};
   if isempty(filename), 
      errordlg2('Error: Empty file name', 'Error');
      return;
   end;
   options = { 'filetype' listtype{res{2}} 'format' exportfields ...
         'header' fastif(res{5}, 'on', 'off') 'customheader' res{6} };
else
	options = varargin;   
end;

% generate history
% ----------------
if isempty(inputname(1)) % not a variable name -> probably the structure from pop_chanedit
    writelocs(chans, filename, options{:});
   com = sprintf('pop_writelocs( EEG.chanlocs, ''%s'', %s);', filename, vararg2str(options));
else
    if strcmpi(inputname(1), 'chantmp')
        % do not write file (yet)
        com = sprintf('pop_writelocs( chans, ''%s'', %s);', filename, vararg2str(options));
    else
        writelocs(chans, filename, options{:});
        com = sprintf('pop_writelocs( %s, ''%s'', %s);', inputname(1), filename, vararg2str(options));
    end;
end;
