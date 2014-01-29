% pop_chansel() - pop up a graphic interface to select channels
%
% Usage:
%   >> [chanlist] = pop_chansel(chanstruct); % a window pops up
%   >> [chanlist channames strallnames] = ...
%                        pop_chansel(chanstruct, 'key', 'val', ...);
%
% Inputs:
%   chanstruct     - channel structure. See readlocs()
%
% Optional input:
%   'withindex'      - ['on'|'off'] add index to each entry. May also a be 
%                      an array of indices
%   'select'         - selection of channel. Can take as input all the
%                      outputs of this function.
%   'selectionmode' - selection mode 'multiple' or 'single'. See listdlg2().
%
% Output:
%   chanlist  - indices of selected channels
%   channames - names of selected channels
%   strallnames - all channel names concatenated
%
% Author: Arnaud Delorme, CNL / Salk Institute, 3 March 2003

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 3 March 2003 Arnaud Delorme, Salk Institute, arno@salk.edu
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

% $Log: pop_chansel.m,v $
% Revision 1.19  2006/03/29 00:49:38  scott
% txt
%
% Revision 1.18  2006/01/10 00:42:58  arno
% fixing scrolling more than 60 channels
%
% Revision 1.17  2005/09/27 22:03:07  arno
% fix argument 'withindex' must be numeric
%
% Revision 1.16  2005/07/28 15:47:51  arno
% allow using input indices
%
% Revision 1.15  2004/11/10 17:34:46  arno
% add selection mode
%
% Revision 1.14  2004/11/10 17:27:07  arno
% debug last
%
% Revision 1.13  2004/11/10 16:48:02  arno
% NEW CHANNEL SELECTOR
%
% Revision 1.12  2004/11/10 16:09:45  arno
% nothing
%
% Revision 1.11  2003/08/05 18:20:05  arno
% same
%
% Revision 1.10  2003/08/05 18:17:49  arno
% assign default arguments
%
% Revision 1.9  2003/05/14 18:12:32  arno
% typo
%
% Revision 1.8  2003/04/16 00:25:58  arno
% also generate a string with all channels names
%
% Revision 1.7  2003/04/16 00:16:43  arno
% returning channel names
%
% Revision 1.6  2003/03/05 18:53:32  arno
% handle empty entries
%
% Revision 1.5  2003/03/05 18:46:49  arno
% debug for numerical channels
%
% Revision 1.4  2003/03/05 18:34:50  arno
% same
%
% Revision 1.3  2003/03/05 18:33:27  arno
% handling cancel
%
% Revision 1.2  2003/03/04 15:06:27  roberto
% no change
%
% Revision 1.1  2003/03/03 19:32:31  arno
% Initial revision
%

function [chanlist,chanliststr, allchanstr] = pop_chansel(chans, varargin); 
    
    if nargin < 1
        help pop_chansel;
        return;
    end;
    if isempty(chans), disp('Empty input'); return; end;
    if ~iscell(chans), error('Can only process cell array'); end;
    chanlist    = [];
    chanliststr = {};
    allchanstr  = '';
    
    g = finputcheck(varargin, { 'withindex'     {  'integer' 'string' } { [] {'on' 'off'} }   'off';
                                'select'        { 'cell' 'string' 'integer' } [] [];
                                'selectionmode' 'string' { 'single' 'multiple' } 'multiple'});
    if isstr(g), error(g); end;
    if ~isstr(g.withindex), chan_indices = g.withindex; g.withindex = 'on';
    else                    chan_indices = 1:length(chans);
    end;
    
    % convert selection to integer
    % ----------------------------
    if isstr(g.select) & ~isempty(g.select)
        g.select = parsetxt(g.select);
    end;
    if iscell(g.select) & ~isempty(g.select)
        if isstr(g.select{1})
            tmplower = lower( chans );
            for index = 1:length(g.select)
                matchind = strmatch(lower(g.select{index}), tmplower, 'exact');
                if ~isempty(matchind), g.select{index} = matchind;
                else error( [ 'Cannot find ''' g.select{index} '''' ] );
                end;
            end;
        end;
        g.select = [ g.select{:} ];
    end;
    if ~isnumeric( g.select ), g.select = []; end;
    
    % add index to channel name
    % -------------------------
	tmpstr = {chans};
    if isnumeric(chans{1})
        tmpstr = [ chans{:} ];
        tmpfieldnames = cell(1, length(tmpstr));
        for index=1:length(tmpstr), 
            if strcmpi(g.withindex, 'on')
                tmpfieldnames{index} = [ num2str(chan_indices(index)) '  -  ' num2str(tmpstr(index)) ]; 
            else
                tmpfieldnames{index} = num2str(tmpstr(index)); 
            end;
        end;
    else
        tmpfieldnames = chans;
        if strcmpi(g.withindex, 'on')
            for index=1:length(tmpfieldnames), 
                tmpfieldnames{index} = [ num2str(chan_indices(index)) '  -  ' tmpfieldnames{index} ]; 
            end;
        end;
    end;
    [chanlist,tmp,chanliststr] = listdlg2('PromptString',strvcat('(use shift|Ctrl to', 'select several)'), ...
                'ListString', tmpfieldnames, 'initialvalue', g.select, 'selectionmode', g.selectionmode);       
    allchanstr = chans(chanlist);
    
    % get concatenated string (if index)
    % -----------------------
    if strcmpi(g.withindex, 'on')
        if isnumeric(chans{1})
            chanliststr = num2str(allchanstr);
        else
            chanliststr = '';
            for index = 1:length(allchanstr)
                chanliststr = [ chanliststr allchanstr{index} ' ' ];
            end;
            chanliststr = chanliststr(1:end-1);
        end;
    end;
       
    return;
    % old version
    % -----------
    updatefields = [ 'tmpdata = get(gcf, ''userdata'');' ...
                     'tmpobj = findobj(gcf, ''tag'', ''list2'');' ...
                     'set(tmpobj, ''string'', strvcat(tmpdata{2}));' ...
                     'clear tmpobj tmpdata;' ];
    addfieldcom = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                    'tmpobj = findobj(gcf, ''tag'', ''list1'');' ...
                    'if strmatch(  tmpdata{1}{get(tmpobj, ''value'')}, tmpdata{2} ),' ...
                    '   clear tmpobj tmpdata; return;' ...
                    'end;' ...
                    'tmpdata{2}{end+1} = tmpdata{1}{get(tmpobj, ''value'')};' ...
                    'set(gcbf, ''userdata'', tmpdata);' ...
                    updatefields ];
    rmfieldcom  = [ 'tmpdata = get(gcbf, ''userdata'');' ...
                    'tmpobj = findobj(gcbf, ''tag'', ''list2'');' ...
                    'if get(tmpobj, ''value'') == length(tmpdata{2}) &' ...
                    '   length(tmpdata{2}) > 1,' ...
                    '   set(tmpobj, ''value'', get(tmpobj, ''value'')-1);' ...
                    'end;' ...
                   'try, tmpdata{2}(get(tmpobj, ''value'')) = [];' ...
                    'catch, end;' ...
                    'set(gcbf, ''userdata'', tmpdata);' ...
                    updatefields ];                      

   channelnames = strvcat({chans.labels});
   geometry = { [1] [1 4 1] [1 1] [1 4 1] };
   listui = { ...
         { 'style' 'text' 'string' 'Select channels to remove' } ...
         { } { 'style' 'pushbutton' 'string' '-> ADD TO LIST' 'callback' addfieldcom 'userdata' 'setfield' } { } ...
         { 'style' 'listbox' 'tag' 'list1' 'string' channelnames 'userdata' 'setfield' } ...
         { 'style' 'listbox' 'tag' 'list2' 'string' '' 'userdata' 'setfield' } ...
         { } { 'style' 'pushbutton' 'string' 'REMOVE FROM LIST <-' 'callback' rmfieldcom 'userdata' 'setfield' } { } ...
	};
   
   [outparams userdat] = inputgui(geometry, listui, 'pophelp(''pop_chansel'');', ...
                                  'Select channels -- pop_chansel()', { {chans.labels} {} }, 'normal', [ 1 1 10 1 1]);
   
   % output
   % ------
   if isempty(userdat), chanlist = []; return; end;
   chanliststr = userdat{2};
   chanlist = [];
   for index = 1:length(chanliststr)
       i = strmatch (chanliststr{index}, channelnames, 'exact');
       chanlist  = [chanlist i];
   end;
   [chanlist indices] = sort(chanlist);
   chanliststr = chanliststr(indices);
   
   % generate all channel name string
   if ~isempty(chanliststr)
       allchanstr = chanliststr{1};
       for index = 2:length(chanliststr)
           allchanstr = [ allchanstr ' ' chanliststr{index} ];
       end;
   end;
