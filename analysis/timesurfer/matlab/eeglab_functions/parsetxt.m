% parsetxt() - parse text input into cell array
%
% Usage: >> cellarray = parsetxt( txt, delims );
%
% Inputs: 
%    txt    - input text
%    delims - optional char array of delimiters (default: [' ' ',' 9]);
%
% Note: commas, and simple quotes are ignored
%
% Author: Arnaud Delorme, CNL / Salk Institute, 18 April 2002

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) 18 April 2002 Arnaud Delorme, SCCN/INC/UCSD, scott@sccn.ucsd.edu
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

% $Log: parsetxt.m,v $
% Revision 1.2  2006/03/02 23:45:37  arno
% nothing
%
% Revision 1.1  2002/04/05 17:39:45  jorn
% Initial revision
%
% 04-03-02 reprogrammed the function with findstr -ad

function cellarray = parsetxt(txt, delims);

if nargin < 1
	help parsetxt;
	return;
end;	

if nargin < 2
    delims = [' ' ',' 9 '"' '''' ];
end;
    
cellarray = {};
tmptxt = '';
for index =1:length(txt)
    if ~isempty(findstr(txt(index), delims))
        if ~isempty(tmptxt), cellarray = { cellarray{:}, tmptxt }; end;
        tmptxt = '';
    else
        tmptxt = [ tmptxt txt(index) ];
    end;
end;        
if ~isempty(tmptxt)
    cellarray = { cellarray{:}, tmptxt };
end;
return;
