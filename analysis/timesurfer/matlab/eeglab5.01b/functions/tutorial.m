% tutorial() - Bring up the ICA / electrophysiology toolbox tutorial
%              in a browser window (see docopt.m in the toolbox dir).
%              Tutorial URL: http://www.sccn.ucsd.edu/tutorial/
%              Download: See http://www.sccn.ucsd.edu/ica.html
%
% Authors: Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD, 12/29/00

% Copyright (C) 12/29/00 Scott Makeig & Tzyy-Ping Jung, SCCN/INC/UCSD,
% scott@sccn.ucsd.edu
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

% $Log: tutorial.m,v $
% Revision 1.4  2002/11/15 15:55:10  arno
% debugging for windows
%
% Revision 1.3  2002/11/15 15:48:21  arno
% simplifying
%
% Revision 1.2  2002/08/14 17:03:04  arno
% updating tutdir
%
% Revision 1.1  2002/04/05 17:36:45  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 

icadefs % load icadefs.m globals including TUTDIR

%TUTDIR = which('eeglab');
%TUTDIR = TUTDIR(1:findstr(TUTDIR, 'eeglab')-1);

%if exist([TUTDIR 'index.html'])
%   eval(['web file://' TUTDIR 'index.html' ]);
%else
%   fprintf('ICA Matlab Toolbox Tutorial not found in the toolbox directory.\n');
   fprintf('Opening the toolbox www site ...\n\n');
   web(TUTORIAL_URL);
%end

