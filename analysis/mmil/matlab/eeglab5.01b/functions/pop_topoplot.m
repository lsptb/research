% pop_topoplot() - Plot scalp map(s) in a figure window. If number of input
%                  arguments is less than 3, pop up an interactive query window.
%                  Makes (possibly repeated) calls to topoplot().
% Usage:
%   >> pop_topoplot( EEG); % pops up a parameter query window
%   >> pop_topoplot( EEG, typeplot, items, title, plotdip, options...); % no pop-up
%
% Inputs:
%   EEG      - Input EEG dataset structure (see >> help eeglab)
%   typeplot - 1-> Plot channel ERP maps; 0-> Plot component maps {default:1}.
%
% Commandline inputs also set in pop-up window:
%   items    - [array] If typeplot==1 (ERP maps), within-epoch latencies 
%              (ms) at which to plot the maps. If typeplot==0 (component
%              maps), component indices to plot. In this case,
%              negative map indices -> invert map polarity, or 
%              NaN -> leave a blank subplot. (Ex: [1 -3 NaN 4])
%   title    - Plot title.
%   rowscols - Vector of the form [m,n] giving [rows, cols] per page.
%              If the number of maps exceeds m*n, multiple figures 
%              are produced {default|0 -> one near-square page}.
%   plotdip  - [0|1] plot associated dipole(s) for scalp map if present
%              in dataset.
%   options  - optional topoplot() arguments. Separate using commas. 
%              Example 'style', 'straight'. See >> help topoplot
%              for further details. {default: none}
%
% Note:
%   A new figure is created automatically only when the pop_up window is 
%   called or when more than one page of maps are plotted. Thus, this 
%   command may be used to draw topographic maps in a figure sub-axis.
%
% Author: Arnaud Delorme, CNL / Salk Institute, 2001
%
% See also: topoplot(), eeglab()

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

% $Log: pop_topoplot.m,v $
% Revision 1.65  2005/09/05 21:11:43  scott
% made 'electrodes' 'off' the default for > 64 channels -sm
%
% Revision 1.64  2005/03/31 02:52:39  arno
% dipole plot
%
% Revision 1.63  2005/03/18 18:02:24  arno
% fix plot component dipole for spherical model
%
% Revision 1.62  2005/03/10 18:46:30  arno
% nothing
%
% Revision 1.61  2005/03/10 17:48:40  arno
% converting MNI coordinates back to spherical
%
% Revision 1.60  2005/03/05 00:05:52  arno
% same
%
% Revision 1.59  2005/03/05 00:04:42  arno
% adding chaninfo
%
% Revision 1.58  2005/02/07 18:01:47  scott
% adjust pop window text
%
% Revision 1.57  2004/11/09 19:22:47  arno
% testing if curpos empty
%
% Revision 1.56  2004/09/15 19:05:14  hilit
% undo previous change
%
% Revision 1.55  2004/09/15 18:58:41  hilit
%
% Revision 1.54  2004/09/13 19:31:00  arno
% reverting version 1.49
%
% Revision 1.49  2004/09/13 18:21:02  arno
% debug ploting for Matlab 7
%
% Revision 1.48  2004/09/13 16:29:17  arno
% remove debug message
%
% Revision 1.47  2004/09/10 21:07:23  hilit
% more changes
%
% Revision 1.46  2004/09/10 21:02:15  hilit
% corrected a problem with the current axes
%
% Revision 1.45  2004/09/10 17:35:07  hilit
% added the current axes setting to try and fix the plotting bug in Matlab 7
%
% Revision 1.44  2004/08/31 19:32:16  scott
% still trying - but plot 2 is drawn on plot 1 (Matlab7) ?? -sm
%
% Revision 1.43  2004/08/31 18:59:49  scott
% debug same
%
% Revision 1.42  2004/08/31 18:57:42  scott
% added figure(curfig) calls to keep plotting straight in Matlab 7.0.0 -sm
%
% Revision 1.41  2004/03/18 00:32:19  arno
% remove skirt
%
% Revision 1.40  2004/02/22 23:54:00  scott
% making electrodes,on the default option
%
% Revision 1.39  2004/02/22 23:48:49  scott
% made srhink,skirt the default option
%
% Revision 1.38  2004/02/17 22:48:39  arno
% add dipsphere
%
% Revision 1.37  2004/01/30 15:59:23  arno
% topoplot msg displayed only for 1st plot
%
% Revision 1.36  2003/12/11 16:16:42  arno
% debug history generation
%
% Revision 1.35  2003/11/06 02:12:57  arno
% put RV to title
%
% Revision 1.34  2003/11/06 01:54:14  arno
% adding dipole plot
%
% Revision 1.33  2003/11/05 20:35:54  arno
% plot dipoles
%
% Revision 1.32  2003/10/29 01:06:08  arno
% remove debuging message
%
% Revision 1.31  2003/10/11 02:26:01  arno
% modify default title
%
% Revision 1.30  2003/10/11 02:20:31  arno
% shift colorbar to left for 1 plot
%
% Revision 1.29  2003/10/11 02:11:42  arno
% debuging for empty option
%
% Revision 1.28  2003/08/07 20:49:15  arno
% option 'masksurf' to speed up display
%
% Revision 1.27  2003/07/17 23:20:16  scott
% formatting
%
% Revision 1.26  2003/07/03 16:55:07  arno
% header
%
% Revision 1.25  2003/07/03 16:52:44  arno
% same
%
% Revision 1.24  2003/07/03 16:51:16  arno
% header text
%
% Revision 1.23  2003/05/12 22:27:56  arno
% verbose option off
%
% Revision 1.22  2003/03/12 03:20:41  arno
% help button
%
% Revision 1.21  2003/03/06 02:20:39  arno
% use vararg2str
%
% Revision 1.20  2002/10/26 20:26:50  arno
% help msg
%
% Revision 1.19  2002/10/11 22:15:59  arno
% header
%
% Revision 1.18  2002/08/27 00:38:21  arno
% cl more optimal auto location
%
% Revision 1.17  2002/08/22 17:38:39  arno
% correct title
%
% Revision 1.16  2002/08/22 17:01:14  arno
% debugging for 1 plot
%
% Revision 1.15  2002/08/20 00:05:38  arno
% adding test for plotting a large number of components
%
% Revision 1.14  2002/08/19 22:20:21  arno
% 10 -> strvcat
%
% Revision 1.13  2002/08/19 22:15:02  arno
% text for components
%
% Revision 1.12  2002/08/17 22:30:28  scott
% ERP/components
%
% Revision 1.11  2002/08/17 22:26:05  scott
% typo
%
% Revision 1.10  2002/08/17 22:24:15  scott
% tyop
%
% Revision 1.9  2002/08/17 22:21:02  scott
% help msg and menu text
%
% Revision 1.8  2002/08/17 22:16:19  scott
% editing help msg
%
% Revision 1.7  2002/08/17 22:15:51  scott
% *** empty log message ***
%
% Revision 1.6  2002/08/13 17:49:37  arno
% debug color
%
% Revision 1.5  2002/08/12 16:31:20  arno
% inputdlg2
%
% Revision 1.4  2002/08/12 01:46:44  arno
% color
%
% Revision 1.3  2002/08/11 22:22:19  arno
% color
%
% Revision 1.2  2002/04/18 18:26:14  arno
% typo can not
%
% Revision 1.1  2002/04/05 17:32:13  jorn
% Initial revision
%

% 01-25-02 reformated help & license -ad 
% 02-15-02 text interface editing -sm & ad 
% 02-16-02 added axcopy -ad & sm
% 03-18-02 added title -ad & sm

function com = pop_topoplot( EEG, typeplot, arg2, topotitle, rowcols, varargin);

com = '';
if nargin < 1
   help pop_topoplot;
   return;
end;   
if nargin < 2   
   typeplot = 1;
end;
if typeplot == 0 & isempty(EEG.icasphere)
   disp('Error: no ICA data for this set, first run ICA'); return;
end;   
if isempty(EEG.chanlocs)
   disp('Error: cannot plot topography without channel location file'); return;
end;   

if nargin < 3
	% which set to save
	% -----------------
	if typeplot
		txtwhat2plot1 = 'Plotting ERP scalp maps at these latencies';
        txtwhat2plot2 = sprintf('(range: %d to %d ms, NaN -> empty):', ...
                                round(EEG.xmin*1000), round(EEG.xmax*1000));
        editwhat2plot = [''];
	else
		txtwhat2plot1 = 'Component numbers';
		txtwhat2plot2 = '(negate index to invert component polarity; NaN -> empty subplot; Ex: -1 NaN 3)';
        editwhat2plot = ['1:' int2str(size(EEG.icaweights,1))];
 	end;	
        if EEG.nbchan > 64, 
            elecdef = ['''electrodes'', ''on''']; 
        else, 
            elecdef = ['''electrodes'', ''on''']; 
        end;
    uilist = { { 'style'   'text'     'string'    txtwhat2plot1 } ...
               { 'style'   'edit'     'string'    editwhat2plot } ...
               { 'style'   'text'     'string'    txtwhat2plot2 } ...
               { } ...
               { 'style'   'text'     'string'    'Plot title' } ...
               { 'style'   'edit'     'string'    fastif(~isempty(EEG.setname), [EEG.setname], '') } ...
               { 'style'   'text'     'string'    'Plot geometry (rows,col.); [] -> near square' } ...
               { 'style'   'edit'     'string'    '[]' } ...
               { 'style'   'text'     'string'    'Plot associated dipole(s) (if present)' } ...
               { 'style'   'checkbox' 'string'    '' } { } ...
               { } ...
               { 'style'   'text'     'string'    [ '-> Additional topoplot()' fastif(typeplot,'',' (and dipole)') ...
                                                  ' options (see Help)' ] } ...
               { 'style'   'edit'     'string'    elecdef } };
    uigeom = { [1.5 1] [1] [1] [1.5 1] [1.5 1] [1.55 0.2 0.8] [1] [1] [1] };
    if typeplot
        uilist(9:11) = [];
        uigeom(6) = [];
    end;
    guititle = fastif( typeplot, 'Plot ERP scalp maps in 2-D -- pop_topoplot()', ...
                       'Plot component scalp maps in 2-D -- pop_topoplot()');
    
    result = inputgui( uigeom, uilist, 'pophelp(''pop_topoplot'')', guititle, [], 'normal');
	if length(result) == 0 return; end;
	
    % reading first param
    % -------------------
    arg2   	     = eval( [ '[' result{1} ']' ] );
	if length(arg2) > EEG.nbchan
		tmpbut = questdlg2(...
                  ['This involves drawing ' int2str(length(arg2)) ' plots. Continue ?'], ...
                         '', 'Cancel', 'Yes', 'Yes');
		if strcmp(tmpbut, 'Cancel'), return; end;
	end;
    if isempty(arg2), error('Nothing to plot; enter parameter in first edit box'); end;
    
    % reading other params
    % --------------------
	topotitle   = result{2};
	rowcols     = eval( [ '[' result{3} ']' ] );
    if typeplot
        plotdip = 0;
        try, options      = eval( [ '{ ' result{4} ' }' ]);
        catch, error('Invalid scalp map options'); end;
    else
        plotdip     = result{4};
        try, options      = eval( [ '{ ' result{5} ' }' ]);
        catch, error('Invalid scalp map options'); end;
    end;        
    if length(arg2) == 1, 
      figure; curfig=gcf; 
      try, icadefs; 
         set(curfig, 'color', BACKCOLOR); 
      catch, end; 
    end;
else
    if ~isempty(varargin) & isnumeric(varargin{1})
        plotdip = varargin{1};
        varargin = varargin(2:end);
    else
        plotdip = 0;
    end;
    options = varargin;
end;

% additional options
% ------------------
options    = { options{:} 'masksurf' 'on' };
outoptions = { options{:} }; % for command

nbgraph = size(arg2(:),1);
if ~exist('topotitle')
    topotitle = '';
end;    
if ~exist('rowcols') | isempty(rowcols) | rowcols == 0
    rowcols(2) = ceil(sqrt(nbgraph));
    rowcols(1) = ceil(nbgraph/rowcols(2));
end;    

SIZEBOX = 150;

fprintf('Plotting...\n');
if isempty(EEG.chanlocs)
	fprintf('Error: set has no channel location file\n');
	return;
end;

% determine the scale for plot of different times (same scales)
% -------------------------------------------------------------
if typeplot
	SIGTMP = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
	pos = round( (arg2/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
	indexnan = find(isnan(pos));
	nanpos = find(isnan(pos));
	pos(nanpos) = 1;
	SIGTMPAVG = mean(SIGTMP(:,pos,:),3);
	SIGTMPAVG(:, nanpos) = NaN;
	maxlim = max(SIGTMPAVG(:));
	minlim = min(SIGTMPAVG(:));
	maplimits = [ -max(maxlim, -minlim) max(maxlim, -minlim)];
end;
	
% plot the graphs
% ---------------
counter = 1;
countobj = 1;
allobj = zeros(1,1000);
curfig = gcf;
for index = 1:size(arg2(:),1)
	if nbgraph > 1
        if mod(index, rowcols(1)*rowcols(2)) == 1
            if index> 1, figure(curfig); a = textsc(0.5, 0.05, topotitle); set(a, 'fontweight', 'bold'); end;
        	curfig = figure;
			pos = get(curfig,'Position');
			posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
			posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
			set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
			try, icadefs; set(curfig, 'color', BACKCOLOR); catch, end;
       end;    
		curax = subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1);
        set(curax, 'visible', 'off')
    end;

	% add dipole location if present
    % ------------------------------
    options = outoptions;
    dipoleplotted = 0;
    if plotdip & typeplot == 0
        if isfield(EEG, 'dipfit') & isfield(EEG.dipfit, 'model')
            if length(EEG.dipfit.model) >= index
                %curpos = EEG.dipfit.model(arg2(index)).posxyz/EEG.dipfit.vol.r(end);
                curpos = EEG.dipfit.model(arg2(index)).posxyz;
                curmom = EEG.dipfit.model(arg2(index)).momxyz;
                try, 
                      select = EEG.dipfit.model(arg2(index)).select;
                catch select = 0;
                end;
                if strcmpi(EEG.dipfit.coordformat, 'MNI') % from MNI to sperical coordinates
                    transform = pinv( sph2spm );                    
                    tmpres = transform * [ curpos(1,:) 1 ]'; curpos(1,:) = tmpres(1:3);
                    tmpres = transform * [ curmom(1,:) 1 ]'; curmom(1,:) = tmpres(1:3);
                    try, tmpres = transform * [ curpos(2,:) 1 ]'; curpos(2,:) = tmpres(1:3); catch, end;
                    try, tmpres = transform * [ curmom(2,:) 1 ]'; curmom(2,:) = tmpres(1:3); catch, end;
                end;
                curpos = curpos / 85;
                if ~isempty(curpos)
                    if size(curpos,1) > 1 & length(select) == 2
                        options = { options{:} 'dipole' [ curpos(:,1:2) curmom(:,1:3) ] };
                        dipoleplotted = 1;
                    else
                        if any(curpos(1,:) ~= 0)
                            options = { options{:} 'dipole' [ curpos(1,1:2) curmom(1,1:3) ] };
                            dipoleplotted = 1;
                        end;
                    end;
                end;
                if nbgraph ~= 1
                    options = {  options{:} 'dipscale' 0.6 };
                end;
                %options = { options{:} 'dipsphere' max(EEG.dipfit.vol.r) };
            end;
        end;
    end;
    
	% plot scalp map
    % --------------
    if index == 1
        addopt = { 'verbose', 'on' };
    else 
        addopt = { 'verbose', 'off' };
    end;
    %fprintf('Printing to figure %d.\n',curfig);
    if isfield(EEG, 'chaninfo'), options = { options{:} 'chaninfo' EEG.chaninfo }; end;
    if ~isnan(arg2(index))
		if typeplot
            if nbgraph > 1, axes(curax); end;
            tmpobj = topoplot( SIGTMPAVG(:,index), EEG.chanlocs, 'maplimits', maplimits, addopt{:}, options{:});
			if nbgraph == 1, 
                 figure(curfig); if nbgraph > 1, axes(curax); end;
                 title( [ 'Latency ' int2str(arg2(index)) ' ms from ' topotitle] );
			else 
                 figure(curfig); if nbgraph > 1, axes(curax); end; 
                 title([int2str(arg2(index)) ' ms']);
			end;
		else
            if arg2(index) < 0
                 figure(curfig);  if nbgraph > 1, axes(curax); end;
                 tmpobj = topoplot( -EEG.icawinv(:, -arg2(index)), EEG.chanlocs, addopt{:}, options{:} );
            else
                 figure(curfig);  if nbgraph > 1, axes(curax); end;
                 tmpobj = topoplot( EEG.icawinv(:, arg2(index)), EEG.chanlocs, addopt{:}, options{:} );
            end;    			
			if nbgraph == 1, texttitle = [ 'IC ' int2str(arg2(index)) ' from ' topotitle];
			else             texttitle = ['' int2str(arg2(index))];
			end;
            if dipoleplotted, texttitle = [ texttitle ' (' num2str(EEG.dipfit.model(arg2(index)).rv*100,2) '%)']; end;
            figure(curfig);  if nbgraph > 1, axes(curax); end; title(texttitle);
		end;
        allobj(countobj:countobj+length(tmpobj)-1) = tmpobj;
        countobj = countobj+length(tmpobj);
		drawnow;
		axis square; 
		if index == size(arg2(:),1)
	        if nbgraph == 1
                clim = get(gca, 'clim');
                pos = get(gca,'position');
                q = [pos(1) pos(2) 0 0];
                s = [pos(3) pos(4) pos(3) pos(4)];
                col = colormap;
                ax = axes('position', [0.95 0 .05 1].*s+q);
                cbar(ax,[1:64],clim);
	        else 
                cbar('vert');
            end;
	    end;	   
    else
    axis off
    end;
end;
if nbgraph> 1, 
   figure(curfig); a = textsc(0.5, 0.05, topotitle); 
   set(a, 'fontweight', 'bold'); 
end;
if nbgraph== 1, 
   com = 'figure;'; 
end;
set(allobj(1:countobj-1), 'visible', 'on');

figure(curfig);
axcopy(curfig, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); clear postmp;');

com = [com sprintf('pop_topoplot(%s,%d, %s);', ...
                   inputname(1), typeplot, vararg2str({arg2 topotitle rowcols plotdip outoptions{:} }))];
return;

		
