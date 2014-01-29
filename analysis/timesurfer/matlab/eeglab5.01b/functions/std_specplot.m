% std_specplot() - visualizes component cluster spectra, either mean spectra for 
%                  all requested clusters in the same figure, with spectra for 
%                  different conditions (if any) plotted in different colors, 
%                  or spectra for each specified cluster in a separate figure 
%                  for each condition,  showing the cluster component spectra plus 
%                  the mean cluster spectrum (in bold). The spectra can be 
%                  plotted only if component spectra have been computed and 
%                  saved with the EEG datasets in Matlab files "[datasetname].icaspec" 
%                  using pop_preclust() or std_preclust(). Called by pop_clustedit(). 
%                  Calls std_readspec() and internal function std_plotcompspec()
% Usage:    
%              >> [STUDY] = std_specplot(STUDY, ALLEEG, key1, val1, key2, val2, ...);  
% Inputs:
%   STUDY      - STUDY structure comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - vector of EEG dataset structures for the dataset(s) in the STUDY, 
%                typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [int vector] -> cluster numbers to plot.
%                       'all' -> plot all clusters in STUDY.
%                {default: 'all'}.
%   'comps'    - [int vector] -> indices of cluster components to plot.
%                       'all' -> plot all the components in the cluster 
%                {default: 'all'}.
%   'mode'     - ['centroid'|'comps'] plotting mode. In 'centroid' mode, the average 
%                spectra of the requested clusters are plotted in the same figure, 
%                with spectra for  different conditions (if any) plotted in different 
%                colors. In 'comps' mode, spectra for each specified cluster are 
%                plotted in separate figures (per condition), each containing the
%                cluster component spectra plus the mean cluster spectrum in bold.
%                {default: 'centroid'}. Note that this option is irrelevant when 
%                component indices are provided as input.
%   'figure'   - ['on'|'off'] for the 'centroid' mode option, 'on' plots in a new 
%                figure, while 'off'  plots in the current figure. {default: 'on'}
%
% Outputs:
%   STUDY      - the input STUDY set structure modified with the plotted cluster 
%                mean spectra, to allow quick replotting (unless the cluster means 
%                already exists in the STUDY).  
%
%   Example:
%            >> [STUDY] = std_specplot(STUDY,ALLEEG, 'clusters', 2, 'mode', 'comps');
%               % Plot component spectra for Cluster 2 plus the mean cluster spectrum 
%               % (in bold). 
%
%  See also  pop_clustedit(), pop_preclust()
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

% $Log: std_specplot.m,v $
% Revision 1.19  2006/03/28 15:38:31  scott
% help msg
%
% Revision 1.18  2006/03/28 14:54:26  arno
% cell array format
%
% Revision 1.17  2006/03/21 15:36:33  arno
% new .sets format
%
% Revision 1.16  2006/03/14 03:10:48  scott
% help msg
%
% Revision 1.15  2006/03/10 18:23:37  arno
% rename variables
%
% Revision 1.14  2006/03/10 18:09:28  arno
% remove wait bar
%
% Revision 1.13  2006/03/10 16:28:07  arno
% reading spectrum
%
% Revision 1.12  2006/03/10 00:21:18  arno
% removing reference to etc field
%
% Revision 1.11  2006/03/09 18:13:21  arno
% spectrum read
%
% Revision 1.10  2006/03/08 21:04:03  arno
% typo
%
% Revision 1.9  2006/03/08 21:00:40  arno
% rename func
%
% Revision 1.8  2006/03/08 20:21:26  arno
% rename func
%
% Revision 1.7  2006/03/07 18:45:19  arno
% allow plotting parent cluster
%
% Revision 1.6  2006/02/16 19:51:45  arno
% header
%
% Revision 1.5  2006/02/16 19:47:48  arno
% set xlimits and move std_plotcompspec.m inside
%
% Revision 1.4  2006/02/16 19:01:33  arno
% plotting both conditions on the same figure
%
% Revision 1.3  2006/02/16 00:01:21  arno
% fixing axis limits
%
% Revision 1.2  2006/02/15 22:59:20  arno
% adding scaling etc...
%

function STUDY = std_specplot(STUDY, ALLEEG,  varargin)

icadefs; % read EEGLAB defaults

% Set default values
cls = []; % plot all clusters in STUDY
mode = 'centroid'; % plot clusters centroid 
figureon = 1; % plot on a new figure
for k = 3:2:nargin
    switch varargin{k-2}
        case 'clusters'
            if isnumeric(varargin{k-1})
                cls = varargin{k-1};
                if isempty(cls)
                    cls = 2:length(STUDY.cluster);
                end
            else
                if isstr(varargin{k-1}) & strcmpi(varargin{k-1}, 'all')
                    cls = 2:length(STUDY.cluster);
                else
                    error('std_specplot(): ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'comps'
            STUDY = std_plotcompspec(STUDY, ALLEEG,  cls, varargin{k-1});
            return;
        case 'mode' % Plotting mode 'centroid' / 'comps'
            mode = varargin{k-1};
        case 'figure'
            if strcmpi(varargin{k-1},'off')
                figureon = 0;
            end
    end
end

% select clusters to plot
% -----------------------
if isempty(cls)
    tmp =[];
    cls = 2:length(STUDY.cluster); % plot all clusters in STUDY
    for k = 1: length(cls)
        % don't include 'Notclust' clusters
        if ~strncmpi('Notclust',STUDY.cluster(cls(k)).name,8) & ~strncmpi('ParentCluster',STUDY.cluster(cls(k)).name,13)
            tmp = [tmp cls(k)];
        end
    end
    cls = tmp;
end;

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond =1;
end
% Plot all the components in the cluster ('comps' mode)
if strcmpi(mode, 'comps')         
    for clus = 1: length(cls) % For each cluster requested
        len = length(STUDY.cluster(cls(clus)).comps);
        if ~isfield(STUDY.cluster(cls(clus)).centroid,'spec')
            STUDY = std_centroid(STUDY,ALLEEG, cls(clus) , 'spec');
        end
        
        % figure properties
        % -----------------
        figure
        rowcols(2) = ceil(sqrt(Ncond)); 
        rowcols(1) = ceil((Ncond)/rowcols(2));
        pos = get(gcf, 'position');
        magnif = 2.5/sqrt(Ncond);
        set(gcf, 'position', [ pos(1)+15 pos(2)+15 pos(3)*magnif pos(4)/rowcols(2)*rowcols(1)*magnif ]);
        orient tall
        set(gcf,'Color', BACKCOLOR);        
        try
            clusnval = std_clustread(STUDY, ALLEEG, cls(clus),'spec',1:Ncond);
        catch,
            warndlg2([ 'Not all component spectra found: aborting'] , ['Abort - Plot spectra'] );   
            return;
        end
        
        for n = 1:Ncond
            handl(n) = sbplot(rowcols(1),rowcols(2),n);
            ave_spec = STUDY.cluster(cls(clus)).centroid.spec{n};
            f = STUDY.cluster(cls(clus)).centroid.spec_freqs;
            for index = 1:size(clusnval.spec,2)
                plot(f,clusnval.spec{n, index},'color', [0.5 0.5 0.5]);
                hold on
            end;
            plot(f,ave_spec,'k','linewidth',2);
            xlabel('Frequency [Hz]');
            ylabel('Power [dB]');
            title(['Spectra, '  STUDY.cluster(cls(clus)).name ', ' STUDY.condition{n} ', ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) 'Ss']);
            if n == 1
                ylimits = get(gca,'YLim');
            else
                tmp = get(gca,'YLim');
                ylimits(1) = min(tmp(1),ylimits(1) );
                ylimits(2) = max(tmp(2),ylimits(2) );
            end
            if n == Ncond %set all condition figures to be on the same scale
                for condi = 1: Ncond
                    axes(handl(condi));
                    axis([f(1) f(end)  ylimits(1)  ylimits(2) ]);
                    axcopy;
                end
            end
        end % finished one condition
    end % finished all requested clusters 
end % Finished 'comps' mode plot option
       
% Plot clusters mean spec
if strcmpi(mode, 'centroid') 
    len = length(cls);
    rowcols(2) = ceil(sqrt(len)); rowcols(1) = ceil((len)/rowcols(2));
    if figureon
        try 
            % optional 'CreateCancelBtn', 'delete(gcbf); error(''USER ABORT'');', 
            h_wait = waitbar(0,['Computing spectra ...'], 'Color', BACKEEGLABCOLOR,'position', [300, 200, 300, 48]);
        catch % for Matlab 5.3
            h_wait = waitbar(0,['Computing spectra ...'],'position', [300, 200, 300, 48]);
        end
        figure
	end
    color_codes = {'b', 'r', 'g', 'c', 'm', 'y', 'k','b--', 'r--', 'g--', 'c--', 'm--', 'y--', 'k--','b-.', 'r-.', 'g-.', 'c-.', 'm-.', 'y-.', 'k-.'};
    orient tall
    
    min_spec = Inf;
    max_spec = -Inf;
    for k = 1:len % Go through the clusters
        if ~isfield(STUDY.cluster(cls(k)).centroid,'spec')
            STUDY = std_centroid(STUDY,ALLEEG, cls(k) , 'spec');
        end
        if len ~= 1
            sbplot(rowcols(1),rowcols(2),k) ; 
        end
        hold on;
        for n = 1:Ncond
            if k == 1
                leg_color{n} = [STUDY.condition{n}];
            end
            ave_spec = STUDY.cluster(cls(k)).centroid.spec{n};
            f = STUDY.cluster(cls(k)).centroid.spec_freqs;
            plot(f,ave_spec,color_codes{n},'linewidth',2);
            
            min_spec = min(min(ave_spec), min_spec);
            max_spec = max(max(ave_spec), max_spec);
            
            if n == Ncond
                a = [ STUDY.cluster(cls(k)).name ', '  num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ];
                title(a);
                set(gcf,'Color', BACKCOLOR);
                set(gca,'UserData', leg_color);
                set(gcf,'UserData', leg_color);
                xlim([f(1) f(end)]);
                axcopy(gcf, 'leg_color = get(gca,''''UserData'''') ; legend(leg_color); xlabel(''''Frequency [Hz]'''');ylabel(''''Power [dB]'''') ;');
                if figureon
                    waitbar(k/len,h_wait);
                end
            end
            if (all(rowcols == 1) | (mod(k, rowcols(2)) == 1 & (floor(k/rowcols(2)) == rowcols(1)-1))) & n == Ncond
                xlabel('Frequency [Hz]');
                ylabel('Power [dB]');
                if len ~= 1
                    maintitle = ['Mean cluster spectra across all conditions'];
                    a = textsc(maintitle, 'title'); 
                    set(a, 'fontweight', 'bold'); 
                else
                    a = [ STUDY.cluster(cls(k)).name ' spectra, '  num2str(length(unique(STUDY.cluster(cls(k)).sets(1,:)))) 'Ss' ];
                    title(a);
                end
                diff_spec = max_spec-min_spec;
                ylim( [ min_spec - 0.1*diff_spec max_spec + 0.1*diff_spec]);
                set(gcf,'Color', BACKCOLOR);
                legend(leg_color);
            end
        end % finished the different conditions
    end % finished all clusters 
    if figureon
        delete(h_wait);
    end
end % finished 'centroid' plot mode

% std_plotcompspec() - Commandline function, to visualizing cluster component spectra. 
%                   Displays the spectra of specified cluster components with the cluster mean 
%                   spectra on separate figures, using one figure for all conditions. 
%                   The spectra can be visualized only if component spectra     
%                   were calculated and saved in the EEG datasets in the STUDY.
%                   These can be computed during pre-clustering using the GUI-based function
%                   pop_preclust() or the equivalent commandline functions eeg_createdata() 
%                   and eeg_preclust(). A pop-function that calls this function is pop_clustedit().
% Usage:    
%                   >> [STUDY] = std_plotcompspec(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster
%                                                      (as in std_specplot). {default: 'all'}.
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster
%                     spectrum mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= 'all';  
%                         >> [STUDY] = std_plotcompspec(STUDY,ALLEEG, cluster, comps);
%                    Plots all components of cluster 4, calls std_specplot() . 
%
%  See also  pop_clustedit, pop_preclust, eeg_createdata, std_specplot         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 07, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_plotcompspec(STUDY, ALLEEG, cls, varargin)
icadefs;

if ~exist('cls')
    error('std_plotcompspec(): you must provide a cluster number as an input.');
end
if isempty(cls)
   error('std_plotcompspec(): you must provide a cluster number as an input.');
end
if nargin == 3 % no components indices were given
    % Default plot all components of the cluster
    [STUDY] = std_specplot(STUDY, ALLEEG, 'clusters', cls, 'mode', 'comps');
    return
else
    comp_ind = varargin{1}; 
end

Ncond = length(STUDY.condition);
if Ncond == 0
    Ncond =1;
end
for ci = 1 : length(comp_ind) %for each comp
    rowcols(2) = ceil(sqrt(Ncond)); rowcols(1) = ceil((Ncond)/rowcols(2));
    comp = STUDY.cluster(cls).comps(comp_ind(ci));     
    figure
    orient tall
    set(gcf,'Color', BACKCOLOR);
    subject  = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).subject;
    if Ncond >1
        textsc(['Spectra, ' subject ' / IC' num2str(comp) ', ' STUDY.cluster(cls).name ],'title');
    end
    for n = 1:Ncond  %for each cond
        abset = STUDY.datasetinfo(STUDY.cluster(cls).sets(n,comp_ind(ci))).index;
        sbplot(rowcols(1),rowcols(2),n),
        hold on
        if Ncond  > 1
            a = [ 'IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name ', ' STUDY.condition{n} ];
        else
            a = [ 'Spectra, IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name ];
        end
        if ~isfield(STUDY.cluster(cls).centroid,'spec')
            STUDY = std_centroid(STUDY,ALLEEG, cls, 'spec');
        end
        [spec, f] = std_readspec(ALLEEG, abset, comp, STUDY.preclust.specclustfreqs);
        plot(f,spec, 'c');
        ave_spec = STUDY.cluster(cls).centroid.spec{n};
        plot(f,ave_spec,'k','linewidth',2);
        xlim([f(1) f(end)]);
        xlabel('Frequency [Hz]');
        ylabel('Power [dB]');
        title(a);
        axcopy
    end
end

