% std_dipplot() - Commandline function to plot cluster component dipoles. Dipoles for each
%                 named cluster is displayed in a separate figure. To view all the clustered 
%                 components in the STUDY on the same figure (in a separate subplot), all 
%                 STUDY clusters must be requested.
%                 To visualize dipoles, they first must be stored in the EEG dataset structures
%                 using dipfit(). Only components that have dipole locations will be displayed,
%                 along with the cluster mean dipole (in red). 
% Usage:    
%                 >> [STUDY] = std_dipplot(STUDY, ALLEEG, clusters);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in 
%                the STUDY. ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%
% Optional inputs:
%   'clusters' - [numeric vector]  -> specific cluster numbers to plot.
%                  'all'         -> plot all clusters in STUDY.
%                  {default: 'all'}.
%   'comps'    - [numeric vector]  -> indices of the cluster components to plot.
%                  'all'  -> plot all the components in the cluster {default: 'all'}.
%   'mode'     - ['joined'|'apart'] Display all requested cluster on one 
%                  figure or separate figures. {default: 'joined'}
%                  'joined' -> plot all 'clusters' in one figure (without the gui).
%                  'apart'   -> plot each cluster in a separate figure. Note
%                  that this parameter has no effect if the 'comps' option is used.
%   'figure'   - ['on'|'off'] plots on a new figure ('on')  or plots on current
%                  figure ('off'). If figure 'off' does not display gui controls,
%                  Useful for incomporating one cluster dipplot into a 
%                  complex figure. {default: 'on'}. 
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster 
%                  mean dipole, to allow quick replotting (unless cluster means 
%                  already exists in the STUDY).  
%   Example:
%      >> [STUDY] = std_dipplot(STUDY,ALLEEG, 'clusters', 5, 'mode', 'apart', 'figure', 'off');
%                   % Plot cluster-5 component dipoles (in blue), plus ther mean dipole (in red), 
%                   % on an exisiting (gui-less) figure. 
%
%  See also  pop_clustedit(), dipplot()        
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 08, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_dipplot(STUDY, ALLEEG, varargin)

% Set default values
cls = []; % plot all clusters in STUDY
mode = 'apart'; % plot clusters on separate figures  
figureon = 1; % plot on a new figure
opt_dipplot = {'normlen', 'on', 'pointout', 'on', 'verbose', 'off', 'dipolelength', 0};
%, 'spheres', 'on'
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
                    error('std_dipplot: ''clusters'' input takes either specific clusters (numeric vector) or keyword ''all''.');
                end
            end
        case 'comps'
            STUDY = std_plotcompdip(STUDY, ALLEEG,  cls, varargin{k-1});
            return;
        case 'mode' % Plotting mode 'apart' / 'joined'
            mode = varargin{k-1};
            if strcmpi(mode, 'centroid') | strcmpi(mode, 'comps')
                if length(cls) == 1, mode = 'apart'; else mode = 'joined'; end;
            end;
         case 'figure'
            if strcmpi(varargin{k-1},'off') 
                opt_dipplot{end + 1} = 'gui';
                opt_dipplot{end + 1} = 'off';
                figureon = 0;
            else
                opt_dipplot{end + 1} = 'projimg';
                opt_dipplot{end + 1} = 'on';
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

if strcmpi(mode, 'apart')  % case each cluster on a separate figure
    for clus = 1: length(cls) % For each cluster requested
        if length(STUDY.cluster(cls(clus)).comps) > 0  % check there are comps in cluster
            max_r = 0;
            clear cluster_dip_models;
            len = length(STUDY.cluster(cls(clus)).comps);
            ndip = 0;
            dip_ind = [];
            if ~isfield(STUDY.cluster(cls(clus)).centroid,'dipole')
                STUDY = std_centroid(STUDY,ALLEEG, cls(clus) , 'dipole');
            end
            for k = 1:len
                abset   = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).index;
                subject = STUDY.datasetinfo(STUDY.cluster(cls(clus)).sets(1,k)).subject;
               if ~isfield(ALLEEG(abset), 'dipfit')
                   warndlg2(['No dipole information available in dataset ' ALLEEG(abset).filename ' , abort plotting'], 'Aborting plot dipoles');
                   return;
               end
               comp = STUDY.cluster(cls(clus)).comps(k); 
               cluster_dip_models(k) = ALLEEG(abset).dipfit.model(comp);
               if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                   if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                       load('-mat', ALLEEG(abset).dipfit.hdmfile);
                       max_r = max(max_r, max(vol.r));
                   else % old version of dipfit
                       max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
                   end
               end
               comp_to_disp{k} = ['IC' num2str(comp) ', ' subject];
               if ~isempty(cluster_dip_models(k).posxyz)
                   ndip = ndip +1;
                   dip_ind = [dip_ind k];
               end
           end % finished going over cluster comps
           cluster_dip_models(end + 1) = STUDY.cluster(cls(clus)).centroid.dipole;
           
           % additional options
           % ------------------
           dip_color = cell(1,ndip+1);
           dip_color(1:ndip) = {'b'};
           dip_color(end) = {'r'};
           options = opt_dipplot;
           options{end+1} =  'mri';
           options{end+1} =  ALLEEG(abset).dipfit.mrifile;
           options{end+1} =  'coordformat';
           options{end+1} =  ALLEEG(abset).dipfit.coordformat;
           options{end+1} =  'dipnames';
           options{end+1} = {comp_to_disp{dip_ind } [STUDY.cluster(cls(clus)).name ' mean']};
           options{end+1} = 'color';
           options{end+1} = dip_color;
           if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
               options{end+1} = 'sphere';
               options{end+1} = max_r;
           else
               options{end+1} = 'meshdata';
               options{end+1} = ALLEEG(abset).dipfit.hdmfile;
           end
           if ndip < 6 % less than 6 dipoles, plot projected lines 
               options{end+1} = 'projlines';
               options{end+1} = 'on';
           end
           
           if figureon
               dipplot(cluster_dip_models, options{:});
               fig_h = gcf;
               set(fig_h,'Name', [STUDY.cluster(cls(clus)).name ' - ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) ...
                       ' sets - ' num2str(length(STUDY.cluster(cls(clus)).comps)) ' components (' num2str(ndip) ' dipoles)' ],'NumberTitle','off');
           else
               dipplot(cluster_dip_models, options{:},'view', [0.5 -0.5 0.5],'spheres', 'on', 'projlines', 'on');
               for gind = 1:length(options) % remove the 'gui' 'off' option
                   if isstr(options{gind}) 
                       if strfind(options{gind}, 'gui')
                           break;
                       end
                   end
               end
               options(gind:gind+1) = [];
               dipinfo.dipmod =  cluster_dip_models;
               dipinfo.op = options;
               diptitle = [STUDY.cluster(cls(clus)).name ', ' num2str(length(unique(STUDY.cluster(cls(clus)).sets(1,:)))) ' sets -' ...
                   num2str(length(STUDY.cluster(cls(clus)).comps)) ' components (' num2str(ndip) ' dipoles)' ];
               dipinfo.title = diptitle;
               set(gcf, 'UserData', dipinfo);
               set(gca,'UserData', dipinfo);
               rotate3d off;
               axcopy(gcf, ['dipinfo = get(gca, ''''UserData''''); dipplot(dipinfo.dipmod, dipinfo.op{:}); set(gcf, ''''Name'''', dipinfo.title,''''NumberTitle'''',''''off''''); ']);
           end
        end % finished the if condition that cluster isn't empty
    end % finished going over requested clusters
end

if strcmpi(mode, 'joined')  % case all clusters are plotted in the same figure (must be a new figure)
    N = length(cls);
    rowcols(2) = ceil(sqrt(N)); % Number of rows in the subplot figure.
    rowcols(1) = ceil(N/rowcols(2));
    fig_h = figure; 
    orient tall
    set(fig_h,'Color', 'black');
    set(fig_h,'Name', 'All clusters dipoles','NumberTitle','off');
    set(fig_h, 'resize','off');
    for l = 1:N
        len = length(STUDY.cluster(cls(l)).comps);
        max_r = 0;
        clear cluster_dip_models;
        if ~isfield(STUDY.cluster(cls(l)).centroid,'dipole')
            STUDY = std_centroid(STUDY,ALLEEG, cls(l) , 'dipole');
        end
        for k = 1: len
            abset = STUDY.datasetinfo(STUDY.cluster(cls(l)).sets(1,k)).index;
           if ~isfield(ALLEEG(abset), 'dipfit')
               warndlg2(['No dipole information available in dataset ' num2str(abset) ' , abort plotting'], 'Aborting plot dipoles');
               return;
           end
           comp = STUDY.cluster(cls(l)).comps(k); 
           cluster_dip_models(k) = ALLEEG(abset).dipfit.model(comp);
           if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
                if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                   load('-mat', ALLEEG(abset).dipfit.hdmfile);
                   max_r = max(max_r, max(vol.r));
               else % old version of dipfit
                   max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
               end
           end
        end % finished going over cluster comps
        cluster_dip_models(end + 1) = STUDY.cluster(cls(l)).centroid.dipole;
        dip_color = cell(1,length(cluster_dip_models));
        dip_color(1:end-1) = {'b'};
        dip_color(end) = {'r'};
        options = opt_dipplot;
        options{end + 1} =  'gui';
        options{end + 1} =  'off';
        options{end+1} =  'mri';
        options{end+1} =  ALLEEG(abset).dipfit.mrifile;
        options{end+1} =  'coordformat';
        options{end+1} =  ALLEEG(abset).dipfit.coordformat;
        options{end+1} = 'color';
        options{end+1} = dip_color;
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
            options{end+1} = 'sphere';
            options{end+1} = max_r;
        else
            options{end+1} = 'meshdata';
            options{end+1} = ALLEEG(abset).dipfit.hdmfile;
        end
        subplot(rowcols(1),rowcols(2),l) , 
        dipplot(cluster_dip_models, options{:});
        diptitle = [STUDY.cluster(cls(l)).name ', ' num2str(length(unique(STUDY.cluster(cls(l)).sets(1,:)))) 'Ss'];
        title(diptitle, 'Color', 'white');
        % Complex axcopy
        %if l == 1
        %    for gind = 1:length(options) % remove the 'gui' 'off' option
        %        if isstr(options{gind}) 
        %            if strfind(options{gind}, 'gui')
        %                break;
        %            end
        %        end
        %    end
        %    options(gind:gind+1) = [];
        %end
        %dipinfo.dipmod =  cluster_dip_models;
        %dipinfo.op = options;
        %dipinfo.title = diptitle;
        %set(gcf, 'UserData', dipinfo);
        %set(gca,'UserData', dipinfo);
        %axcopy(gcf, ['dipinfo = get(gca, ''''UserData''''); dipplot(dipinfo.dipmod, dipinfo.op{:}); set(gcf, ''''Name'''', dipinfo.title,''''NumberTitle'''',''''off'''');']);
   end %finished going over all clusters
   set(fig_h, 'resize','on');
end % finished case of 'all' clusters
% std_plotcompdip() - Commandline function, to visualizing cluster components dipoles. 
%                   Displays the dipoles of specified cluster components with the cluster mean 
%                   dipole on separate figures. 
%                   To visualize dipoles they first must be stored in the EEG dataset structures
%                   using dipfit(). Only components that have a dipole locations will be displayed,
%                   along with the cluster mean dipole in red. 
% Usage:    
%                   >> [STUDY] = std_plotcompdip(STUDY, ALLEEG, cluster, comps);  
% Inputs:
%   STUDY      - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG     - global EEGLAB vector of EEG structures for the dataset(s) included in the STUDY. 
%                     ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   cluster     - single cluster number.  
%
% Optional inputs:
%   comps      - [numeric vector]  -> indices of the cluster components to plot.
%                       'all'                       -> plot all the components in the cluster {default: 'all'}. 
%
% Outputs:
%   STUDY    - the input STUDY set structure modified with plotted cluster
%                     dipole mean, to allow quick replotting (unless cluster mean 
%                     already existed in the STUDY).  
%
%   Example:
%                         >> cluster = 4; comps= 1;  
%                         >> [STUDY] = std_plotcompdip(STUDY,ALLEEG, cluster, comps);
%                    Plots component 1 dipole in blue with the cluster 4 mean dipole in red. 
%
%  See also  pop_clustedit, dipfit, std_dipplot         
%
% Authors:  Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005

%123456789012345678901234567890123456789012345678901234567890123456789012

% Copyright (C) Hilit Serby, SCCN, INC, UCSD, June 08, 2005, hilit@sccn.ucsd.edu
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

function STUDY = std_plotcompdip(STUDY, ALLEEG, cls, varargin)
if ~exist('cls')
    error('std_plotcompdip: you must provide a cluster number as an input.');
end
if isempty(cls)
   error('std_plotcompdip: you must provide a cluster number as an input.');
end
if nargin == 3 % no components indices were given
    % Default plot all components of the cluster
    [STUDY] = std_dipplot(STUDY, ALLEEG, 'clusters', cls);
    return
else
    comp_ind = varargin{1}; 
end

for ci = 1:length(comp_ind)
    abset = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).index;
    comp = STUDY.cluster(cls).comps(comp_ind(ci));
    subject = STUDY.datasetinfo(STUDY.cluster(cls).sets(1,comp_ind(ci))).subject;
    if ~isfield(ALLEEG(abset), 'dipfit')
        warndlg2(['No dipole information available in dataset ' num2str(abset) ' , abort plotting'], 'Aborting plot dipoles');
        return;
    end
    if ~isfield(STUDY.cluster(cls).centroid,'dipole')
        STUDY = std_centroid(STUDY,ALLEEG, cls , 'dipole');
    end
    comp_to_disp = ['IC' num2str(comp) ' / ' subject];
    cluster_dip_models = ALLEEG(abset).dipfit.model(comp);
    cluster_dip_models(2) = STUDY.cluster(cls).centroid.dipole;
    if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
        if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
            load('-mat', ALLEEG(abset).dipfit.hdmfile);
            max_r = max(vol.r);
        else
            max_r = max(ALLEEG(abset).dipfit.vol.r);
        end
        dipplot(cluster_dip_models, 'sphere', max_r, 'mri', ALLEEG(abset).dipfit.mrifile,'coordformat', ALLEEG(abset).dipfit.coordformat , ...
           'normlen' ,'on', 'pointout' ,'on','color', {'b', 'r'}, 'dipnames', {comp_to_disp [ STUDY.cluster(cls).name ' mean' ] },...
            'spheres', 'on', 'verbose', 'off');
    else
       dipplot(cluster_dip_models, 'meshdata', ALLEEG(abset).dipfit.hdmfile, 'mri', ALLEEG(abset).dipfit.mrifile,'coordformat', ALLEEG(abset).dipfit.coordformat , ...
          'normlen' ,'on', 'pointout' ,'on','color', {'b', 'r'}, 'dipnames', {comp_to_disp [STUDY.cluster(cls).name ' mean']}, ...
          'spheres', 'on', 'verbose', 'off');
    end
    fig_h = gcf;
    set(fig_h,'Name', [ 'IC' num2str(comp) ' / ' subject ', ' STUDY.cluster(cls).name],'NumberTitle','off');
end
        
