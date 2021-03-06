function [cfg] = sourceplot(cfg, data)

% SOURCEPLOT plots functional source reconstruction data on slices or on a
% surface, optionally as an overlay on anatomical MRI data, where
% statistical data can be used to determine the opacity of the mask.
%
% Use as:
%   sourceplot(cfg, data)
%
% The data can contain functional data, anatomical MRI data and statistical data,
% interpolated onto the same grid.
%
% The configuration should contain:
%   cfg.method        = 'slice',   plots the data on a number of slices in the same plane
%                       'ortho',   plots the data on three orthogonal slices
%                       'surface', plots the data on a 3D brain surface
%
%   cfg.anaparameter  = string, field in data with the anatomical data (default = 'anatomy' if present in data)
%   cfg.funparameter  = string, field in data with the functional parameter of interest (default = [])
%   cfg.maskparameter = string, field in the data to be used for opacity masking of fun data (default = [])
%                        If values are between 0 and 1, zero is fully transparant and one is fully opaque.
%                        If values in the field are not between 0 and 1 they will be scaled depending on the values
%                        of cfg.opacitymap and cfg.opacitylim (see below)
%                        You can use masking in several ways, f.i.
%                        - use outcome of statistics to show only the significant values and mask the insignificant
%                          NB see also cfg.opacitymap and cfg.opacitylim below
%                        - use the functional data itself as mask, the highest value (and/or lowest when negative)
%                          will be opaque and the value closest to zero transparent
%                        - Make your own field in the data with values between 0 and 1 to control opacity directly
%
% The following parameters can be used in all methods:
%   cfg.downsample    = downsampling for resolution reduction, integer value (default = 1) (orig: from surface)
%
% The following parameters can be used for the functional data:
%   cfg.funcolormap   = colormap for functional data, see COLORMAP (default = 'auto')
%                       'auto', depends structure funparameter, or on funcolorlim
%                         - funparameter: only positive values, or funcolorlim:'zeromax' -> 'hot'
%                         - funparameter: only negative values, or funcolorlim:'minzero' -> 'cool'
%                         - funparameter: both pos and neg values, or funcolorlim:'maxabs' -> 'jet'
%                         - funcolorlim: [min max] if min & max pos-> 'hot', neg-> 'cool', both-> 'jet'
%   cfg.funcolorlim   = color range of the functional data (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(funparameter)) to +max(abs(funparameter))
%                        'zeromax', from 0 to max(abs(funparameter))
%                        'minzero', from min(abs(funparameter)) to 0
%                        'auto', if funparameter values are all positive: 'zeromax',
%                          all negative: 'minzero', both possitive and negative: 'maxabs'
%   cfg.colorbar      = 'yes' or 'no' (default = 'yes')
%
% The following parameters can be used for the masking data:
%   cfg.opacitymap    = opacitymap for mask data, see ALPHAMAP (default = 'auto')
%                       'auto', depends structure maskparameter, or on opacitylim
%                         - maskparameter: only positive values, or opacitylim:'zeromax' -> 'rampup'
%                         - maskparameter: only negative values, or opacitylim:'minzero' -> 'rampdown'
%                         - maskparameter: both pos and neg values, or opacitylim:'maxabs' -> 'vdown'
%                         - opacitylim: [min max] if min & max pos-> 'rampup', neg-> 'rampdown', both-> 'vdown'
%                         - NB. to use p-values use 'rampdown' to get lowest p-values opaque and highest transparent
%   cfg.opacitylim    = range of mask values to which opacitymap is scaled (default = 'auto')
%                        [min max]
%                        'maxabs', from -max(abs(maskparameter)) to +max(abs(maskparameter))
%                        'zeromax', from 0 to max(abs(maskparameter))
%                        'minzero', from min(abs(maskparameter)) to 0
%                        'auto', if maskparameter values are all positive: 'zeromax',
%                          all negative: 'minzero', both possitive and negative: 'maxabs'
%
% The folowing parameters apply for ortho-plotting
%   cfg.location      = location of cut, (default = 'auto')
%                        'auto', 'center' if only anatomy, 'max' if functional data
%                        'min' and 'max' position of min/max funparameter
%                        'center' of the brain
%                        [x y z], coordinates in voxels or head, see cfg.locationcoordinates
%   cfg.locationcoordinates = coordinate system used in cfg.location, 'head' or 'voxel' (default = 'head')
%                              'head', headcoordinates from anatomical MRI
%                              'voxel', voxelcoordinates
%   cfg.crosshair     = 'yes' or 'no' (default = 'yes')
%   cfg.axis          = 'on' or 'off' (default = 'on')
%   cfg.interactive   = 'yes' or 'no' (default = 'no')
%                        in interactive mode cusor click determines location of cut
%
% The folowing parameters apply for slice-plotting
%   cfg.nslices       = number of slices, (default = 20)
%   cfg.slicerange    = range of slices in data, (default = 'auto')
%                       'auto', full range of data
%                       [min max], coordinates of first and last slice in voxels
%   cfg.slicedim      = dimension to slice 1 (x-axis) 2(y-axis) 3(z-axis) (default = 3)
%
% The folowing parameters apply for surface-plotting
%   cfg.surffile       = string, file that contains the surface (default = 'single_subj_T1.mat')
%                        'single_subj_T1.mat' contains a triangulation that corresponds with the
%                         SPM anatomical template in MNI coordinates
%   cfg.surfinflated   = string, file that contains the inflated surface (default = [])
%   cfg.surfdownsample = number (default = 1, i.e. no downsampling)
%   cfg.projmethod     = projection method, how functional volume data is projected onto surface
%                        'nearest', 'sphere_avg', 'sphere_weighteddistance'
%   cfg.sphereradius   = maximum distance from each voxel to the surface to be
%                        included in the sphere projection methods, expressed in mm
%   cfg.distmat        = precomputed distance matrix (default = [])
%   cfg.camlight       = 'yes' or 'no' (default = 'yes')
%
% TODO have to be built in:
%   cfg.marker        = [Nx3] array defining N marker positions to display (orig: from sliceinterp)
%   cfg.markersize    = radius of markers (default = 5)
%   cfg.markercolor   = [1x3] marker color in RGB (default = [1 1 1], i.e. white) (orig: from sliceinterp)
%   cfg.title         = optional title (default is '') (orig: from sliceinterp)
%   TT lookup function
%   white background option

% undocumented TODO
%   slice in all directions
%   surface also optimal when inside present
%   anamask possibilities

% Copyright (C) 2007, Robert Oostenveld
%
% $Log: sourceplot.m,v $
% Revision 1.51  2007/11/26 09:08:57  roboos
% give instructive error if nothing is selected to be plotted
%
% Revision 1.50  2007/06/06 08:47:19  ingnie
% added possibility to turn lights off by cfg.camlight
%
% Revision 1.49  2007/05/06 12:10:08  roboos
% switched from interp_ungridded to interp_gridded, which is much faster
% added cfg.surfinflated for inflated surface
% use surface curvature (if available) for color
%
% Revision 1.48  2007/05/02 15:23:39  roboos
% only changes in whitespace
%
% Revision 1.47  2007/05/02 09:20:41  ingnie
% fixed bug; make local dim variable after downsampling
%
% Revision 1.46  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.45  2007/03/30 16:53:03  ingnie
% added checkinput; only proceed when data is allowed datatype, added error in method 'slice' when no ana and no fun, fixed typo in error
%
% Revision 1.44  2007/03/20 07:24:17  roboos
% comverted from DOS ascii to UNIX ascii
%
% Revision 1.43  2007/03/08 17:04:53  ingnie
% added colorbar option for surfaceplot
%
% Revision 1.42  2007/03/08 16:26:05  ingnie
% fixed bug in surface plotting; added hasfun and hasmsk around colormap and alphamap setting
%
% Revision 1.41  2007/02/22 18:02:59  ingnie
% added intelligent mask which prevents outside from being plotted if no maskparam specified, fixed bug that appeared when mask islogical, only transparent fun when hasana
%
% Revision 1.40  2007/02/16 16:12:02  ingnie
% added cla, fixed bug with opacitymap
%
% Revision 1.39  2007/02/08 16:47:50  ingnie
% Whole new implementation; combining sliceinterp, surfaceplot and sourceplot_old
%

%% checkdata see below!!! %%

% set the common defaults
if ~isfield(cfg, 'method'),              cfg.method = 'ortho';              end %added ingnie
if ~isfield(cfg, 'anaparameter'),
  if isfield(data, 'anatomy'),
    cfg.anaparameter = 'anatomy';
  else
    cfg.anaparameter = [];
  end
end % added ingnie

% all methods
if ~isfield(cfg, 'funparameter'),        cfg.funparameter = [];             end
if ~isfield(cfg, 'maskparameter'),       cfg.maskparameter = [];            end
if ~isfield(cfg, 'downsample'),          cfg.downsample = 1;                end
if ~isfield(cfg, 'title'),               cfg.title = '';                    end
if ~isfield(cfg, 'marker'),              cfg.marker = [];                   end
if ~isfield(cfg, 'markersize'),          cfg.markersize = 5;                end
if ~isfield(cfg, 'markercolor'),         cfg.markercolor = [1,1,1];         end

% set the common defaults for the functional data
if ~isfield(cfg, 'funcolormap'),         cfg.funcolormap = 'auto';          end  % was cfg.colormap
if ~isfield(cfg, 'funcolorlim'),         cfg.funcolorlim = 'auto';          end; % was min max iso range

% set the common defaults for the statistical data
if ~isfield(cfg, 'opacitymap'),         cfg.opacitymap = 'auto';            end; % was maskmap
if ~isfield(cfg, 'opacitylim'),         cfg.opacitylim = 'auto';            end; % added ingnie, iso whole bunch of other cfg

% set the defaults per method
% ortho
if ~isfield(cfg, 'location'),            cfg.location = 'auto';              end
if ~isfield(cfg, 'locationcoordinates'), cfg.locationcoordinates = 'head';   end
if ~isfield(cfg, 'crosshair'),           cfg.crosshair = 'yes';              end
if ~isfield(cfg, 'colorbar'),            cfg.colorbar  = 'yes';              end
if ~isfield(cfg, 'axis'),                cfg.axis   = 'on';                  end
if ~isfield(cfg, 'interactive'),         cfg.interactive   = 'no';           end
% slice
if ~isfield(cfg, 'nslices');            cfg.nslices = 20;                    end
if ~isfield(cfg, 'slicedim');           cfg.slicedim = 3;                    end
if ~isfield(cfg, 'slicerange');         cfg.slicerange = 'auto';             end
% surface
if ~isfield(cfg, 'downsample'),         cfg.downsample     = 1;              end
if ~isfield(cfg, 'surfdownsample'),     cfg.surfdownsample = 1;              end
if ~isfield(cfg, 'surffile'),           cfg.surffile = 'single_subj_T1.mat'; end % use a triangulation that corresponds with the collin27 anatomical template in MNI coordinates
if ~isfield(cfg, 'surfinflated'),       cfg.surfinflated = [];               end
if ~isfield(cfg, 'sphereradius'),       cfg.sphereradius = [];               end
if ~isfield(cfg, 'distmat'),            cfg.distmat = [];                    end
if ~isfield(cfg, 'camlight'),           cfg.camlight = 'yes';                end

if isequal(cfg.method,'surface')
  if ~isfield(cfg, 'projmethod'),         error('specify cfg.projmethod');     end
end

% for backward compatibility
if strcmp(cfg.location, 'interactive')
  cfg.location = 'auto';
  cfg.interactive = 'yes';
end

%%%%%%%
if isstr(data)
  % read the anatomical MRI data from file
  filename = data;
  fprintf('reading MRI from file\n');
  data = read_fcdc_mri(filename);
end

% convert the coordinates along the axes (i.e. xgrid/ygrid/zgrid) into a homogenous transformation matrix
data = grid2transform(data);

% check if the input data is valid for this function
data = checkdata(data, 'datatype', {'source', 'volume'}, 'feedback', 'yes');

% select the functional and the mask parameter
cfg.funparameter  = parameterselection(cfg.funparameter, data);
cfg.maskparameter = parameterselection(cfg.maskparameter, data);
% only a single parameter should be selected
try, cfg.funparameter  = cfg.funparameter{1};  end
try, cfg.maskparameter = cfg.maskparameter{1}; end

% downsample all volumes
tmpcfg = [];
tmpcfg.parameter  = {cfg.funparameter, cfg.maskparameter, 'anatomy'};
tmpcfg.downsample = cfg.downsample;
data = volumedownsample(tmpcfg, data);

%%% make the local variables:
dim = data.dim;

% put data.xgrid, data.ygrid, data.zgrid in data
% TODO Question ingnie WHY LIKE THIS? robert: goede vraag...
data.xgrid = 1:dim(1);
data.ygrid = 1:dim(2);
data.zgrid = 1:dim(3);

%%% handle anaparameter
if isequal(cfg.anaparameter,'anatomy')
  if isfield(data, 'anatomy')
    hasana = 1;
    mri8  = isa(data.anatomy, 'uint8');
    mri16 = isa(data.anatomy, 'uint16');
    % convert integers to single precision float if neccessary
    if mri8 || mri16
      fprintf('converting anatomy to double\n');
      ana = double(data.anatomy);
    else
      ana = data.anatomy;
    end
  else
    warning('no anatomical volume present, not plotting anatomy\n')
    hasana = 0;
  end
elseif isempty(cfg.anaparameter);
  hasana = 0;
  fprintf('not plotting anatomy\n');
else
  warning('do not understand cfg.anaparameter, not plotting anatomy\n')
  hasana = 0;
end

%%% funparameter
%% has fun?
if ~isempty(cfg.funparameter)
  if issubfield(data, cfg.funparameter)
    hasfun = 1;
    fun = getsubfield(data, cfg.funparameter);
  else
    error('cfg.funparameter not found in data'); % added ingnie
  end
else
  hasfun = 0;
  fprintf('no functional parameter\n');
end
%% handle fun
if hasfun
  % determine scaling min and max (fcolmin fcolmax) and funcolormap
  funmin = min(fun(:));
  funmax = max(fun(:));
  % smart lims: make from auto other string
  if isequal(cfg.funcolorlim,'auto') % ingnie included all new scaling options
    if sign(funmin)>-1 && sign(funmax)>-1
      cfg.funcolorlim = 'zeromax';
    elseif sign(funmin)<1 && sign(funmax)<1
      cfg.funcolorlim = 'minzero';
    else
      cfg.funcolorlim = 'maxabs';
    end
  end
  if ischar(cfg.funcolorlim)
    % limits are given as string
    if isequal(cfg.funcolorlim,'maxabs')
      fcolmin = -max(abs([funmin,funmax]));
      fcolmax =  max(abs([funmin,funmax]));
      if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'jet'; end;
    elseif isequal(cfg.funcolorlim,'zeromax')
      fcolmin = 0;
      fcolmax = funmax;
      if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'hot'; end;
    elseif isequal(cfg.funcolorlim,'minzero')
      fcolmin = funmin;
      fcolmax = 0;
      if isequal(cfg.funcolormap,'auto'); cfg.funcolormap = 'cool'; end;
    else
      error('do not understand cfg.funcolorlim');
    end
  else
    % limits are numeric
    fcolmin = cfg.funcolorlim(1); % ingnie changed to lim iso min
    fcolmax = cfg.funcolorlim(2);
    % smart colormap
    if isequal(cfg.funcolormap,'auto') % ingnie added smart colormap choice
      if sign(fcolmin) == -1 && sign(fcolmax) == 1
        cfg.funcolormap = 'jet';
      else
        if fcolmin < 0
          cfg.funcolormap = 'cool';
        else
          cfg.funcolormap = 'hot';
        end
      end
    end
  end %if ischar
  clear funmin funmax;
  % ingnie removed fun scaling, is not done anymore, now colormap is scaled iso data
  % ensure that the functional data is real
  if ~isreal(fun)
    fprintf('taking absolute value of complex data\n');
    fun = abs(fun);
  end
end % handle fun

%%% maskparameter
%% has mask?
if ~isempty(cfg.maskparameter)
  if issubfield(data, cfg.maskparameter)
    if ~hasfun
      error('you can not have a mask without functional data')
    else
      hasmsk = 1;
      msk = getsubfield(data, cfg.maskparameter);
      if islogical(msk) %otherwise sign() not posible
        msk = double(msk);
      end
    end
  else
    error('cfg.maskparameter not found in data'); %added ingnie
  end
else
  hasmsk = 0;
  fprintf('no masking parameter\n');
end
%% handle mask
if hasmsk
  % determine scaling and opacitymap
  mskmin = min(msk(:));
  mskmax = max(msk(:));
  % determine the opacity limits and the opacity map
  % smart lims: make from auto other string, or equal to funcolorlim if funparameter == maskparameter
  if isequal(cfg.opacitylim,'auto')
    if isequal(cfg.funparameter,cfg.maskparameter)
      cfg.opacitylim = cfg.funcolorlim;
    else
      if sign(mskmin)>-1 && sign(mskmax)>-1
        cfg.opacitylim = 'zeromax';
      elseif sign(mskmin)<1 && sign(mskmax)<1
        cfg.opacitylim = 'minzero';
      else
        cfg.opacitylim = 'maxabs';
      end
    end
  end
  if ischar(cfg.opacitylim)
    % limits are given as string
    switch cfg.opacitylim
      case 'zeromax'
        opacmin = 0;
        opacmax = mskmax;
        if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'rampup'; end;
      case 'minzero'
        opacmin = mskmin;
        opacmax = 0;
        if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'rampdown'; end;
      case 'maxabs'
        opacmin = -max(abs([mskmin, mskmax]));
        opacmax =  max(abs([mskmin, mskmax]));
        if isequal(cfg.opacitymap,'auto'), cfg.opacitymap = 'vdown'; end;
      otherwise
        error('incorrect specification of cfg.opacitylim');
    end
  else
    % limits are numeric
    opacmin = cfg.opacitylim(1);
    opacmax = cfg.opacitylim(2);
    if isequal(cfg.opacitymap,'auto')
      if sign(opacmin)>-1 && sign(opacmax)>-1
        cfg.opacitymap = 'rampup';
      elseif sign(opacmin)<1 && sign(opacmax)<1
        cfg.opacitymap = 'rampdown';
      else
        cfg.opacitymap = 'vdown';
      end
    end
  end % handling opacitylim and opacitymap
  clear mskmin mskmax;
end % ingnie removed mask scaling, is not done anymore, now opacitymap is scaled iso data

%% intelligence: make mask to prevent outside fun from being plotted
if hasfun && isfield(data,'inside') && ~hasmsk
  hasmsk = 1;
  msk = zeros(dim);
  cfg.opacitymap = 'rampup';
  opacmin = 0;
  opacmax = 1;
  % make intelligent mask
  if hasana
    msk(data.inside) = 0.5; %so anatomy is visible
  else
    msk(data.inside) = 1;
  end;
end;

%% set color and opacity mapping for this figure
if hasfun
  cfg.funcolormap = colormap(cfg.funcolormap);
  colormap(cfg.funcolormap);
end
if hasmsk
  cfg.opacitymap  = alphamap(cfg.opacitymap);
  alphamap(cfg.opacitymap);
end

%%% ensure that they are all 3D volumes
if hasana, ana = reshape(ana, dim); end;
if hasfun, fun = reshape(fun, dim); end;
if hasmsk, msk = reshape(msk, dim); end;

%%% determine what has to be plotted, depends on method
if isequal(cfg.method,'ortho')
  if ~isstr(cfg.location)
    if strcmp(cfg.locationcoordinates, 'head')
      % convert the headcoordinates location into voxel coordinates
      loc = inv(data.transform) * [cfg.location(:); 1];
      loc = round(loc(1:3));
    elseif strcmp(cfg.locationcoordinates, 'voxel')
      % the location is already in voxel coordinates
      loc = round(cfg.location(1:3));
    else
      error('you should specify cfg.locationcoordinates');
    end
  else
    if isequal(cfg.location,'auto')
      if hasfun
        if isequal(cfg.funcolorlim,'maxabs');
          loc = 'max';
        elseif isequal(cfg.funcolorlim, 'zeromax');
          loc = 'max';
        elseif isequal(cfg.funcolorlim, 'minzero');
          loc = 'min';
        else %if numerical
          loc = 'max';
        end
      else
        loc = 'center';
      end;
    else
      loc = cfg.location;
    end
  end

  % determine the initial intersection of the cursor (xi yi zi)
  if isstr(loc) && strcmp(loc, 'min')
    if isempty(cfg.funparameter)
      error('cfg.location is min, but no functional parameter specified');
    end
    [minval, minindx] = min(fun(:));
    [xi, yi, zi] = ind2sub(dim, minindx);
  elseif isstr(loc) && strcmp(loc, 'max')
    if isempty(cfg.funparameter)
      error('cfg.location is max, but no functional parameter specified');
    end
    [maxval, maxindx] = max(fun(:));
    [xi, yi, zi] = ind2sub(dim, maxindx);
  elseif isstr(loc) && strcmp(loc, 'center')
    xi = round(length(data.xgrid)/2);
    yi = round(length(data.ygrid)/2);
    zi = round(length(data.zgrid)/2);
  elseif ~isstr(loc)
    xi = nearest(data.xgrid, loc(1));
    yi = nearest(data.ygrid, loc(2));
    zi = nearest(data.zgrid, loc(3));
  end

  %% do the actual plotting %%

  interactive_flag = 1; % it happens at least once
  while(interactive_flag)
    interactive_flag = strcmp(cfg.interactive, 'yes');

    xi = round(xi); xi = max(xi, 1); xi = min(xi, dim(1));
    yi = round(yi); yi = max(yi, 1); yi = min(yi, dim(2));
    zi = round(zi); zi = max(zi, 1); zi = min(zi, dim(3));

    if interactive_flag
      fprintf('\n');
      fprintf('click with mouse button to reposition the cursor\n');
      fprintf('press q on keyboard to quit interactive mode\n');
    end

    ijk = [xi yi zi 1]';
    xyz = data.transform * ijk;
    if hasfun
      val = fun(xi, yi, zi);
      fprintf('voxel %d, indices [%d %d %d], location [%.1f %.1f %.1f], value %f\n', sub2ind(dim, xi, yi, zi), ijk(1:3), xyz(1:3), val);
    else
      fprintf('voxel %d, indices [%d %d %d], location [%.1f %.1f %.1f]\n', sub2ind(dim, xi, yi, zi), ijk(1:3), xyz(1:3));
    end

    %     if strcmp(cfg.TTlookup, 'yes')  %FIXME make this work
    %       lab = TTatlas_lookup(tlrc, mni2tal(xyz(1:3)), cfg.TTqueryrange);
    %       if isempty(lab)
    %         fprintf('Talairach-Tournoux labels: not found\n');
    %       else
    %         fprintf('Talairach-Tournoux labels: ')
    %         fprintf('%s', lab{1});
    %         for i=2:length(lab)
    %           fprintf(', %s', lab{i});
    %         end
    %         fprintf('\n');
    %       end
    %     end

    % make vols and scales, containes volumes to be plotted (fun, ana, msk) %added ingnie
    vols = {};
    if hasana; vols{1} = ana; scales{1} = []; end; % needed when only plotting ana
    if hasfun; vols{2} = fun; scales{2} = [fcolmin fcolmax]; end;
    if hasmsk; vols{3} = msk; scales{3} = [opacmin opacmax]; end;

    if isempty(vols)
      % this seems to be a problem that people often have
      error('no anatomy is present and no funcctional data is selected, please check your cfg.funparameter');
    end

    subplot(2,2,1);
    [vols2D] = handle_ortho(vols, [xi yi zi], 2, dim);
    plot2D(vols2D, scales);
    xlabel('i'); ylabel('k'); axis(cfg.axis);
    if strcmp(cfg.crosshair, 'yes'), crosshair([xi zi]); end

    subplot(2,2,2);
    [vols2D] = handle_ortho(vols, [xi yi zi], 1, dim);
    plot2D(vols2D, scales);
    xlabel('j'); ylabel('k'); axis(cfg.axis);
    if strcmp(cfg.crosshair, 'yes'), crosshair([yi zi]); end

    subplot(2,2,3);
    [vols2D] = handle_ortho(vols, [xi yi zi], 3, dim);
    plot2D(vols2D, scales);
    xlabel('i'); ylabel('j'); axis(cfg.axis);
    if strcmp(cfg.crosshair, 'yes'), crosshair([xi yi]); end

    if strcmp(cfg.colorbar,  'yes'),
      if hasfun
        % vectorcolorbar = linspace(fcolmin, fcolmax,length(cfg.funcolormap));
        % imagesc(vectorcolorbar,1,vectorcolorbar);colormap(cfg.funcolormap);
        subplot(2,2,4);
        % use a normal Matlab coorbar, attach it to the invisible 4th subplot
        caxis([fcolmin fcolmax]);
        hc = colorbar;
        set(hc, 'YLim', [fcolmin fcolmax]);
        set(gca, 'Visible', 'off');
      else
        warning('no colorbar possible without functional data')
      end
    end
    drawnow;

    if interactive_flag
      try, [d1, d2, key] = ginput(1); catch, key='q'; end
      if key=='q'
        break;
      else
        % update the view to a new position
        l1 = get(get(gca, 'xlabel'), 'string');
        l2 = get(get(gca, 'ylabel'), 'string');
        switch l1,
          case 'i'
            xi = d1;
          case 'j'
            yi = d1;
          case 'k'
            zi = d1;
        end
        switch l2,
          case 'i'
            xi = d2;
          case 'j'
            yi = d2;
          case 'k'
            zi = d2;
        end
      end
    end % if interactive_flag

  end % while interactive_flag

elseif isequal(cfg.method,'surface')

  % read the triangulated cortical surface from file
  tmp = load(cfg.surffile, 'bnd');
  surf = tmp.bnd;
  if isfield(surf, 'transform'),
    % compute the surface vertices in head coordinates
    surf.pnt = warp_apply(surf.transform, surf.pnt);
  end

  % downsample the cortical surface
  if cfg.surfdownsample > 1
    if ~isempty(cfg.surfinflated)
      error('downsampling the surface is not possible in combination with an inflated surface');
    end
    fprintf('downsampling surface from %d vertices\n', size(surf.pnt,1));
    [surf.tri, surf.pnt] = reducepatch(surf.tri, surf.pnt, 1/cfg.surfdownsample);
  end

  % these are required
  if ~isfield(data, 'transform'),
    data.transform = eye(4);
  end
  if ~isfield(data, 'inside')
    data.inside = true(dim);
  end

  fprintf('%d voxels in functional data\n', prod(dim));
  fprintf('%d vertices in cortical surface\n', size(surf.pnt,1));

  if hasfun 
    [interpmat, cfg.distmat] = interp_gridded(data.transform, fun, surf.pnt, 'projmethod', cfg.projmethod, 'distmat', cfg.distmat, 'sphereradius', cfg.sphereradius, 'inside', data.inside);
    % interpolate the functional data
    val = interpmat * fun(data.inside(:));
    if hasmsk
      % also interpolate the opacity mask
      maskval = interpmat * msk(data.inside(:));
    end
  end

  if ~isempty(cfg.surfinflated)
    % read the inflated triangulated cortical surface from file
    tmp = load(cfg.surfinflated, 'bnd');
    surf = tmp.bnd;
    if isfield(surf, 'transform'),
      % compute the surface vertices in head coordinates
      surf.pnt = warp_apply(surf.transform, surf.pnt);
    end
  end

  %------do the plotting
  cortex_light = [0.781 0.762 0.664];
  cortex_dark  = [0.781 0.762 0.664]/2;
  if isfield(surf, 'curv')
    % the curvature determines the color of gyri and sulci
    color = surf.curv(:) * cortex_light + (1-surf.curv(:)) * cortex_dark;
  else
    color = repmat(cortex_light, size(surf.pnt,1), 1);
  end
  h1 = patch('Vertices', surf.pnt, 'Faces', surf.tri, 'FaceVertexCData', color , 'FaceColor', 'interp');
  set(h1, 'EdgeColor', 'none');
  axis   off;
  axis vis3d;
  axis equal;

  h2 = patch('Vertices', surf.pnt, 'Faces', surf.tri, 'FaceVertexCData', val , 'FaceColor', 'interp');
  set(h2, 'EdgeColor', 'none');
  if hasmsk
    set(h2, 'FaceVertexAlphaData', maskval);
    set(h2, 'FaceAlpha',          'interp');
    set(h2, 'AlphaDataMapping',   'scaled');
    alim(gca, [opacmin opacmax]);
  end
  caxis(gca,[fcolmin fcolmax]);

  lighting gouraud
  if hasfun
    colormap(cfg.funcolormap);
  end
  if hasmsk
    alphamap(cfg.opacitymap);
  end

  if strcmp(cfg.camlight,'yes')
    camlight
  end

  if strcmp(cfg.colorbar,  'yes'),
    if hasfun
      % use a normal Matlab coorbar
      hc = colorbar;
      set(hc, 'YLim', [fcolmin fcolmax]);
    else
      warning('no colorbar possible without functional data')
    end
  end

elseif isequal(cfg.method,'slice')
  % white BG => mskana

  %% TODO: HERE THE FUNCTION THAT MAKES TO SLICE DIMENSION ALWAYS THE THIRD
  %% DIMENSION, AND ALSO KEEP TRANSFORMATION MATRIS UP TO DATE
  % zoiets
  %if hasana; ana = shiftdim(ana,cfg.slicedim-1); end;
  %if hasfun; fun = shiftdim(fun,cfg.slicedim-1); end;
  %if hasmsk; msk = shiftdim(msk,cfg.slicedim-1); end;
  %if hasmskana; mskana = shiftdim(mskana,cfg.slicedim-1); end;
  %%%%% select slices
  if ~isstr(cfg.slicerange)
    ind_fslice = cfg.slicerange(1);
    ind_lslice = cfg.slicerange(2);
  elseif isequal(cfg.slicerange, 'auto')
    if hasfun %default
      if isfield(data,'inside')
        ind_fslice = min(find(max(max(data.inside,[],1),[],2)));
        ind_lslice = max(find(max(max(data.inside,[],1),[],2)));
      else
        ind_fslice = min(find(~isnan(max(max(fun,[],1),[],2))));
        ind_lslice = max(find(~isnan(max(max(fun,[],1),[],2))));
      end
    elseif hasana %if only ana, no fun
      ind_fslice = min(find(max(max(ana,[],1),[],2)));
      ind_lslice = max(find(max(max(ana,[],1),[],2)));
    else
      error('no functional parameter and no anatomical parameter, can not plot');
    end
  else
    error('do not understand cfg.slicerange');
  end
  ind_allslice = linspace(ind_fslice,ind_lslice,cfg.nslices);
  ind_allslice = round(ind_allslice);
  % make new ana, fun, msk, mskana with only the slices that will be plotted (slice dim is always third dimension)
  if hasana; new_ana = ana(:,:,ind_allslice); clear ana; ana=new_ana; clear new_ana; end;
  if hasfun; new_fun = fun(:,:,ind_allslice); clear fun; fun=new_fun; clear new_fun; end;
  if hasmsk; new_msk = msk(:,:,ind_allslice); clear msk; msk=new_msk; clear new_msk; end;
  %if hasmskana; new_mskana = mskana(:,:,ind_allslice); clear mskana; mskana=new_mskana; clear new_mskana; end;

  % update the dimensions of the volume
  if hasana; dim=size(ana); else dim=size(fun); end;

  %%%%% make "quilts", that contain all slices on 2D patched sheet
  % Number of patches along sides of Quilt (M and N)
  % Size (in voxels) of side of patches of Quilt (m and n)
  m = dim(1);
  n = dim(2);
  M = ceil(sqrt(dim(3)));
  N = ceil(sqrt(dim(3)));
  num_patch = N*M;
  if cfg.slicedim~=3
    error('only supported for slicedim=3');
  end
  num_slice = (dim(cfg.slicedim));
  num_empt = num_patch-num_slice;
  % put empty slides on ana, fun, msk, mskana to fill Quilt up
  if hasana; ana(:,:,end+1:num_patch)=0; end;
  if hasfun; fun(:,:,end+1:num_patch)=0; end;
  if hasmsk; msk(:,:,end+1:num_patch)=0; end;
  %if hasmskana; mskana(:,:,end:num_patch)=0; end;
  % put the slices in the quilt
  for iSlice = 1:num_slice
    xbeg = floor((iSlice-1)./M);
    ybeg = mod(iSlice-1, M);
    if hasana
      quilt_ana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(ana(:,:,iSlice));
    end
    if hasfun
      quilt_fun(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(fun(:,:,iSlice));
    end
    if hasmsk
      quilt_msk(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(msk(:,:,iSlice));
    end
    %     if hasmskana
    %       quilt_mskana(ybeg.*m+1:(ybeg+1).*m, xbeg.*n+1:(xbeg+1).*n)=squeeze(mskana(:,:,iSlice));
    %     end
  end
  % make vols and scales, containes volumes to be plotted (fun, ana, msk) %added ingnie
  if hasana; vols2D{1} = quilt_ana; scales{1} = []; end; % needed when only plotting ana
  if hasfun; vols2D{2} = quilt_fun; scales{2} = [fcolmin fcolmax]; end;
  if hasmsk; vols2D{3} = quilt_msk; scales{3} = [opacmin opacmax]; end;

  plot2D(vols2D, scales);

  axis off

  if strcmp(cfg.colorbar,  'yes'),
    if hasfun
      % use a normal Matlab coorbar
      hc = colorbar;
      set(hc, 'YLim', [fcolmin fcolmax]);
    else
      warning('no colorbar possible without functional data')
    end
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% handle_ortho makes an overlay of 3D anatomical, functional and probability
% volumes. The three volumes must be scaled between 0 and 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vols2D] = handle_ortho(vols, indx, slicedir, dim);

% put 2Dvolumes in fun, ana and msk
if length(vols)>=1 && isempty(vols{1}); hasana=0; else ana=vols{1}; hasana=1; end;
if length(vols)>=2
  if isempty(vols{2}); hasfun=0; else fun=vols{2}; hasfun=1; end;
else hasfun=0; end
if length(vols)>=3
  if isempty(vols{3}); hasmsk=0; else msk=vols{3}; hasmsk=1; end;
else hasmsk=0; end

% select the indices of the intersection
xi = indx(1);
yi = indx(2);
zi = indx(3);

% select the slice to plot
if slicedir==1
  yi = 1:dim(2);
  zi = 1:dim(3);
elseif slicedir==2
  xi = 1:dim(1);
  zi = 1:dim(3);
elseif slicedir==3
  xi = 1:dim(1);
  yi = 1:dim(2);
end

% cut out the slice of interest
if hasana; ana = squeeze(ana(xi,yi,zi)); end;
if hasfun; fun = squeeze(fun(xi,yi,zi)); end;
if hasmsk; msk = squeeze(msk(xi,yi,zi)); end;

%put fun, ana and msk in vols2D
if hasana; vols2D{1} = ana; end;
if hasfun; vols2D{2} = fun; end;
if hasmsk; vols2D{3} = msk; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot2D plots a two dimensional plot, used in ortho and slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot2D(vols2D, scales);
cla;
% put 2D volumes in fun, ana and msk
hasana = length(vols2D)>0 && ~isempty(vols2D{1});
hasfun = length(vols2D)>1 && ~isempty(vols2D{2});
hasmsk = length(vols2D)>2 && ~isempty(vols2D{3});

% the transpose is needed for displaying the matrix using the Matlab image() function
if hasana; ana = vols2D{1}'; end;
if hasfun; fun = vols2D{2}'; end;
if hasmsk; msk = vols2D{3}'; end;


if hasana
  % scale anatomy between 0 and 1
  fprintf('scaling anatomy\n');
  amin = min(ana(:));
  amax = max(ana(:));
  ana = (ana-amin)./(amax-amin);
  clear amin amax;
  % convert anatomy into RGB values
  ana = cat(3, ana, ana, ana);
  ha = imagesc(ana);
end
hold on

if hasfun
  hf = imagesc(fun);
  caxis(scales{2});
  % apply the opacity mask to the functional data
  if hasmsk
    % set the opacity
    set(hf, 'AlphaData', msk)
    set(hf, 'AlphaDataMapping', 'scaled')
    alim(scales{3});
  elseif hasana
    set(hf, 'AlphaData', 0.5)
  end
end

axis equal
axis tight
axis xy

