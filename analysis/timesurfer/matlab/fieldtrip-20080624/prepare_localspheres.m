function [vol, cfg] = prepare_localspheres(cfg, mri);

% PREPARE_LOCALSPHERES creates a MEG volume conductor model with a sphere
% for every sensor. You can also use it to create a single sphere
% model that is fitted to the MRI or to the head shape points.
%
% Use as
%   [vol, cfg] = prepare_localspheres(cfg)
%   [vol, cfg] = prepare_localspheres(cfg, seg)
%
% The input configuration should contain
%   cfg.grad         = structure with gradiometer definition, or
%   cfg.gradfile     = filename containing gradiometer definition
%   cfg.headshape    = filename containing headshape, or Nx3 matrix with surface points
%   cfg.radius       = number, which points to select for each channel (default = 7 cm)
%   cfg.baseline     = number, baseline of axial/planar gradiometer (default = 5 cm)
%   cfg.feedback     = 'yes' or 'no' (default = 'yes')
%   cfg.singlesphere = 'yes' or 'no', fit only a single sphere (default = 'no')
%
% The following options are relevant if you use a segmented MRI
%   cfg.smooth      = 'no' or the FWHM of the gaussian kernel in voxels (default = 'no')
%   cfg.mriunits    = 'mm' or 'cm' (default = 'mm')
%   cfg.sourceunits = 'mm' or 'cm' (default = 'cm')
%   cfg.threshold   = 0.5, relative to the maximum value in the segmentation
%
% This function implements
%   Huang MX, Mosher JC, Leahy RM.
%   A sensor-weighted overlapping-sphere head model and exhaustive head model comparison for MEG
%   Phys Med Biol. 1999 Feb;44(2):423-40

% TODO cfg.spheremesh  should be renamed consistently with other mesh generation cfgs
%
% Undocumented local options:
% cfg.spheremesh, number of points that is placed on the brain surface (default 4000)
% cfg.maxradius

% Copyright (C) 2005-2006, Jan-Mathijs Schoffelen & Robert Oostenveld
%
% $Log: prepare_localspheres.m,v $
% Revision 1.18  2008/05/13 15:37:24  roboos
% switched to using read_data/header instead of the read_fcdc_data/header wrapper functions
%
% Revision 1.17  2007/08/06 09:20:14  roboos
% added support for bti_hs
%
% Revision 1.16  2007/07/26 08:02:10  roboos
% remove double vertices in headshape
% changed some indentation
%
% Revision 1.15  2006/10/09 15:25:00  roboos
% added option to fit only a single sphere to the MRI or headshape points
%
% Revision 1.14  2006/08/16 10:52:12  marsie
% fixed use of spm_smooth()
%
% Revision 1.13  2006/08/01 10:32:01  marsie
% fixed bug in using spm_smooth
%
% Revision 1.12  2006/07/27 08:29:38  roboos
% use spm_smooth instead of spm_conv, updated documentation
%
% Revision 1.11  2006/06/08 07:50:52  roboos
% updated the conversion between the source and MRI units (support mm,cm,dm,m for both)
%
% Revision 1.10  2006/06/07 09:34:19  roboos
% changed checktoolbox into hastoolbox
%
% Revision 1.9  2006/04/20 09:58:34  roboos
% updated documentation
%
% Revision 1.8  2006/04/18 19:04:35  roboos
% changed hard-coded 10000 vertex points for brain surface into cfg.spheremesh with default value of 4000
%
% Revision 1.7  2006/04/05 15:07:37  roboos
% return the configuration as second argument, store the headshape points in the cfg
%
% Revision 1.6  2006/04/03 10:47:13  roboos
% pre-allocate memory to hold the result, ensure that the dimensions are correct (thanks to Hanneke)
%
% Revision 1.5  2006/03/21 09:45:07  roboos
% some irrelevant changes in the code flow to deal with presence/absence of segmented MRI input
%
% Revision 1.4  2006/03/09 08:18:34  roboos
% added code to use a segmented anatomical MRI, i.e. smooth and do edge detection
% with options similar to prepare_dipole_grid in gray matter
% updated help
%
% Revision 1.3  2006/02/06 21:22:28  roboos
% removed obvious dependency on CTF for reading hdr.grad
%
% Revision 1.2  2006/01/25 10:08:27  roboos
% added the complete implementation of this function, sofar it had been empty
% the implementation is based on code from jansch, which again was partially based on roboos' code
%
% Revision 1.1  2005/11/03 11:15:45  roboos
% new implementation
%

% set the defaults
if ~isfield(cfg, 'radius'),        cfg.radius = 8.5;        end
if ~isfield(cfg, 'maxradius'),     cfg.maxradius = 20;      end
if ~isfield(cfg, 'baseline'),      cfg.baseline = 5;        end
if ~isfield(cfg, 'feedback'),      cfg.feedback = 'yes';    end
if ~isfield(cfg, 'smooth');        cfg.smooth    = 5;       end % in voxels
if ~isfield(cfg, 'mriunits');      cfg.mriunits = 'mm';     end
if ~isfield(cfg, 'sourceunits'),   cfg.sourceunits = 'cm';  end
if ~isfield(cfg, 'threshold'),     cfg.threshold = 0.5;     end % relative
if ~isfield(cfg, 'spheremesh'),    cfg.spheremesh = 4000;   end
if ~isfield(cfg, 'singlesphere'),  cfg.singlesphere = 'no'; end

if nargin>1
  % obtain the head shape from the segmented MRI
  if isfield(cfg, 'headshape')
    error('cfg.headshape should not be used in combination with a segmented mri')
  end
  seg = zeros(mri.dim);
  if isfield(mri, 'gray')
    fprintf('including gray matter in segmentation for brain compartment\n')
    seg = seg | (mri.gray>(cfg.threshold*max(mri.gray(:))));
  end
  if isfield(mri, 'white')
    fprintf('including white matter in segmentation for brain compartment\n')
    seg = seg | (mri.white>(cfg.threshold*max(mri.white(:))));
  end
  if isfield(mri, 'csf')
    fprintf('including CSF in segmentation for brain compartment\n')
    seg = seg | (mri.csf>(cfg.threshold*max(mri.csf(:))));
  end
  if ~strcmp(cfg.smooth, 'no'),
    % check whether the required SPM2 toolbox is available
    hastoolbox('spm2', 1);
    fprintf('smoothing the segmentation with a %d-pixel FWHM kernel\n',cfg.smooth);
    seg = double(seg);
    spm_smooth(seg, seg, cfg.smooth);
  end
  % threshold for the last time
  seg = (seg>(cfg.threshold*max(seg(:))));
  % determine the center of gravity of the segmented brain
  xgrid = 1:mri.dim(1);
  ygrid = 1:mri.dim(2);
  zgrid = 1:mri.dim(3);
  [X, Y, Z] = ndgrid(xgrid, ygrid, zgrid);
  ori(1) = mean(X(seg));
  ori(2) = mean(Y(seg));
  ori(3) = mean(Z(seg));
  pnt = triangulate_seg(seg, cfg.spheremesh, ori);
  pnt(:,4) = 1;
  pnt = (mri.transform * pnt')';
  % convert the MRI surface points into the same units as the source/gradiometer
  scale = 1;
  switch cfg.sourceunits
    case 'mm'
      scale = scale * 1000;
    case 'cm'
      scale = scale * 100;
    case 'dm'
      scale = scale * 10;
    case 'm'
      scale = scale * 1;
    otherwise
      error('unknown physical dimension in cfg.sourceunits');
  end
  switch cfg.mriunits
    case 'mm'
      scale = scale / 1000;
    case 'cm'
      scale = scale / 100;
    case 'dm'
      scale = scale / 10;
    case 'm'
      scale = scale / 1;
    otherwise
      error('unknown physical dimension in cfg.mriunits');
  end
  if scale~=1
    fprintf('converting MRI surface points from %s into %s\n', cfg.sourceunits, cfg.mriunits);
    shape = pnt(:,1:3) * scale;
  else
    shape = pnt(:,1:3);
  end
  fprintf('placed %d points on the brain surface\n', length(shape));
elseif ischar(cfg.headshape) && filetype(cfg.headshape, 'ctf_shape')
  % read the headshape from file
  shape = read_ctf_shape(cfg.headshape);
  shape = shape.pnt;
elseif ischar(cfg.headshape) && filetype(cfg.headshape, '4d_hs')
  % read the headshape from file
  shape = read_bti_hs(cfg.headshape);
else
  % use the headshape points that are specified in the configuration
  shape = cfg.headshape;
end % nargin

% remove double vertices
shape = unique(shape, 'rows');

% remember the headshape in the configuration, mainly for debugging
cfg.headshape = shape;

% read the gradiometer definition from file or copy it from the configuration
if isfield(cfg, 'gradfile')
  hdr = read_header(cfg.gradfile);
  grad = hdr.grad;
else
  grad = cfg.grad;
end

Nshape = size(shape,1);
Nchan  = size(grad.tra, 1);

% set up an empty figure
if strcmp(cfg.feedback, 'yes')
  figure
  hold on
  axis equal
  axis vis3d
  axis off
  drawnow
end

% plot all channels and headshape points
if strcmp(cfg.feedback, 'yes')
  cla
  plot3(grad.pnt(:,1), grad.pnt(:,2), grad.pnt(:,3), 'b.');	% all coils
  plot3(shape(:,1), shape(:,2), shape(:,3), 'g.');
  drawnow
end

% fit a single sphere to all headshape points
ini = mean(shape,1);
[single_o, single_r] = fit_sphere(shape, ini);
fprintf('single sphere,   %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', Nshape, single_o(1), single_o(2), single_o(3), single_r);

if strcmp(cfg.singlesphere, 'yes')
  % only return a single sphere
  vol.r = single_r;
  vol.o = single_o;
  return;
end

% start with an empty structure that will hold the results
vol.r = zeros(Nchan,1);    % radius of every sphere
vol.o = zeros(Nchan,3);    % origin of every sphere
vol.label = cell(Nchan,1); % corresponding gradiometer channel label for every sphere

for chan=1:Nchan
  coilsel = find(grad.tra(chan,:)~=0);
  allpnt  = grad.pnt(coilsel, :);	% position of all coils belonging to this channel
  allori  = grad.ori(coilsel, :);	% orientation of all coils belonging to this channel

  if strcmp(cfg.feedback, 'yes')
    cla
    plot3(grad.pnt(:,1), grad.pnt(:,2), grad.pnt(:,3), 'b.');	% all coils
    plot3(allpnt(:,1), allpnt(:,2), allpnt(:,3), 'r*');		% this channel in red
  end

  % determine the average position and orientation of this channel
  thispnt = mean(allpnt,1);
  [u, s, v] = svd(allori);
  thisori = v(:,1)';
  if dot(thispnt,thisori)<0
    % the orientation should be outwards pointing
    thisori = -thisori;
  end

  % compute the distance from every coil along this channels orientation
  dist = [];
  for i=1:length(coilsel)
    dist(i) = dot((allpnt(i,:)-thispnt), thisori);
  end

  [m, i] = min(dist);
  % check whether the minimum difference is larger than a typical distance
  if abs(m)>(cfg.baseline/4)
    % replace the position of this channel by the coil that is the closest to the head (axial gradiometer)
    % except when the center of the channel is approximately just as good (planar gradiometer)
    thispnt = allpnt(i,:);
  end

  % find the headshape points that are close to this channel
  dist = sqrt(sum((shape-repmat(thispnt,Nshape,1)).^2, 2));
  shapesel = find(dist<cfg.radius);
  if strcmp(cfg.feedback, 'yes')
    plot3(shape(shapesel,1), shape(shapesel,2), shape(shapesel,3), 'g.');
    drawnow
  end

  % fit a sphere to these headshape points
  if length(shapesel)>10
    ini = mean(shape(shapesel,:),1) - cfg.radius * thisori;
    [o, r] = fit_sphere(shape(shapesel,:), ini);
    fprintf('channel = %s, %5d surface points, center = [%4.1f %4.1f %4.1f], radius = %4.1f\n', grad.label{chan}, length(shapesel), o(1), o(2), o(3), r);
  else
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end

  if r > cfg.maxradius
    fprintf('channel = %s, not enough surface points, using all points\n', grad.label{chan});
    o = single_o;
    r = single_r;
  end

  % add this sphere to the volume conductor
  vol.o(chan,:)   = o;
  vol.r(chan)     = r;
  vol.label{chan} = grad.label{chan};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function that fits a sphere to a cloud of points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [center, radius] = fit_sphere(pnt, initial);
% start the dipole fitting with default options
options = optimset('TolFun',1e-10,...
  'TypicalX',norm(range(pnt))/100,...
  'Display','none');
warning off
center = fminunc(@fit_sphere_error, initial, options, pnt);
warning on
[err, radius] = fit_sphere_error(center, pnt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% helper function that computes the goal function to be minimized
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [err, radius] = fit_sphere_error(center, pnt);
dist   = sqrt(sum(((pnt - repmat(center, size(pnt,1), 1)).^2), 2));
radius = mean(dist);
err    = std(dist,1);

