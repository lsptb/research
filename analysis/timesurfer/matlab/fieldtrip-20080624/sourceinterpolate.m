function [interp] = sourceinterpolate(cfg, functional, anatomical);

% SOURCEINTERPOLATE reslices and interpolates a source reconstruction
% or a statistical distribution as an overlay onto an anatomical MRI.
%
% The source volume and the anatomical volume should be expressed in the
% same coordinate sytem, i.e. either both in CTF coordinates (NAS/LPA/RPA)
% or both in SPM coordinates (AC/PC). The output volume will contain a
% resliced source and anatomical volume that can be plotted together with
% SOURCEPLOT or SLICEINTERP, or that can be written to file using SOURCEWRITE.
%
% Use as
%   [interp] = sourceinterpolate(cfg, source, mri)   or
%   [interp] = sourceinterpolate(cfg, stat, mri)
% where
%   source is the output of SOURCEANALYSIS
%   stat   is the output of SOURCESTATISTICS
%   mri    is the output of READ_FCDC_MRI or the filename of a MRI
% and cfg is a structure with any of the following fields
%   cfg.parameter     = string, default is 'all'
%   cfg.interpmethod  = 'linear', 'cubic', 'nearest' or 'spline'
%   cfg.sourceunits   = 'mm' or 'cm' (default is 'cm')
%   cfg.mriunits      = 'mm' or 'cm' (default is 'mm')
%   cfg.downsample    = integer number (default = 1, i.e. no downsampling)
%
% See also SOURCEANALYSIS, SOURCESTATISTICS, READ_FCDC_MRI

% Undocumented options
%   cfg.keepinside = 'yes' (default) or 'no'
%   cfg.voxelcoord = 'yes' (default) or 'no' determines whether the
%   downsampled output anatomical MRI will have the x/y/zgrid converted or
%   the homogenous transformation matrix

% Copyright (C) 2003-2007, Robert Oostenveld
%
% $Log: sourceinterpolate.m,v $
% Revision 1.43  2007/04/18 10:30:41  roboos
% switched to interpolation in voxel indeices, since in some cases interpn had problems (esp when the coordinate axes were strongly rotated)
% apply the cm->mm scaling to the functional transfortmation matrix prior to computing voxel positions
% allow the interpolation of one anatomy onto another, i.e. do not overwrite anatomy if it was interpolated
%
% Revision 1.42  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.41  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.40  2007/02/27 15:24:57  ingnie
% sel from voxel index to logical index, ~sel volumes NaNs iso zeros (except inside)
%
% Revision 1.39  2006/07/27 08:30:01  roboos
% updated documentation
%
% Revision 1.38  2006/06/08 07:50:52  roboos
% updated the conversion between the source and MRI units (support mm,cm,dm,m for both)
%
% Revision 1.37  2006/03/09 08:19:44  roboos
% try not to interpolate cell-array volumes, e.g. mom, lf, csd
%
% Revision 1.36  2006/02/24 16:44:41  roboos
% switched to fixvolume function, with inside as 3d volume
% not neccessary any more to treat inside/outside seperately
% not neccessary any more to reshape volumes
% not neccessary any more to change from xgrid/etc into transform
%
% Revision 1.35  2006/02/07 22:21:47  roboos
% changed fprintf feedback into progress() function, added cfg.feedback
%
% Revision 1.34  2006/01/31 09:48:39  jansch
% fixed bug regarding last change
%
% Revision 1.33  2006/01/31 09:45:32  jansch
% changed else cfg.keepinside='no' into an elseif to prevent crashes
%
% Revision 1.32  2006/01/31 09:27:58  roboos
% in case of keepinside: instead of setting outside to [], remove it from structure
%
% Revision 1.31  2006/01/30 13:55:19  roboos
% changed the order in which the initial preparatory steps are taken
% related to inside-conversion, parameterselection and grid2transform
%
% Revision 1.30  2006/01/24 21:28:35  roboos
% removes obsolete tmpcfg.voxelcoord for downsample (default is now yes)
% updated the local voxelcoords() subfunction
%
% Revision 1.29  2006/01/05 13:38:23  roboos
% Changed the use of the parameterselection function which now always returns a cell-array.
% Use the private/grid2transform function to ensure that both the
% anatomical and functional volumes are described using the homogenous
% transformation matrix and not with xgrid/ygrid/zgrid.
%
% Revision 1.28  2005/11/10 10:02:14  roboos
% minor change in whitespace and help
%
% Revision 1.27  2005/10/14 15:48:57  roboos
% add inside to cfg.parameter if keepinside=yes and parameter~=all
%
% Revision 1.26  2005/08/29 10:18:16  roboos
% perform downsamplevolume with restricted copy of configuration and explicitely on parameter=anatomy
%
% Revision 1.25  2005/08/19 17:10:35  roboos
% completely new implementation from scratch. It now also supports interpolation of
% any volume to any volume, i.e. the function input does not have to be aligned with
% the axes of the coordinate system.
%
% Revision 1.24  2005/06/29 12:46:29  roboos
% the following changes were done on a large number of files at the same time, not all of them may apply to this specific file
% - added try-catch around the inclusion of input data configuration
% - moved cfg.version, cfg.previous and the assignment of the output cfg to the end
% - changed help comments around the configuration handling
% - some changes in whitespace
%
% Revision 1.23  2005/06/10 09:13:58  jansch
% fixed bug which was introduced during the last revision
%
% Revision 1.22  2005/06/03 08:59:59  roboos
% added support for interpolating already interpolated low-res stats structures back onto the high-res anatomy
% fixed bug in scaling of source positions (it always scaled them 10x, regardless of cfg.sourceunits)
%
% Revision 1.21  2005/05/17 17:50:38  roboos
% changed all "if" occurences of & and | into && and ||
% this makes the code more compatible with Octave and also seems to be in closer correspondence with Matlab documentation on shortcircuited evaluation of sequential boolean constructs
%
% Revision 1.20  2005/02/09 13:09:52  roboos
% explicitly emptied source.outside if keepinside=yes
%
% Revision 1.19  2005/02/08 11:45:00  roboos
% implemented cfg.mriunits (default mm)
% implemented cfg.keepinside (default yes)
% removed obsolete cfg.parameter
% changed handling of all volumes that should be interpolated (now through list with subfields)
% cleaned up code and help
%
% Revision 1.18  2004/10/13 14:11:30  roboos
% changed cfg.previous, now consistent over functions with try-statement
% and will also work if the input already had a cfg.previous
%
% Revision 1.17  2004/09/21 09:37:24  jansch
% fixed small bug due to yesterday's changes
%
% Revision 1.16  2004/09/20 09:32:15  roboos
% changed implementation of downsampling to use the external downsamplevolume function
%

%% checkdata see below!!! %%

% set the defaults
if ~isfield(cfg, 'parameter'),    cfg.parameter    = 'all';     end
if ~isfield(cfg, 'interpmethod'); cfg.interpmethod = 'linear';  end
if ~isfield(cfg, 'downsample');   cfg.downsample   = 1;         end
if ~isfield(cfg, 'sourceunits');  cfg.sourceunits  = 'cm';      end
if ~isfield(cfg, 'mriunits');     cfg.mriunits     = 'mm';      end
if ~isfield(cfg, 'voxelcoord'),   cfg.voxelcoord   = 'yes';     end
if ~isfield(cfg, 'keepinside'),   cfg.keepinside   = 'yes';     end
if ~isfield(cfg, 'feedback'),     cfg.feedback     = 'text';    end

if strcmp(cfg.keepinside, 'yes')
  % add inside to the list of parameters
  if ~iscell(cfg.parameter),
    cfg.parameter = {cfg.parameter 'inside'};
  else
    cfg.parameter(end+1) = {'inside'};
  end
end

if ischar(anatomical)
  % read the anatomical MRI data from file
  fprintf('reading MRI from file\n');
  filename   = anatomical;
  anatomical = read_fcdc_mri(filename);
end

% ensure that the structure correctly describes a volume
functional = fixvolume(functional);
anatomical = fixvolume(anatomical);

% check if the input data is valid for this function
functional = checkdata(functional, 'datatype', 'source', 'feedback', 'yes');
anatomical = checkdata(anatomical, 'datatype', 'volume', 'feedback', 'yes');

% select the parameters that should be interpolated
cfg.parameter = parameterselection(cfg.parameter, functional);

% downsample the anatomical volume
tmpcfg = [];
tmpcfg.downsample = cfg.downsample;
tmpcfg.parameter  = 'anatomy';
anatomical = volumedownsample(tmpcfg, anatomical);

% collect the functional volumes that should be converted
vol_name = {};
vol_data = {};
for i=1:length(cfg.parameter)
  if ~iscell(getsubfield(functional, cfg.parameter{i}))
    vol_name{end+1} = cfg.parameter{i};
    vol_data{end+1} = getsubfield(functional, cfg.parameter{i});
  else
    fprintf('not interpolating %s, since it is not a scalar field\n', cfg.parameter{i});
  end
end

% convert the source/functional data into the same units as the anatomical MRI
s = 1;
switch cfg.sourceunits
case 'mm'
  s = s / 1000;
case 'cm'
  s = s / 100;
case 'dm'
  s = s / 10;
case 'm'
  s = s / 1;
otherwise
  error('unknown physical dimension in cfg.sourceunits');
end
switch cfg.mriunits
case 'mm'
  s = s * 1000;
case 'cm'
  s = s * 100;
case 'dm'
  s = s * 10;
case 'm'
  s = s * 1;
otherwise
  error('unknown physical dimension in cfg.mriunits');
end

if s~=1
  fprintf('converting functional data from %s into %s\n', cfg.sourceunits, cfg.mriunits);
  functional.transform = scale([s s s]) * functional.transform;
end

% compute the position of each voxel in both volumes, expressed in headcoordinates
[fx, fy, fz] = voxelcoords(functional);
[ax, ay, az] = voxelcoords(anatomical);
% convert the anatomical voxel positions into voxel indices into the functional volume
pos = [ax(:) ay(:) az(:)];
pos = warp_apply(inv(functional.transform), pos);
ax = reshape(pos(:,1), anatomical.dim);
ay = reshape(pos(:,2), anatomical.dim);
az = reshape(pos(:,3), anatomical.dim);
clear pos

% estimate the subvolume of the anatomy that is spanned by the functional volume
minfx = 1;
minfy = 1;
minfz = 1;
maxfx = functional.dim(1);
maxfy = functional.dim(2);
maxfz = functional.dim(3);
sel = ax(:)>=minfx & ...
      ax(:)<=maxfx & ...
      ay(:)>=minfy & ...
      ay(:)<=maxfy & ...
      az(:)>=minfz & ...
      az(:)<=maxfz;
fprintf('selecting subvolume of %.1f%%\n', 100*sum(sel)./prod(anatomical.dim));

% start with an empty output structure
interp = [];

% reslice and interpolate all functional volumes
for i=1:length(vol_name)
  fprintf('reslicing and interpolating %s\n', vol_name{i});
  fv = double(vol_data{i});
  av = zeros(anatomical.dim);
  av( sel) = my_interpn(fv, ax(sel), ay(sel), az(sel), cfg.interpmethod, cfg.feedback);
  % av( sel) = my_interpn(fx, fy, fz, fv, ax(sel), ay(sel), az(sel), cfg.interpmethod, cfg.feedback);
  av(~sel) = nan;
  interp = setsubfield(interp, vol_name{i}, av);
end

if isfield(interp, 'inside')
  % convert back to a logical volume
  interp.inside(~sel) = 0; % these values were previously set to NaN (see 10 lines above)
  interp.inside       = abs(interp.inside-1)<=10*eps;
end

% add the other parameters to the output
interp.dim       = anatomical.dim;
interp.transform = anatomical.transform;
if ~any(strcmp(cfg.parameter, 'anatomy'))
  % copy the anatomy into the functional data
  interp.anatomy   = anatomical.anatomy;
end

% add version information to the configuration
try
  % get the full name of the function
  cfg.version.name = mfilename('fullpath');
catch
  % required for compatibility with Matlab versions prior to release 13 (6.5)
  [st, i] = dbstack;
  cfg.version.name = st(i);
end
cfg.version.id = '$Id: sourceinterpolate.m,v 1.43 2007/04/18 10:30:41 roboos Exp $';
% remember the configuration details of the input data
cfg.previous = [];
try, cfg.previous{1} = functional.cfg; end
try, cfg.previous{2} = anatomical.cfg; end
% remember the exact configuration details in the output
interp.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION this function computes the location of all voxels in head
% coordinates in a memory efficient manner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, z] = voxelcoords(volume)
dim       = volume.dim;
transform = volume.transform;
if isfield(volume, 'xgrid')
  xgrid = volume.xgrid;
  ygrid = volume.ygrid;
  zgrid = volume.zgrid;
else
  xgrid = 1:dim(1);
  ygrid = 1:dim(2);
  zgrid = 1:dim(3);
end
npix = prod(dim(1:2));  % number of voxels in a single slice
nvox = prod(dim);
x = zeros(dim);
y = zeros(dim);
z = zeros(dim);
X = zeros(1,npix);
Y = zeros(1,npix);
Z = zeros(1,npix);
E = ones(1,npix);
% determine the voxel locations per slice
for i=1:dim(3)
  [X(:), Y(:), Z(:)] = ndgrid(xgrid, ygrid, zgrid(i));
  tmp = transform*[X; Y; Z; E];
  x((1:npix)+(i-1)*npix) = tmp(1,:);
  y((1:npix)+(i-1)*npix) = tmp(2,:);
  z((1:npix)+(i-1)*npix) = tmp(3,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBFUNCTION for memory efficient interpolation
% the only reason for this function is that it does the interpolation in smaller chuncks
% this prevents memory problems that I often encountered here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [av] = my_interpn(fx, fy, fz, fv, ax, ay, az, interpmethod, feedback);
function [av] = my_interpn(fv, ax, ay, az, interpmethod, feedback);
num = prod(size(ax));       % total number of voxels
blocksize = floor(num/20);  % number of voxels to interpolate at once, split it into 20 chuncks
lastblock = 0;              % boolean flag for while loop
sel = 1:blocksize;          % selection of voxels that are interpolated, this is the first chunck
progress('init', feedback, 'interpolating');
while (1)
  progress(sel(1)/num, 'interpolating %.1f%%\n', 100*sel(1)/num);
  if sel(end)>num
    sel = sel(1):num;
    lastblock = 1;
  end
  av(sel) = interpn(fv, ax(sel), ay(sel), az(sel), interpmethod);
  if lastblock
    break
  end
  sel = sel + blocksize;
end
progress('close');

