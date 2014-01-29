function [normalise] = volumenormalise(cfg, interp)

% VOLUMENORMALISE normalises anatomical and functional volume data
% to a template anatomical MRI.
%
% Use as
%   [volume] = volumenormalise(cfg, volume)
%
% The input volume should be the result from SOURCEINTERPOLATE.
% Alternatively, the input can contain a single anatomical MRI that
% was read with READ_FCDC_MRI, or you can specify a filename of an
% anatomical MRI.
%
% Configuration options are:
%   cfg.template    = string with filename (default = '/home/common/matlab/spm2/templates/T1.mnc')
%   cfg.parameter   = cell-array with the functional data which has to
%                     be normalised, can be 'all'
%   cfg.downsample  = integer number (default = 1, i.e. no downsampling)
%   cfg.coordinates = 'spm, 'ctf' or empty for interactive (default = [])
%   cfg.name        = string for output filename
%   cfg.write       = 'no' (default) or 'yes', writes the segmented volumes to SPM2
%                     compatible analyze-file, with the suffix
%                     _anatomy for the anatomical MRI volume
%                     _param   for each of the functional volumes
%   cfg.nonlinear   = 'yes' (default) or 'no', estimates a nonlinear transformation
%                     in addition to the linear affine registration. If a reasonably
%                     accurate normalisation is sufficient, a purely linearly transformed
%                     image allows for 'reverse-normalisation', which might come in handy
%                     when for example a region of interest is defined on the normalised
%                     group-average.

% Undocumented local options:
%   cfg.keepintermediate = 'yes' or 'no'
%   cfg.intermediatename = prefix of the the coregistered images and of the
%                          original images in the original headcoordinate system
%   cfg.spmparams        = one can feed in parameters from a prior normalisation

% Copyright (C) 2004-2006, Jan-Mathijs Schoffelen
%
% $Log: volumenormalise.m,v $
% Revision 1.15  2007/04/03 15:37:07  roboos
% renamed the checkinput function to checkdata
%
% Revision 1.14  2007/03/30 17:05:40  ingnie
% checkinput; only proceed when input data is allowed datatype
%
% Revision 1.13  2006/07/27 08:29:09  roboos
% updated documentation
%
% Revision 1.12  2006/07/25 10:07:06  roboos
% Changed the method to check for spm2, add automatically to path if possible.
% Moved the interactive questions regarding coordinate system to a subfunction.
% Not only determine input source coordinates, but also template coordinates.
% Start with a well defined "initial" transform, and end with a "final" transform, both are stored in the output cfg.
% Use a hard-coded initial transform for ctf2spm instead of align_ctf2spm function.
% Renamed the internal variable names for the volumes for consistency with SPM.
% Removed old log comments.
%
% Revision 1.11  2006/07/19 12:22:52  jansch
% implemented smoothing for functional volumes

%% checkdata see below!!! %%

% check the availability of the required SPM2 toolbox
hastoolbox('spm2', 1);

% set the defaults
if ~isfield(cfg,'template'),         cfg.template = '/home/common/matlab/spm2/templates/T1.mnc'; end;
if ~isfield(cfg,'parameter'),        cfg.parameter = 'all';                       end;
if ~isfield(cfg,'downsample'),       cfg.downsample = 1;                          end;
if ~isfield(cfg,'write'),            cfg.write = 'no';                            end;
if ~isfield(cfg,'keepinside'),       cfg.keepinside = 'yes';                      end;
if ~isfield(cfg,'keepintermediate'), cfg.keepintermediate = 'no';                 end;
if ~isfield(cfg,'coordinates'),      cfg.coordinates = [];                        end;
if ~isfield(cfg,'initial'),          cfg.initial = [];                            end;
if ~isfield(cfg,'nonlinear'),        cfg.nonlinear = 'yes';                       end;
if ~isfield(cfg,'smooth'),           cfg.smooth    = 'no';                        end;

if strcmp(cfg.keepinside, 'yes')
  % add inside to the list of parameters
  if ~iscell(cfg.parameter),
    cfg.parameter = {cfg.parameter 'inside'};
  else
    cfg.parameter(end+1) = {'inside'};
  end
end

if ~isfield(cfg,'intermediatename')
  cfg.intermediatename = tempname;
end

if ~isfield(cfg,'name') && strcmp(cfg.write,'yes')
  error('you must specify the output filename in cfg.name');
end

if isempty(cfg.template),
  error('you must specify a template anatomical MRI');
end

if ~isfield(interp,'anatomy'),
  error('no anatomical information available, this is required for normalisation');
end

if isfield(cfg, 'coordinates')
  source_coordinates = cfg.coordinates;
end

% the template anatomy should always be stored in a SPM-compatible file
if filetype(cfg.template, 'analyze_hdr') || filetype(cfg.template, 'analyze_img') || filetype(cfg.template, 'minc')
  % based on the filetype assume that the coordinates correspond with MNI/SPM convention
  template_coordinates = 'spm';
end

% the source anatomy can be in a file that is not understood by SPM (e.g. in the
% native CTF *.mri format), therefore start by reading it into memory
if ischar(interp),
  fprintf('reading source MRI from file\n');
  filename = interp;
  interp   = read_fcdc_mri(filename);
  if filetype(filename, 'ctf_mri')
    % based on the filetype assume that the coordinates correspond with CTF convention
    source_coordinates = 'ctf';
  end
end

% ensure that the structure correctly describes a volume
interp = fixvolume(interp);

% check if the input data is valid for this function
interp = checkdata(interp, 'datatype', 'volume', 'feedback', 'yes');

% select the parameters that should be normalised
cfg.parameter = parameterselection(cfg.parameter, interp);

% the anatomy should always be normalised as the first volume
cfg.parameter(strcmp(cfg.parameter, 'anatomy')) = [];
cfg.parameter = {'anatomy' cfg.parameter{:}};

% downsample the volume
tmpcfg            = [];
tmpcfg.downsample = cfg.downsample;
tmpcfg.parameter  = cfg.parameter;
tmpcfg.smooth     = cfg.smooth;
interp = volumedownsample(tmpcfg, interp);

if isempty(source_coordinates)
  source_coordinates = determine_coordinates('input');       % use interactive helper-function
end
if isempty(template_coordinates)
  template_coordinates = determine_coordinates('template');  % use interactive helper-function
end

% ensure that the source volume is approximately aligned with the template
if strcmp(source_coordinates, 'ctf') && strcmp(template_coordinates, 'spm')
  fprintf('converting input coordinates from CTF into approximate SPM coordinates\n');
  initial = [
    0.0000   -1.0000    0.0000    0.0000
    0.9987    0.0000   -0.0517  -32.0000
    0.0517    0.0000    0.9987  -54.0000
    0.0000    0.0000    0.0000    1.0000 ];
elseif strcmp(source_coordinates, template_coordinates)
  fprintf('not converting input coordinates\n');
  initial = eye(4);
else
  error('cannot determine the approximate alignmenmt of the source coordinates with the template');
end

% apply the approximate transformatino prior to passing it to SPM
interp.transform = initial * interp.transform;

warning off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here the normalisation starts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create an spm-compatible header for the anatomical volume data
VF = volumewrite_spm([cfg.intermediatename,'_anatomy.img'], interp.anatomy, interp.transform);

% create an spm-compatible file for each of the functional volumes
for parlop=2:length(cfg.parameter)  % skip the anatomy
  tmp  = cfg.parameter{parlop};
  data = reshape(getsubfield(interp, tmp), interp.dim);
  tmp(find(tmp=='.')) = '_';
  volumewrite_spm([cfg.intermediatename,'_' tmp '.img'], data, interp.transform);
end

% read the template anatomical volume
if strcmp(char(cfg.template(end-2:end)),'mnc'),
  VG    = spm_vol_minc(cfg.template);
elseif strcmp(char(cfg.template(end-2:end)),'img'),
  VG    = spm_vol(cfg.template);
else
  error('Unknown template');
end

fprintf('performing the normalisation\n');
% do spatial normalisation according to these steps
% step 1: read header information for template and source image
% step 2: compute transformation parameters
% step 3: write the results to a file with prefix 'w'

if ~isfield(cfg, 'spmparams') && strcmp(cfg.nonlinear, 'yes'),
  fprintf('warping the invdividual anatomy to the template anatomy\n');
  % compute the parameters by warping the individual anatomy
  VF        = spm_vol([cfg.intermediatename,'_anatomy.img']);
  params    = spm_normalise(VG,VF);
elseif ~isfield(cfg, 'spmparams') && strcmp(cfg.nonlinear, 'no'),
  fprintf('warping the invdividual anatomy to the template anatomy, using only linear transformations\n');
  % compute the parameters by warping the individual anatomy
  VF         = spm_vol([cfg.intermediatename,'_anatomy.img']);
  flags.nits = 0; %put number of non-linear iterations to zero
  params     = spm_normalise(VG,VF,[],[],[],flags);
else
  fprintf('using the parameters specified in the configuration\n');
  % use the externally specified parameters
  params = cfg.spmparams;
end
flags.vox = [cfg.downsample,cfg.downsample,cfg.downsample];
files     = {};

% determine the affine source->template coordinate transformation
final = VG.mat * inv(params.Affine) * inv(VF.mat) * initial;

% remember the normalisation parameters in the configuration
cfg.spmparams = params;
cfg.initial   = initial;
cfg.final     = final;

% apply the normalisation parameters to each of the volumes
for parlop=1:length(cfg.parameter)
  fprintf('creating normalised analyze-file for %s\n', cfg.parameter{parlop});
  tmp = cfg.parameter{parlop};
  tmp(find(tmp=='.')) = '_';
  files{parlop} = sprintf('%s_%s.img', cfg.intermediatename, tmp);
  [p, f, x] = fileparts(files{parlop});
  wfiles{parlop} = fullfile(p, ['w' f x]);
end
spm_write_sn(char(files),params,flags);  % this creates the 'w' prefixed files

normalise = [];

% read the normalised results from the 'w' prefixed files
V = spm_vol(char(wfiles));
for vlop=1:length(V)
  normalise = setsubfield(normalise, cfg.parameter{vlop}, spm_read_vols(V(vlop)));
end

normalise.transform = V(1).mat;
normalise.dim       = size(normalise.anatomy);

if isfield(normalise, 'inside')
  % convert back to a logical volume
  normalise.inside  = abs(normalise.inside-1)<=10*eps;
end

% flip and permute the dimensions to align the volume with the headcoordinate axes
normalise = align_ijk2xyz(normalise);

if strcmp(cfg.write,'yes')
  % create an spm-compatible file for each of the normalised volumes
  for parlop=1:length(cfg.parameter)  % include the anatomy
    tmp  = cfg.parameter{parlop};
    data = reshape(getsubfield(normalise, tmp), normalise.dim);
    tmp(find(tmp=='.')) = '_';
    volumewrite_spm([cfg.name,'_' tmp '.img'], data, normalise.transform);
  end
end

if strcmp(cfg.keepintermediate,'no')
  % remove the intermediate files
  for flop=1:length(files)
    [p, f, x] = fileparts(files{flop});
    delete(fullfile(p, [f, '.*']));
    [p, f, x] = fileparts(wfiles{flop});
    delete(fullfile(p, [f, '.*']));
  end
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
cfg.version.id = '$Id: volumenormalise.m,v 1.15 2007/04/03 15:37:07 roboos Exp $';
% remember the configuration details of the input data
try, cfg.previous = interp.cfg; end
% remember the exact configuration details in the output
normalise.cfg = cfg;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTION that asks a few questions to determine the coordinate system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [coordinates] = determine_coordinates(str);
fprintf('Please determine the coordinate system of the %s anatomical volume.\n', str);
% determine in which direction the nose is pointing
nosedir = [];
while isempty(nosedir)
  response = input('Is the nose pointing in the positive X- or Y-direction? [x/y] ','s');
  if strcmp(response, 'x'),
    nosedir = 'positivex';
  elseif strcmp(response, 'y'),
    nosedir = 'positivey';
  end
end
% determine where the origin is
origin  = [];
while isempty(origin)
  response = input('Is the origin on the Anterior commissure or between the Ears? [a/e] ','s');
  if strcmp(response, 'e'),
    origin = 'interauricular';
  elseif strcmp(response, 'a'),
    origin = 'ant_comm';
  end
end
% determine the coordinate system of the MRI volume
if strcmp(origin, 'interauricular') && strcmp(nosedir, 'positivex')
  coordinates = 'ctf';
elseif strcmp(origin, 'ant_comm') && strcmp(nosedir, 'positivey')
  coordinates = 'spm';
end
