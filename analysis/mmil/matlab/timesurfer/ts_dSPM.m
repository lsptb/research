function ts_dSPM(avg_data,subjname,varargin)
%function ts_dSPM(avg_data,subjname,[options]);
%
% Usage:
%  ts_dSPM(avg_data,subjname, 'key1', value1,...);
%
% Required Input:
%  avg_data - average data structure
%    (see ts_process_fif_data and ts_avg_fif_data)
%  subjname - freesurfer subject name (must be found in $SUBJECTS_DIR)
%
% Optional parameters:
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  'prefix' - prefix of all output files
%    {default: 'dSPM'}
%  'rootoutdir' - root output directory (several subdirectories will be created)
%    {default: pwd} (current working directory)
%  'conditions' - vector of condition numbers to analyze
%    {default: []} (if empty, use all conditions found in avg_data.averages)
%  'badchans'    - vector of bad channel indices
%     { default: [] }
%  'badchanfile' - name of text file containing bad channel labels
%    {default: []}
%  'usegrad_flag' - [1|0] Toggle use of gradiometer data, if available
%     {default: 1}
%  'usemag_flag' - [1|0] Toggle use of magnitometer data, if available
%     {default: 1}
%  'useEEG_flag' - [1|0] Toggle use of EEG data, if available
%     {default: 1}
%  'baseline_flag'  - [1|0] Toggle baseline correction of input data
%     { default: 1 }
%  'baseline_start   - start time of baseline period (msec)
%     if ncov_type = 1, use this time window for baseline correction
%       and noise covariance and scaling factor calculation from average data
%     relative to trigger onset; negative times occur before trigger
%     { default: -Inf } (start at beginning of prestimulus period)
%  'baseline_end'   - end time of baseline period (msec)
%     { default: 0 } (end at trigger onset)
%  'ssp_projmat' - Signal Space Projection matrix (nchan x nchan)
%     will be applied to both the input data and the forward gain matrix
%     {default: []}
%
% Optional parameters for forward model:
%  'overwrite_forward_flag' - [1|0] Toggle overwriting of forward solution
%     if 0, will used pre-calculated forward gain matrix if it exists
%     {default: 0}
%  'lh_dip_file' - name of left hemi freesurfer dipole file (from tksurfer)
%    {default: 'bem/lh_white.dip'} (relative to subject dir)
%  'rh_dip_file' - name of right hemi freesurfer dipole file (from tksurfer)
%    {default: 'bem/rh_white.dip'} (relative to subject dir)
%  'lh_dip_info' - matrix of left hemi dipole information (6 x ndips)
%    6 columns for x,y,z coordinates and nx,ny,nz normal vector
%    if specified, will ignore lh_dip_file
%    {default: []}
%  'rh_dip_info' - matrix of right hemi dipole information (6 x ndips)
%    6 columns for x,y,z coordinates and nx,ny,nz normal vector
%    if specified, will ignore rh_dip_file
%    {default: []}
%  'lh_dec_file' - name of left hemi freesurfer dec file (from tksurfer)
%     {default: 'bem/lh_white_7.dec'} (relative to subject dir)
%     dec file contains 0's and 1's indicating subset of dipoles to use
%  'rh_dec_file' - name of right hemi freesurfer dec file (from tksurfer)
%     {default: 'bem/rh_white_7.dec'} (relative to subject dir)
%  'lh_dec_dips' - vector of 0's and 1's indicating subset of left hemi dipoles
%     length must match ndips (nvertices)
%     if specified, will ignore lh_dec_file
%     {default: []}
%  'rh_dec_dips' - vector of 0's and 1's indicating subset of right hemi dipoles
%     if specified, will ignore rh_dec_file
%     {default: []}
%  'bem_flag' - [1|0] Toggle use of boundary element method (BEM)
%     if 0, use spherical shell model instead
%     {default: 1}
%  'nlayers' - number of layers for BEM or spherical shell model
%     For MEG only, can be 1 or 3; with EEG (alone or with MEG), must be 3
%     {default: 3}
%  'bem_surf_files' - cell array of file names for FreeSurfer tri files
%    (text files specifying triangles and faces) with MRI-derived skull
%     and scalp surfaces
%    This should contain three filenames (each relative to subject dir):
%      (1) inner skull
%      (2) outer skull
%      (3) outer scalp
%    Optionally, for single-shell MEG-only gain matrix, a single filename
%     for the inner skull surface can be specified here.
%    Files must exist if bem_flag=1
%    {default: {'bem/inner_skull4.tri'...
%               'bem/outer_skull4.tri'...
%               'bem/outer_scalp4.tri'}     }
%  'cen_sph' - center of sphere (mm) - vector of 3 numbers (x,y,z)
%     ignored if bem_flag=1
%     {default: []}
%  'radii' - vector of three values specifying radii (mm) of spherical shell layers:
%     (1) inner skull (2) outer skull (3) outer scalp
%     Required if useEEG_flag=1 and bem_flag=0 (spherical shell model)
%     {default: []}
%  'conductivities' - vector of three values specifying conductivity of:
%     (1) brain (2) skull  (3) scalp
%     {default: [0.3 0.012 0.3]}
%  'conduct_scalefact' - conductivity scaling factor
%     an optimized value is necessary for integration of MEG and EEG data
%     {default: 1}
%  'refEEG_coords' - reference EEG electrode coordinates (e.g. [0.004,0.020,0.001])
%     if specified, forward solution for these coordinates is subtracted 
%       from forward solution for all other EEG electrodes
%     for data originally in fif format, these coordinates are already
%       stored in the "loc" matrix for each EEG electrode
%     {default: []}
%
% Optional parameters for coordinate transformation:
%  'trans' - 4x4 matrix specifying mri2head transformation
%     {default: []}
%  'transfile' - text file containing 4x4 mri2head transformation
%     {default: []}
%  'alignment_fif' fif file containing 4x4 mri2head transformation
%     note on reading from fif's:
%      This requires Uutela's fifaccess toolbox which only runs on 32 bit Matlab
%      and core dump will result if mri2head transform is missing from
%      specified file
%     {default: []}
%
% Note on transform files:
%   trans, trans_file, and alignment_fiff are different, mutually exclusive
%     ways to specify the mri2head transformation matrix.
%   This transformation matrix is a 4x4 afine transformation that
%     maps locations in the MRI (from dip files) to locations in the
%     MEG/EEG head space defined by the digitized cardinal points
%     (and HPI coils).
%   This transformation matrix can also be found in the avg_data structure.
%     The avg_data structure contains a field called coor_trans.
%     If the subfield "mri2head" exists and is not empty, this will be used.
%   If more than one transformation marix is supplied, the order of dominance is:
%     trans, trans_file, alignment_fiff, avg_data.coor_trans.mri2head
%     i.e., if trans is supplied directly, all other methods of
%     supplying the trans will be ignored.
%
% Optional parameters for noise covariance:
%  'ncov_type' - specify what type of noise covariance matrix to use
%     if 0, use identity matrix
%       (assume uniform white noise, independently scaled for each sensor type)
%     if 1, calculate and use noise covariance matrix from average prestim
%     if 2, use noise covariance matrix calculated from single trials
%       (stored in avg_data.noise.covar)
%    {default: 2}
%  'ncov_conditions' - vector of condition numbers to use for noise covariance
%    matrix calculation (if empty, will use all in 'conditions')
%    {default: []}
%
% Optional parameters for fMRI bias:
%  'lh_sourcecov_wfile' - full pathname of left hemi freesurfer w file
%     file can contain fMRI activations on surface for a priori weighting
%     if not supplied, will used identity matrix as source covariance
%     {default: []}
%  'rh_sourcecov_wfile' - full pathname of right hemi freesurfer w file
%     {default: []}
%  'sourcecov_thresh' - threshold applied to values in sourcecov wfiles
%     suprathreshold vertices (dipoles) will get weighting of sourcecov_maxvar
%     subthreshold vertices (dipoles) will get weighting of sourcecov_minvar
%     {default: 0}
%  'sourcecov_thresh_abs_flag' - [1|0] Toggle absolute value threshold
%     {default: 1}
%  'sourcecov_maxvar' - value assigned to diagonal element for suprathreshold
%      dipoles in source covariance matrix
%     {default: 0.9} (maximum allowed value is 0.999)
%  'sourcecov_minvar' - value assigned to diagonal element for suprathreshold
%      dipoles in source covariance matrix
%     {default: 0.09} (minimum allowed value is 0.001)
%
% Optional parameters for smoothness constraint:
%  'smooth_constr_flag' - [1|0] Toggle smoothness constraint
%     If 1, forces orient_constr_flag = 1
%     If 1, lh/rh_dec_file and lh/rh_dec_dips are ignored
%       Instead, subsampled surface will be generated using matlab function
%         reducepatch so that surface smoothing can be done
%     {default: 0}
%  'smooth_constr_nsmooth' - number of smoothing steps for smoothness constraint
%     Smoothing is done on subsampled surfaces
%     {default: 10}
%  'smooth_constr_subfact' - subsampling factor for creation of subsampled
%     surface file from subject's white freesurfer surfaces
%     Value of 0.02 gives spacing of ~7 mm with ~2500 verts per hemi
%     Min = 0.01, Max = 0.1
%     {default: 0.02}
%
% Optional parameters for inverse:
%  'overwrite_inverse_flag' - [1|0] Toggle overwriting of inverse operator
%     if 0, will used pre-calculated inverse operator otoutdir'.
%     if it exists
%     {default: 0}
%  'orient_constr_flag' - [1|0] Toggle orientation constraint
%     if 1, dipole orientations will be fixed to be perpendicular to cortex
%     if 0, dipole orientations will be free to rotate
%     {default: 0}
%  'orient_tang' - source covariance weight of tangential components
%       when using orientation constraint
%     Value of 0 gives a fixed orientation constraint
%     Value of 1 allows free orientations (equal weighting with normal vector)
%     Intermediate value gives a "loose" orientation constraint
%     Only applicable when orient_constr_flag=1
%     {default: 0}
%  'SNR' - estimated signal-to-noise-ratio (for regularization parameter)
%     {default: 10}
%  'amd_inverse_flag'  - [1|0] Toggle calculation of inverse operator
%     using Anders Dale's original method
%     Otherwise use Matti Hamalainen's pre-whitened inverse operator
%     {default: 1}
%  'calc_scalefacts_flag' - [1|0] Toggle calculation of scaling factors
%     for different channel types (standard deviation of avg baseline)
%     if 1, calculate these scaling factors and apply to data and forward matrix
%     if 0, use default or user specified scaling factors
%     NOTE: ignored if amd_inverse_flag=0
%     {default: 0}
%  'grad_scalefact' - scaling factor applied to gradiometer data
%     purpose of scaling factors is to get data from different channel types
%     into roughly the same scale
%     NOTE: ignored if amd_inverse_flag=0
%     {default: 10^13}
%  'mag_scalefact' - scaling factor applied to magnitometer data
%     NOTE: ignored if amd_inverse_flag=0
%     {default: 10^15}
%  'EEG_scalefact' - scaling factor applied to EEG data
%     NOTE: ignored if amd_inverse_flag=0
%     {default: 10^6}
%  'noisenorm_flag' - [1|0] Noise-normalization (Dale et al., Neuron 2000)
%     If 1, source waveforms are noise sensitivity normalized sqrt(F) or z-stats
%     If 0, source waveforms are in nA*m (but are depth-biased)
%     {default: 1}
%  'depthweight_flag' - [1|0] Depth-weighting (Lin et al., NeuroImage 2006)
%     If 1, source covariance is depth-weighted
%     If 0, source covariance is not depth-weighted
%     {default: 0}
%  'depthweight_p' - depth weighting parameter
%     {default: 0.5}
%  'signed_sources_flag' - [1|0] Toggle output of source waveforms
%       as signed amplitudes along normal vectors
%     If 1, and noisenorm_flag=1, source waveforms are z-stats
%     If 0, and noisenorm_flag=1, source waveforms are sqrt(F-stats)
%       The F-stats are unsigned because they are the sum of the
%        squared components of dipole vector
%     Only available if orient_constr_flag = 1
%     {default: 1}
%
% Optional parameters for output:
%  'overwrite_output_flag' - [1|0] Toggle overwrite of stc, mgh, fif files
%     if they already exist
%     {default: 1}
%  'write_stc_flag' - [1|0] Toggle output of source waveforms as stc files
%     forced to 1 if write_mgh_flag = 1
%     {default: 1}
%  'stc_scalefact' - scale factor multiplier for output stc file waveforms
%     {default: 1}
%  'write_mgh_flag' - [1|0] Toggle additional output of source waveforms
%     as mgh files (displayable with FreeSurfer's tksurfer)
%     forced to 1 if resamp2ico_flag = 1
%     {default: 0}
%  'sparsesmooth' - # of sparse smoothing steps applied before saving mgh file
%     see ts_stc2mgh for details on smoothing
%     {default: 10}
%  'postsmooth - # of normal smoothing steps applied when saving as mgh file
%     {default: 10}
%  'mbmask_flag' - [0|1] whether to mask midbrain & corpus callosum
%     {default: 0}
%  'resamp2ico_flag' - [1|0] Toggle resampling of mgh files to icosahedral sphere
%     necessary preparation for group analysis
%     (displayable on fsaverage subject with FreeSurfer's tksurfer)
%     {default: 0}
%  'icolevel' - icosahedron order number:
%               Order  Number of Vertices
%                 0              12
%                 1              42
%                 2             162
%                 3             642
%                 4            2562
%                 5           10242
%                 6           40962
%                 7          163842
%    {default: 7}
%  'icosmooth' - smoothing steps on surface after sampling to ico
%    {default: 3}
%  'write_fif_flag' - [1|0] Toggle output of fit and residual error
%     sensor waveforms as fif files (displayable with Neuromag's xplotter)
%     {default: 0}
%  'template_fif' - full path name of online average fif file
%     must specify if write_fif_flag=1
%     {default: []}
%
% created:   08/02/06    by Don Hagler
% last mod:  08/05/09    by Don Hagler
%

%% todo: flag to specify whether avg_data.noise.covar is from average or from raw
%% todo: option to pass in noise covar as parameter
%% todo: option to read data from fif file (need new function ts_fif2avg)
%% todo: extra dips (e.g. subcortical)

MIN_SOURCECOV = 0.001;
MAX_SOURCECOV = 0.999;

% scaling factors to get forward matrix into same scale as fif data
GRAD_UNITSFACT = 10^-13; % fT/cm to T/m
MAG_UNITSFACT = 10^-15; % fT to T
EEG_UNITSFACT = 10^-6; % uV to V

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,2)) return; end;
parms = mmil_args2parms(varargin, { ...
  'subjdir',getenv('SUBJECTS_DIR'),[],...
  'rootoutdir',pwd,[],...
  'prefix','dSPM',[],...
  'conditions',[],[],...
  'lh_dip_file','bem/lh_white.dip',[],...
  'rh_dip_file','bem/rh_white.dip',[],...
  'lh_dip_info',[],[],...
  'rh_dip_info',[],[],...
  'lh_dec_file','bem/lh_white_7.dec',[],...
  'rh_dec_file','bem/rh_white_7.dec',[],...
  'lh_dec_dips',[],[],...
  'rh_dec_dips',[],[],...
  'ncov_type',2,[0 2],...
  'identity_ncov_flag',false,[false true],... % for backward compatibility
  'calc_avg_ncov_flag',false,[false true],... % for backward compatibility
  'ncov_conditions',[],[],...
  'calc_scalefacts_flag',false,[false true],...
  'noise_start',[],[-Inf Inf],... % for backward compatibility
  'noise_end',[],[-Inf Inf],... % for backward compatibility
  'baseline_start',-Inf,[-Inf Inf],...
  'baseline_end',0,[-Inf Inf],...
  'baseline_flag',true,[false true],...
  'ssp_projmat',[],[],...
  'SNR',10,[eps Inf],...
  'noisenorm_flag',true,[false true],...
  'depthweight_flag',false,[false true],...
  'depthweight_p',0.5,[0,1],...
  'bem_flag',true,[false true],...
  'radii',[],[],...
  'conductivities',[0.3 0.012 0.3],[],...
  'conduct_scalefact',1,[0 Inf],...
  'nlayers',3,[1 3],...
  'badchans',[],[],...
  'badchanfile',[],[],...
  'usegrad_flag',true,[false true],...
  'usemag_flag',true,[false true],...
  'useEEG_flag',true,[false true],...
  'grad_scalefact',10^13,[-Inf Inf],...
  'mag_scalefact',10^15,[-Inf Inf],...
  'EEG_scalefact',10^6,[-Inf Inf],...
  'overwrite_output_flag',true,[false true],...
  'write_stc_flag',true,[false true],...
  'stc_scalefact',1,[eps Inf],...
  'write_mgh_flag',false,[false true],...
  'sparsesmooth',10,[0 1000],...
  'postsmooth',10,[0 1000],...
  'mbmask_flag',false,sort([false,true]),...
  'resamp2ico_flag',false,[false true],...
  'icolevel',7,[1 7],...
  'icosmooth',3,[0 1000],...
  'write_fif_flag',false,[false true],...
  'template_fif',[],[],...
  'bem_surf_files',...
    {'bem/inner_skull4.tri','bem/outer_skull4.tri','bem/outer_scalp4.tri'},[],...
  'cen_sph',[],[],...
  'trans',[],[],...
  'transfile',[],[],...
  'alignment_fif',[],[],...
  'lh_sourcecov_wfile',[],[],...
  'rh_sourcecov_wfile',[],[],...
  'sourcecov_thresh',0,[-Inf Inf],...
  'sourcecov_thresh_abs_flag',true,[false true],...
  'sourcecov_maxvar',0.9,[0 1],...
  'sourcecov_minvar',0.09,[0 1],...
  'overwrite_forward_flag',false,[false true],...
  'overwrite_inverse_flag',false,[false true],...
  'refEEG_coords',[],[],...
  'amd_inverse_flag',true,[false true],...
  'orient_constr_flag',false,[false true],...
  'orient_tang',0,[0 1],...
  'smooth_constr_flag',false,[false true],...
  'smooth_constr_nsmooth',10,[1,100],...
  'smooth_constr_subfact',0.02,[0.01,0.1],...
  'signed_sources_flag',true,[false true],...
},false);

avg_data  = ts_checkdata_header(avg_data,'conditions',parms.conditions);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

% subjname
parms.subjname = subjname;
if isempty(parms.subjdir)
  parms.subjdir = deblank(getenv('SUBJECTS_DIR'));
  if(isempty(parms.subjdir))
    error('Cannot find SUBJECTS_DIR environment variable');
  end
end;
parms.subjpath = [parms.subjdir '/' parms.subjname];
if ~exist(parms.subjpath,'dir')
  error('subject dir %s not found',parms.subjpath);
end

% conditions to use
if isempty(parms.conditions)
  parms.conditions = 1:length(avg_data.averages);
else
  conditions = 1:length(avg_data.averages);
  parms.conditions = intersect(parms.conditions,conditions);
end;
if isempty(parms.conditions)
  error('no valid conditions specified');
end;
if isempty(parms.ncov_conditions)
  parms.ncov_conditions = parms.conditions;
else
  conditions = 1:length(avg_data.averages);
  parms.ncov_conditions = intersect(parms.ncov_conditions,conditions);
end;
if isempty(parms.ncov_conditions)
  error('no valid ncov conditions specified');
end;
if parms.identity_ncov_flag % for backward compatibility
  parms.ncov_type = 0;
elseif parms.calc_avg_ncov_flag % for backward compatibility
  parms.ncov_type = 1;
end;

% mri2head transformation matrix
T_mri2head = [];
if ~isempty(parms.trans)
  fprintf('%s: using transformation matrix supplied...\n',mfilename);
  T_mri2head = parms.trans;
elseif ~isempty(parms.transfile)
  if ~exist(parms.transfile,'file')
    error('transfile %s not found',parms.transfile);
  end;
  fprintf('%s: reading transformation matrix from %s...\n',...
    mfilename,parms.transfile);
  T_mri2head = ts_read_transfile(parms.transfile);
elseif ~isempty(parms.alignment_fif)
  if ~exist(parms.alignment_fif,'file')
    error('alignment_fif %s not found',parms.alignment_fif);
  end
  fprintf('%s: reading transformation matrix from %s...\n',...
    mfilename,parms.alignment_fif);
  T_mri2head=loadtrans(parms.alignment_fif,'MRI','HEAD');
elseif isfield(avg_data.coor_trans,'mri2head')
  fprintf('%s: using transformation matrix in avg_data structure...\n',mfilename);
  if ~isempty(avg_data.coor_trans.mri2head)
    T_mri2head=avg_data.coor_trans.mri2head;
  end;
end;
if isempty(T_mri2head)
  fprintf('%s: no mri2head transformation was supplied\n',mfilename);
  return;  
end;
parms.trans = T_mri2head;
fprintf('%s: mri2head transformation matrix:\n',mfilename);
disp(parms.trans);

% dipole information
if isempty(parms.lh_dip_info)
  fprintf('%s: loading left hemisphere dipole information...\n',mfilename);
  parms.lh_dip_file = [parms.subjpath '/' parms.lh_dip_file];
  if ~exist(parms.lh_dip_file,'file')
    error('lh_dip_file %s not found',parms.lh_dip_file);
  end
  parms.lh_dip_info = ts_read_dip_file(parms.lh_dip_file);
end;
if isempty(parms.rh_dip_info)
  fprintf('%s: loading right hemisphere dipole information...\n',mfilename);
  parms.rh_dip_file = [parms.subjpath '/' parms.rh_dip_file];
  if ~exist(parms.rh_dip_file,'file')
    error('rh_dip_file %s not found',parms.rh_dip_file);
  end
  parms.rh_dip_info = ts_read_dip_file(parms.rh_dip_file);
end;
parms.num_dips_lh=length(parms.lh_dip_info);
parms.num_dips_rh=length(parms.rh_dip_info);
if(parms.num_dips_lh==0 & parms.num_dips_rh==0)
  error('no dipole info found');
end

% decimated (subset) dipoles
if parms.smooth_constr_flag % subsampled surface instead
  if ~parms.orient_constr_flag
    fprintf('%s: smooth_constr_flag=1, so forcing orient_constr_flag=1',...
      mfilename);
    parms.orient_constr_flag = 1;
  end;
  clear lh_dec_dips;
  matfile=sprintf('%s/matfiles/lh_white_sub%0.2f.mat',...
    parms.rootoutdir,parms.smooth_constr_subfact);
    % assume that every dSPM run will have same subject and surfaces
    % so don't have a prefix-specific matfile
  if ~exist(matfile,'file')
    parms.lh_surf_file = [parms.subjpath '/surf/lh.white'];
    if ~exist(parms.lh_surf_file,'file')
      error('lh_surf_file %s not found',parms.lh_surf_file);
    end;
    % subsample surface, etc.
    fprintf('%s: subsampling surface %s with subfact %0.2f...\n',...
      mfilename,parms.lh_surf_file,parms.smooth_constr_subfact);
    [lh_subsurf,lh_surf,lh_dec_dips] = ...
      ts_subsamp_surf(parms.lh_surf_file,parms.smooth_constr_subfact);
    save(matfile,'lh_surf','lh_subsurf','lh_dec_dips');
  else
    fprintf('%s: loading subsampled surface from %s...\n',...
      mfilename,matfile);
    load(matfile);
  end;
  fprintf('\toriginal:    %d verts, avg dist = %0.1f mm\n',...
    lh_surf.nverts,lh_surf.avgdist);
  fprintf('\tsubsampled:  %d verts, avg dist = %0.1f mm\n',...
    lh_subsurf.nverts,lh_subsurf.avgdist);
  parms.lh_dec_dips = lh_dec_dips;
elseif isempty(parms.lh_dec_dips)
  parms.lh_dec_file = [parms.subjpath '/' parms.lh_dec_file];
  if ~exist(parms.lh_dec_file,'file')
    error('lh_dec_file %s not found',parms.lh_dec_file);
  end;
  fprintf('%s: loading lh dec file %s...\n',mfilename,parms.lh_dec_file);
  parms.lh_dec_dips=ts_read_dec_file(parms.lh_dec_file);
end

if parms.smooth_constr_flag % subsampled surface instead
  clear rh_dec_dips;
  matfile=sprintf('%s/matfiles/rh_white_sub%0.2f.mat',...
    parms.rootoutdir,parms.smooth_constr_subfact);
    % assume that every dSPM run will have same subject and surfaces
    % so don't have a prefix-specific matfile
  if ~exist(matfile,'file')
    parms.rh_surf_file = [parms.subjpath '/surf/rh.white'];
    if ~exist(parms.rh_surf_file,'file')
      error('rh_surf_file %s not found',parms.rh_surf_file);
    end;
    % subsample surface, etc.
    fprintf('%s: subsampling surface %s with subfact %0.2f...\n',...
      mfilename,parms.rh_surf_file,parms.smooth_constr_subfact);
    [rh_subsurf,rh_surf,rh_dec_dips] = ...
      ts_subsamp_surf(parms.rh_surf_file,parms.smooth_constr_subfact);
    save(matfile,'rh_surf','rh_subsurf','rh_dec_dips');
  else
    fprintf('%s: loading subsampled surface from %s...\n',...
      mfilename,matfile);
    load(matfile);
  end;
  fprintf('\toriginal:    %d verts, avg dist = %0.1f mm\n',...
    rh_surf.nverts,rh_surf.avgdist);
  fprintf('\tsubsampled:  %d verts, avg dist = %0.1f mm\n',...
    rh_subsurf.nverts,rh_subsurf.avgdist);
  parms.rh_dec_dips = rh_dec_dips;
elseif isempty(parms.rh_dec_dips)
  parms.rh_dec_file = [parms.subjpath '/' parms.rh_dec_file];
  if ~exist(parms.rh_dec_file,'file')
    error('rh_dec_file %s not found',parms.rh_dec_file);
  end;
  fprintf('%s: loading rh dec file %s...\n',mfilename,parms.rh_dec_file);
  parms.rh_dec_dips=ts_read_dec_file(parms.rh_dec_file);
end
parms.num_dec_dips_lh=length(find(parms.lh_dec_dips));
parms.num_dec_dips_rh=length(find(parms.rh_dec_dips));
if(parms.num_dec_dips_lh==0 & parms.num_dec_dips_rh==0)
  error('no decimated dipoles specified');
end
if length(parms.lh_dec_dips) ~= parms.num_dips_lh
  error('length of lh_dec_dips must match number of lh dips');
end
if length(parms.rh_dec_dips) ~= parms.num_dips_rh
  error('length of rh_dec_dips must match number of rh dips');
end
fprintf('%s: number of selected dipoles: lh=%d rh=%d\n',...
  mfilename,parms.num_dec_dips_lh,parms.num_dec_dips_rh);

% usable channels
parms.grad_chans = find(strncmp('grad',lower({avg_data.sensor_info.typestring}),...
             length('grad')));
parms.mag_chans  = find(strcmp('mag',lower({avg_data.sensor_info.typestring})));
parms.EEG_chans = find(strcmp('eeg',lower({avg_data.sensor_info.typestring})));
if isempty(parms.grad_chans) & parms.usegrad_flag
  fprintf('%s: WARNING: no gradiometers found, setting usegrad_flag = 0\n',...
    mfilename);
  parms.usegrad_flag=0;
end;
if isempty(parms.mag_chans) & parms.usemag_flag
  fprintf('%s: WARNING: no magnetometers found, setting usemag_flag = 0\n',...
    mfilename);
  parms.usemag_flag=0;
end;
if isempty(parms.EEG_chans) & parms.useEEG_flag
  fprintf('%s: WARNING: no EEG channels found, setting useEEG_flag = 0\n',...
    mfilename);
  parms.useEEG_flag=0;
end;
parms.useMEG_flag = (parms.usemag_flag | parms.usegrad_flag);
if ~parms.useMEG_flag & ~parms.useEEG_flag
  error('one of usegrad_flag, usemag_flag, or useEEG_flag must be on');
end;
parms.nchans = length(avg_data.sensor_info);

if ~isempty(parms.ssp_projmat)
  pm_sz = size(parms.ssp_projmat);
  if length(pm_sz)~=2 | any(pm_sz~=[parms.nchans,parms.nchans])
    fprintf('%s: ssp_projmat size:\n',mfilename);
    disp(pm_sz);
    error('ssp_projmat has wrong number of dimensions');
  end;
  parms.ssp_projmat = double(full(parms.ssp_projmat));
end;

if parms.useEEG_flag & ~isempty(parms.refEEG_coords)
  if length(parms.refEEG_coords)~=3
    error('regEEG_coords must contain 3 values (x,y,z)');
  end;
  % replace/add reference coords to loc for each electrode
  for i=1:avg_data.num_sensors
    avg_data.sensor_info(i).loc(1:3,1) = parms.refEEG_coords;
  end;
end;

% bad chans
badchan_i = find(cell2mat({avg_data.sensor_info.badchan})==1);
parms.badchans = union(parms.badchans,badchan_i);
labels = {avg_data.sensor_info.label}; 
if ~isempty(parms.badchanfile)
  badchan_i = ts_read_txt_badchans(parms.badchanfile,labels);
else
  badchan_i = [];
end;
parms.badchans = union(parms.badchans,badchan_i);

% good chans
parms.goodchans = [];
if parms.usegrad_flag
  parms.goodchans = union(parms.goodchans,parms.grad_chans);
end;
if parms.usemag_flag
  parms.goodchans = union(parms.goodchans,parms.mag_chans);
end;
if parms.useEEG_flag
  parms.goodchans = union(parms.goodchans,parms.EEG_chans);
end;
parms.goodchans = setdiff(parms.goodchans,parms.badchans);
if isempty(parms.goodchans)
  error('no good channels left');
end;

% EEG must have 3 layers
if parms.nlayers<3 & parms.useEEG_flag
  error('must have 3 layers (BEM surfacse or spherical shells) for EEG');
end;
% MEG can have 1 or 3 layers
if parms.nlayers~=1 & parms.nlayers~=3
  error('nlayers must be 1 or 3');
end;

% bem vs. sphere
if (parms.bem_flag)
  if iscell(parms.bem_surf_files)
    num_bem_surfs=length(parms.bem_surf_files);
  else
    parms.bem_surf_files={parms.bem_surf_files};
    num_bem_surfs=1;
  end;
  if num_bem_surfs > parms.nlayers % i.e. nlayers = 1
    parms.bem_surf_files = {parms.bem_surf_files{1:parms.nlayers}};
    num_bem_surfs=length(parms.bem_surf_files);
  elseif num_bem_surfs < parms.nlayers
    error('number of bem surf files must match nlayers');
  end;
  % check that files exist and change from relative to absolute path
  for s=1:num_bem_surfs
    fname = parms.bem_surf_files{s};
    fname = [parms.subjpath '/' fname];
    if ~exist(fname,'file')
      error('bem surf %d file %s not found',s,fname);
    end
    parms.bem_surf_files{s}=fname;
  end;
  fprintf('%s: will use BEM forward with surface files:\n',mfilename);
  for s=1:num_bem_surfs
    fprintf('   %d: %s\n',s,parms.bem_surf_files{s});
  end;
else
  if parms.useEEG_flag
    if length(parms.radii) > parms.nlayers
      parms.radii = parms.radii(1:parms.nlayers);
    elseif length(parms.radii) < parms.nlayers
      error('number of radii must match number of layers');
    end;
  end;
  if isempty(parms.cen_sph)
    error('cen_sph must be specified if bem_flag=0');
  end
  if(length(parms.cen_sph)~=3)
    error('cen_sph (center of sphere) incorrectly specified');
  end
  fprintf('%s: spherical shell forward with center of sphere: [%0.1f,%0.1f,%0.1f]\n',...
    mfilename,parms.cen_sph);
  if parms.nlayers==1
    fprintf('%s: with radius: %0.1f',...
      mfilename,parms.radii);
  else
    fprintf('%s: with radii: [%0.1f,%0.1f,%0.1f]\n',...
      mfilename,parms.radii);
  end;
end
if length(parms.conductivities) > parms.nlayers
  parms.conductivities = parms.conductivities(1:parms.nlayers);
elseif length(parms.conductivities) < parms.nlayers
  error('number of conductivities must match number of layers');
end;

% source covariance for smoothness constraint
if parms.smooth_constr_flag
  matfile=sprintf('%s/matfiles/Rsmooth%d.mat',...
    parms.rootoutdir,parms.smooth_constr_nsmooth);
    % assume that every dSPM run will have same subject and surfaces
    % so don't have a prefix-specific matfile
  if ~exist(matfile,'file')
    fprintf('%s: creating source covariance matrix for smoothness constraint...\n',...
      mfilename);
    Rsmooth_lh = ts_smoothconstr_covar(lh_subsurf,parms.smooth_constr_nsmooth);
    Rsmooth_rh = ts_smoothconstr_covar(rh_subsurf,parms.smooth_constr_nsmooth);
    save(matfile,'Rsmooth_lh','Rsmooth_rh');
  else
    load(matfile);
  end;
end;

ndips_lh=3*parms.num_dec_dips_lh;
ndips_rh=3*parms.num_dec_dips_rh;
if ~isempty(parms.lh_sourcecov_wfile) || ~isempty(parms.rh_sourcecov_wfile)
  % source covariance from w files
  parms.sourcecov_thresh = max(parms.sourcecov_thresh,0);
  parms.sourcecov_minvar = max(parms.sourcecov_minvar,MIN_SOURCECOV);
  parms.sourcecov_maxvar = max(parms.sourcecov_maxvar,MIN_SOURCECOV);
  parms.sourcecov_minvar = min(parms.sourcecov_minvar,MAX_SOURCECOV);
  parms.sourcecov_maxvar = min(parms.sourcecov_maxvar,MAX_SOURCECOV);
  if ~isempty(parms.lh_sourcecov_wfile)
    fprintf('%s: loading lh sourcecov wfile %s...\n',...
      mfilename,parms.lh_sourcecov_wfile);
    R_lh=ts_wfile2sourcecov(parms.lh_sourcecov_wfile,...
      'decdips',parms.lh_dec_dips,...
      'thresh',parms.sourcecov_thresh,...
      'threshabs',parms.sourcecov_thresh_abs_flag,...
      'maxvar',parms.sourcecov_maxvar,...
      'minvar',parms.sourcecov_minvar);
    if isempty(R_lh)
      error('cannot read lh sourcecov wfile %s',parms.lh_sourcecov_wfile);
    end;
  else
    R_lh=speye(ndips_lh,ndips_lh)*parms.sourcecov_minvar;
  end
  if ~isempty(parms.rh_sourcecov_wfile)
    fprintf('%s: loading rh sourcecov wfile %s...\n',...
      mfilename,parms.rh_sourcecov_wfile);
    R_rh=ts_wfile2sourcecov(parms.rh_sourcecov_wfile,...
      'decdips',parms.rh_dec_dips,...
      'thresh',parms.sourcecov_thresh,...
      'threshabs',parms.sourcecov_thresh_abs_flag,...
      'maxvar',parms.sourcecov_maxvar,...
      'minvar',parms.sourcecov_minvar);
    if isempty(R_rh)
      error('cannot read rh sourcecov wfile %s',parms.rh_sourcecov_wfile);
      return;
    end;
  else
    R_rh=speye(ndips_rh,ndips_rh)*parms.sourcecov_minvar;
  end
  if parms.smooth_constr_flag
    fprintf('%s: combining smoothness constraint and fMRI bias...\n',...
      mfilename);
    %% todo: test this
    tmp_R_lh = full(diag(R_lh));
    for i=1:length(tmp_R_lh)
      Rsmooth_lh(i,:) = Rsmooth_lh(i,:)*tmp_R_lh(i);
    end;
    tmp_R_rh = full(diag(R_rh));
    for i=1:length(tmp_R_rh)
      Rsmooth_rh(i,:) = Rsmooth_rh(i,:)*tmp_R_rh(i);
    end;
  end;
else
  % source covariance as identity matrix
  if parms.smooth_constr_flag
    R_lh = Rsmooth_lh;
  else
    R_lh = speye(ndips_lh,ndips_lh);
  end;
  if parms.smooth_constr_flag
    R_rh = Rsmooth_rh;
  else
    R_rh = speye(ndips_rh,ndips_rh);
  end;
end;

fprintf('%s: combining left and right source covariance matrices...\n',...
  mfilename);
% combine R_lh and R_rh
if parms.smooth_constr_flag
  tic
  ndips = ndips_lh + ndips_rh;
  R=sparse(ndips,ndips);
  R(1:ndips_lh,1:ndips_lh) = R_lh;
  R(ndips_lh+1:end,ndips_lh+1:end) = R_rh;
  toc  % takes about 17 seconds
else  % diagonalo only
  tic
  ndips = ndips_lh + ndips_rh;
  R=sparse(ndips,ndips);
  for i=1:ndips_lh
    R(i,i)=R_lh(i,i);
  end
  j=1;
  for i=ndips_lh+1:ndips
    R(i,i)=R_rh(j,j);
    j=j+1;
  end
  toc % takes < 1 second
end;

% orient_tang for 2nd and 3rd components
if parms.orient_constr_flag
  fprintf('%s: setting diagonal tangential components of source covariance matrix...\n',...
    mfilename);
  tic
  ndips = ndips_lh + ndips_rh;
  n = ndips/3;
  j = 1;
  for i=1:n
    for k=j+1:j+2
      R(k,k) = R(j,j)*parms.orient_tang;
    end;
    j = j+3;
  end;
  toc
end;

parms.srccovar = R;

% apply projection matrix to data
if ~isempty(parms.ssp_projmat)
  fprintf('%s: applying SSP projection matrix to data...\n',...
    mfilename);
  nconds = length(parms.conditions);
  for c=1:nconds
    cond = parms.conditions(c);
    avg_data.averages(cond).data = parms.ssp_projmat*avg_data.averages(cond).data;
  end;
end;

% for backward compatibility
if ~isempty(parms.noise_start)
  parms.baseline_start = parms.noise_start;
end;
if ~isempty(parms.noise_end)
  parms.baseline_end = parms.noise_end;
end;

% baseline correction
if parms.baseline_flag
  fprintf('%s: subtracting baseline from average data...\n',...
    mfilename);
  nconds = length(parms.conditions);
  for c=1:nconds
    cond = parms.conditions(c);
    time = avg_data.averages(cond).time;
    [tmp,parms.baseline_start_samp] = min(abs(time - parms.baseline_start/1000));
    [tmp,parms.baseline_end_samp] = min(abs(time - parms.baseline_end/1000));
    data = avg_data.averages(cond).data;
    baseline = data(:,parms.baseline_start_samp:parms.baseline_end_samp);
    mean_base=mean(baseline')';
    data=data-mean_base*ones(1,size(data,2));  % correct baseline
    avg_data.averages(cond).data = data;    
  end;
end;

% calculate scaling factors (and noise covariance matrix)
if parms.ncov_type==1 || parms.calc_scalefacts_flag
  fprintf('%s: calculating scaling factors from average data...\n',...
    mfilename);
  C = zeros(parms.nchans);
  total_num_trials = 0;
  nconds = length(parms.ncov_conditions);
  for c=1:nconds
    cond = parms.ncov_conditions(c);
    total_num_trials = total_num_trials + avg_data.averages(cond).num_trials;
    time = avg_data.averages(cond).time;
    [tmp,parms.baseline_start_samp] = min(abs(time - parms.baseline_start/1000));
    [tmp,parms.baseline_end_samp] = min(abs(time - parms.baseline_end/1000));
    data = avg_data.averages(cond).data;
    baseline = data(:,parms.baseline_start_samp:parms.baseline_end_samp);
    mean_base=mean(baseline')';
    data=data-mean_base*ones(1,size(data,2));  % correct baseline
    baseline = data(:,parms.baseline_start_samp:parms.baseline_end_samp);
    num_trials = avg_data.averages(cond).num_trials;
    if isempty(num_trials) | num_trials<2
      num_trials = 2;
    end;
    C=C+diag(var(baseline'))*(num_trials-1);  % diag covariance matrix
  end;
  if isempty(total_num_trials) | total_num_trials<nconds+1
    total_num_trials = nconds+1; % just make it not blow up
    % it's your problem if you don't give any information about trial numbers
  end;
  C = C/(total_num_trials-nconds); % pooled variance across conditions
  scalefacts = sqrt(diag(C));
end;

% calculate scale factors
if parms.calc_scalefacts_flag==1 & parms.amd_inverse_flag==1
  grad_chans = setdiff(parms.grad_chans,parms.badchans);
  mag_chans = setdiff(parms.mag_chans,parms.badchans);
  EEG_chans = setdiff(parms.EEG_chans,parms.badchans);
  if isempty(grad_chans)
    parms.grad_scalefact = 0;
  else
    parms.grad_scalefact = 1/mean(scalefacts(grad_chans));
  end;
  if isempty(mag_chans)
    parms.mag_scalefact = 0;
  else
    parms.mag_scalefact = 1/mean(scalefacts(mag_chans));
  end;
  if isempty(EEG_chans)
    parms.EEG_scalefact = 0;
  else
    parms.EEG_scalefact = 1/mean(scalefacts(EEG_chans));
  end;
elseif ~parms.amd_inverse_flag
  parms.grad_scalefact = 1;
  parms.mag_scalefact = 1;
  parms.EEG_scalefact = 1;
end;
% otherwise use default or user supplied scaling factors
fprintf('%s: scaling factors: grad=%0.2g, mag=%0.2g, EEG=%0.2g\n',...
  mfilename,parms.grad_scalefact,parms.mag_scalefact,parms.EEG_scalefact);

% use calculated noise covariance matrix or get from avg_data.noise.covar
if parms.ncov_type==0
  fprintf('%s: will use identity matrix for noise covariance\n',...
    mfilename);
  C=[];
elseif parms.ncov_type==1
  fprintf('%s: will use noise covariance calculated from average baseline\n',...
    mfilename);
elseif isempty(avg_data.noise.covar)
  fprintf('%s: will use identity matrix for noise covariance\n',...
    mfilename);
  C = [];
else
  fprintf('%s: will use noise covariance matrix from avg_data.noise.covar\n',...
    mfilename);
  C = avg_data.noise.covar;
end;

% apply scale factors to noise covariance matrix
if ~isempty(C)
  scalefacts = zeros(parms.nchans,1);
  scalefacts(parms.grad_chans)=parms.grad_scalefact;
  scalefacts(parms.mag_chans)=parms.mag_scalefact;
  scalefacts(parms.EEG_chans)=parms.EEG_scalefact;
  scale_matrix = scalefacts*scalefacts';
  C = C.*scale_matrix;
  C = C(parms.goodchans,parms.goodchans);
end;
parms.noisecovar = C;

% if writing fifs, need template
if parms.write_fif_flag
  if ~exist(parms.template_fif,'file')
    error('template_fif %s not found',parms.template_fif);
  end;
end;

% create output subdirs
[success,msg,msgid] = mkdir(parms.rootoutdir,'matfiles');
if ~success
  error('unable to create output directory %s/matfiles: %s',...
    parms.rootoutdir,msg);
end;
if(parms.write_fif_flag)
  [success,msg,msgid] = mkdir(parms.rootoutdir,'fifs');
  if ~success
    error('unable to create output directory %s/fifs: %s',...
	    parms.rootoutdir,msg);
  end;
end;
if(parms.resamp2ico_flag)
  parms.write_mgh_flag = 1;
end;
if(parms.write_mgh_flag)
  [success,msg,msgid] = mkdir(parms.rootoutdir,'mghfiles');
  if ~success
    error('unable to create output directory %s/mghfiles',...
	    parms.rootoutdir,msg);
  end;
  parms.write_stc_flag = 1;
end;
if(parms.write_stc_flag)
  [success,msg,msgid] = mkdir(parms.rootoutdir,'stcfiles');
  if ~success
    error('unable to create output directory %s/stcfiles',...
	    parms.rootoutdir,msg);
  end;
end;

% save avg_data as matfile (for future reference)
fprintf('%s: saving avg_data to mat file...\n',mfilename);
matfile=sprintf('%s/matfiles/%s_avg_data.mat',parms.rootoutdir,parms.prefix);
save(matfile,'avg_data');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate forward solution

matfile=sprintf('%s/matfiles/%s_forward.mat',parms.rootoutdir,parms.prefix);
if ~exist(matfile,'file') || parms.overwrite_forward_flag
  fprintf('%s: calculating forward solution gain matrix...\n',mfilename);
  if parms.bem_flag
    bem_matfile = sprintf('%s/matfiles/%s_bem.mat',parms.rootoutdir,parms.prefix);
    if exist(bem_matfile,'file') && parms.overwrite_forward_flag
      delete(bem_matfile);
    end;
    G_xyz=ts_calc_bem_gainmat(avg_data,parms.lh_dip_info,parms.rh_dip_info,...
      parms.lh_dec_dips,parms.rh_dec_dips,parms.trans,...
      bem_matfile,parms.bem_surf_files,parms.conductivities,...
      parms.useMEG_flag,parms.useEEG_flag);
  else
    G_xyz=ts_calc_sph_gainmat(avg_data,parms.lh_dip_info,parms.rh_dip_info,...
      parms.lh_dec_dips,parms.rh_dec_dips,parms.trans,...
      parms.cen_sph,parms.radii,parms.conductivities,...
      parms.useMEG_flag,parms.useEEG_flag);
  end;
  if isempty(G_xyz)
    error('failed to calculate gain matrix');
  end;
  save(matfile, 'G_xyz');
else
  fprintf('%s: loading pre-calculated forward solution...\n',mfilename);
  load(matfile);
end;

% apply ssp projection matrix
if ~isempty(parms.ssp_projmat)
  fprintf('%s: applying SSP projection matrix to gain matrix...\n',mfilename);
  G_xyz = parms.ssp_projmat*G_xyz;
end;

% scale gain matrix
%fprintf('%s: scaling gain matrix...\n',mfilename);
% scale to get to units of data, then scale with scalefacts
% and scale EEG part with conduct_scalefact
scalefacts = zeros(parms.nchans,1);
scalefacts(parms.grad_chans)=GRAD_UNITSFACT*parms.grad_scalefact;
scalefacts(parms.mag_chans)=MAG_UNITSFACT*parms.mag_scalefact;
scalefacts(parms.EEG_chans)=EEG_UNITSFACT*parms.EEG_scalefact*parms.conduct_scalefact;
scale_matrix = scalefacts*ones(1,size(G_xyz,2));
G_xyz = G_xyz.*scale_matrix;

if parms.orient_constr_flag
  % convert xyz to normal and two orthogonal tangential components
  fprintf('%s: converting xyz components to normal and tangentials\n',mfilename);
  [G_norm,G_tang1,G_tang2]=ts_gain_xyz2norm(G_xyz,parms.lh_dip_info,...
    parms.rh_dip_info,parms.lh_dec_dips,parms.rh_dec_dips,parms.trans);
  tmp = zeros(size(G_xyz));
  tmp(:,1:3:end) = G_norm;
  tmp(:,2:3:end) = G_tang1;
  tmp(:,3:3:end) = G_tang2;
  G_xyz = tmp;
end;

% resize to remove unwanted and bad chans
fprintf('%s: resizing gain matrix to remove bad channels...\n',mfilename);
G_xyz = G_xyz(parms.goodchans,:);
[num_sensors,num_sources]=size(G_xyz);
fprintf('%s: size(G_xyz)=[%d,%d]\n',mfilename,num_sensors,num_sources);
num_sources = num_sources/3; % for 3 orientations

% depth weighting
if parms.depthweight_flag
  R = parms.srccovar;
  R_dw = ts_depthweight_covar(G_xyz,parms.depthweight_p);
  if parms.smooth_constr_flag
    fprintf('%s: combining smoothness constraint and fMRI bias...\n',...
      mfilename);
    %% todo: test this
    tmp_R_dw = full(diag(R_dw));
    for i=1:length(tmp_R_dw)
      R(i,:) = R(i,:)*tmp_R_dw(i);
    end;
  else
    R = R.*R_dw;
  end;
  parms.srccovar = R;
end;

% save parms as matfile (for future reference)
fprintf('%s: saving parameters to mat file...\n',mfilename);
matname=sprintf('%s/matfiles/%s_parms.mat',parms.rootoutdir,parms.prefix);
save(matname,'parms');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate inverse operator

matfile=sprintf('%s/matfiles/%s_inverse.mat',parms.rootoutdir,parms.prefix);
if ~exist(matfile,'file') || parms.overwrite_inverse_flag
  fprintf('%s: calculating inverse operator...\n',mfilename);
  if parms.amd_inverse_flag
    [M,nnf] = ts_calc_inverse_operator_amd(G_xyz,parms.SNR,...
      parms.noisecovar,parms.srccovar);
  else
    [M,nnf] = ts_calc_inverse_operator(G_xyz,parms.SNR,...
      parms.noisecovar,parms.srccovar);
  end;
  if isempty(M)
    error('failed to calculate inverse operator');
  end;
  save(matfile, 'M', 'nnf', 'G_xyz');
else
  fprintf('%s: loading pre-calculated inverse operator...\n',mfilename);
  load(matfile);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply inverse to data

% initialize avg_data structure for least-squares (min-norm) fit
fprintf('%s: initializing best fit structure...\n',mfilename);
fit = [];
fit.num_sensors = avg_data.num_sensors;
fit.sfreq = avg_data.sfreq;
fit.sensor_info = avg_data.sensor_info;
fit.coor_trans = avg_data.coor_trans;
fit.noise = avg_data.noise;
for c=1:length(avg_data.averages)
  fit.averages(c).event_code=avg_data.averages(c).event_code;
  fit.averages(c).num_trials=avg_data.averages(c).num_trials;
  fit.averages(c).num_rejects=avg_data.averages(c).num_rejects;
  fit.averages(c).time=avg_data.averages(c).time;
  fit.averages(c).data=zeros(size(avg_data.averages(c).data));
  fit.averages(c).stdev=zeros(size(avg_data.averages(c).stdev));
end;

% initialize avg_data structure for residual error
fprintf('%s: initializing residual error structure...\n',mfilename);
err = [];
err.num_sensors = avg_data.num_sensors;
err.sfreq = avg_data.sfreq;
err.sensor_info = avg_data.sensor_info;
err.coor_trans = avg_data.coor_trans;
err.noise = avg_data.noise;
for c=1:length(avg_data.averages)
  err.averages(c).event_code=avg_data.averages(c).event_code;
  err.averages(c).num_trials=avg_data.averages(c).num_trials;
  err.averages(c).num_rejects=avg_data.averages(c).num_rejects;
  err.averages(c).time=avg_data.averages(c).time;
  err.averages(c).data=zeros(size(avg_data.averages(c).data));
  err.averages(c).stdev=zeros(size(avg_data.averages(c).stdev));
end;

[num_sensors,num_tpoints] = size(avg_data.averages(1).data);
scalefacts = zeros(parms.nchans,1);
scalefacts(parms.grad_chans)=parms.grad_scalefact;
scalefacts(parms.mag_chans)=parms.mag_scalefact;
scalefacts(parms.EEG_chans)=parms.EEG_scalefact;
scale_matrix = scalefacts*ones(1,num_tpoints);
for c=1:length(parms.conditions)
  cond = parms.conditions(c);
  % initialize variables
  S_xyz = zeros(num_tpoints,3*num_sources);
  S = zeros(num_tpoints,num_sources);
  F = zeros(num_tpoints,num_sources);
  Y = zeros(num_tpoints,num_sensors);
  Yfit = zeros(num_tpoints,num_sensors);
  E = zeros(num_tpoints,num_sensors);
  norm_var_E = zeros(num_tpoints,1);

  % apply inverse to data
  fprintf('%s: applying inverse for cond %d...\n',mfilename,cond);
  % first apply scale factors to data
  data = avg_data.averages(cond).data.*scale_matrix;
  Y = data(parms.goodchans,:)';
  S_xyz = (M*Y')';
  
  % calculate amplitudes and do noise normalization
  fprintf('%s: calculating source amplitudes...\n',...
    mfilename);
  % L should be number of trials in average, unless noise covariance
  %   calculated from avg data, in which case it should be 1
  if ismember(parms.ncov_type,[0,1]) || isempty(parms.noisecovar)
    L = 1;
  else
    L = avg_data.averages(cond).num_trials;
  end;
  nnf2L = (nnf.^2)/L;
  nnfL = nnf/L;
  for s=1:num_sources
    j = 3*(s-1) + 1;
    k = j + 2;
    if parms.signed_sources_flag && parms.orient_constr_flag
      S(:,s) = S_xyz(:,j);
      if parms.noisenorm_flag
        % noise normalization
        F(:,s) = S(:,s)/nnfL(j);
      end;
    else
      source = S_xyz(:,j:k);
      S(:,s) = sqrt(sum(source.^2,2));
      if parms.noisenorm_flag
        % noise normalization
        F(:,s) = (S(:,s).^2)/sum(nnf2L(j:k));
      end;
    end;
  end;

  % calculate fit and error
  fprintf('%s: calculating fit and residual error...\n',mfilename);
  % calculate separate error for each channel type
  Yfit = (G_xyz*S_xyz')';
  E = Y-Yfit;
  norm_var_E = 0;
  ntypes = 0;
  if parms.useEEG_flag
    [good_chans,ind,ind_good_eeg] = intersect(parms.EEG_chans,parms.goodchans);
    Yeeg = Y(:,ind_good_eeg);
    Yfiteeg = Yfit(:,ind_good_eeg);
    Eeeg = Yeeg - Yfiteeg;
    var_Eeeg = var(Eeeg,0,2);
    var_Yeeg = var(Yeeg,0,2);
    var_Yfiteeg = var(Yfiteeg,0,2);
    tmp_var_Yeeg = var_Yeeg;
    tmp_var_Yeeg(find(~var_Yeeg))=1;
    norm_var_Eeeg = var_Eeeg./var_Yeeg;
    norm_var_E = norm_var_E + norm_var_Eeeg;
    ntypes = ntypes + 1;
  else
    Yeeg = [];
    Yfiteeg = [];
    Eeeg = [];
    var_Eeeg = [];
    var_Yeeg = [];
    norm_var_Eeeg = [];
  end;
  if parms.usegrad_flag
    [good_chans,ind,ind_good_grad] = intersect(parms.grad_chans,parms.goodchans);
    Ygrad = Y(:,ind_good_grad);
    Yfitgrad = Yfit(:,ind_good_grad);
    Egrad = Ygrad - Yfitgrad;
    var_Egrad = var(Egrad,0,2);
    var_Ygrad = var(Ygrad,0,2);
    tmp_var_Ygrad = var_Ygrad;
    tmp_var_Ygrad(find(~var_Ygrad))=1;
    norm_var_Egrad = var_Egrad./var_Ygrad;
    norm_var_E = norm_var_E + norm_var_Egrad;
    ntypes = ntypes + 1;
  else
    Ygrad = [];
    Yfitgrad = [];
    Egrad = [];
    var_Egrad = [];
    var_Ygrad = [];
    norm_var_Egrad = [];
  end;
  if parms.usemag_flag
    [good_chans,ind,ind_good_mag] = intersect(parms.mag_chans,parms.goodchans);
    Ymag = Y(:,ind_good_mag);
    Yfitmag = Yfit(:,ind_good_mag);
    Emag = Ymag - Yfitmag;
    var_Emag = var(Emag,0,2);
    var_Ymag = var(Ymag,0,2);
    tmp_var_Ymag = var_Ymag;
    tmp_var_Ymag(find(~var_Ymag))=1;
    norm_var_Emag = var_Emag./var_Ymag;
    norm_var_E = norm_var_E + norm_var_Emag;
    ntypes = ntypes + 1;
  else
    Ymag = [];
    Yfitmag = [];
    Emag = [];
    var_Emag = [];
    var_Ymag = [];
    norm_var_Emag = [];
  end;
  norm_var_E = norm_var_E/ntypes;

  % add to fitted average and error structures
  fit.averages(cond).data(parms.goodchans,:)=Yfit';
  err.averages(cond).data(parms.goodchans,:)=E';

  matname = sprintf('%s/matfiles/%s_results_cond%02d.mat',parms.rootoutdir,parms.prefix,cond);
  if ~exist(matname,'file') || parms.overwrite_output_flag
    fprintf('%s: saving results to mat file %s...\n',mfilename,matname);
    save(matname,'F','S','S_xyz',...
      'Y', 'Yfit','E','norm_var_E',...
      'Yeeg','Yfiteeg','Eeeg','var_Eeeg','var_Yeeg','norm_var_Eeeg',...
      'Ygrad','Yfitgrad','Egrad','var_Egrad','var_Ygrad','norm_var_Egrad',...
      'Ymag','Yfitmag','Emag','var_Emag','var_Ymag','norm_var_Emag');
  end;
  if(parms.write_stc_flag)
    % write to lh stc file
    fname = sprintf('%s/stcfiles/%s_cond%02d-lh.stc',parms.rootoutdir,parms.prefix,cond);
    if ~exist(fname,'file') || parms.overwrite_output_flag
      fprintf('%s: saving to stc file %s...\n',mfilename,fname);
      lh_vnums=find(parms.lh_dec_dips==1)-1;
      num_lh_verts = length(lh_vnums);
      time = avg_data.averages(cond).time(1)*1000;
      % write_stc wants sfreq in Hz, t0 in ms
      if parms.signed_sources_flag && parms.orient_constr_flag &&...
         parms.noisenorm_flag
        tmp = F(:,1:num_lh_verts)'*parms.stc_scalefact;
      elseif parms.noisenorm_flag
        tmp = sqrt(F(:,1:num_lh_verts))'*parms.stc_scalefact;
      else
        tmp = S(:,1:num_lh_verts)'*parms.stc_scalefact;
      end;
      ts_write_stc(fname,lh_vnums,avg_data.sfreq,...
        time,tmp);
    end;

    % write to rh stc file
    fname = sprintf('%s/stcfiles/%s_cond%02d-rh.stc',parms.rootoutdir,parms.prefix,cond);
    if ~exist(fname,'file') || parms.overwrite_output_flag
      fprintf('%s: saving to stc file %s...\n',mfilename,fname);
      lh_vnums=find(parms.lh_dec_dips==1)-1;
      num_lh_verts = length(lh_vnums);
      rh_vnums=find(parms.rh_dec_dips==1)-1;
      time = avg_data.averages(cond).time(1)*1000;
      % write_stc wants sfreq in Hz, t0 in ms
      if parms.signed_sources_flag && parms.orient_constr_flag &&...
         parms.noisenorm_flag
        tmp = F(:,num_lh_verts+1:end)'*parms.stc_scalefact;
      elseif parms.noisenorm_flag
        tmp = sqrt(F(:,num_lh_verts+1:end))'*parms.stc_scalefact;
      else
        tmp = S(:,num_lh_verts+1:end)'*parms.stc_scalefact;
      end;
      ts_write_stc(fname,rh_vnums,avg_data.sfreq,...
        time,tmp);
    end;

    if(parms.write_mgh_flag)
      hemilist = {'lh' 'rh'};
      for h=1:length(hemilist)
        hemi = hemilist{h};
        stcfile = sprintf('%s/stcfiles/%s_cond%02d-%s.stc',...
          parms.rootoutdir,parms.prefix,cond,hemi);
        mghfile = sprintf('%s/mghfiles/%s_cond%02d-spsm%d-sm%d-%s.mgh',...
          parms.rootoutdir,parms.prefix,cond,...
          parms.sparsesmooth,parms.postsmooth,...
          hemi);
        if ~exist(mghfile,'file') || parms.overwrite_output_flag
          fprintf('%s: converting %s stc file to mgh file...\n',mfilename,hemi);
          ts_stc2mgh(stcfile,mghfile,parms.subjname,hemi,...
            parms.sparsesmooth,parms.postsmooth,...
            [],[],parms.subjdir,parms.mbmask_flag,parms.overwrite_output_flag);
        end;
        if(parms.resamp2ico_flag)
          infile = mghfile;
          outfile = fs_surf2ico(infile,parms.subjname,...
            'icolevel',parms.icolevel,'smooth_out',parms.icosmooth,...
            'subjdir',parms.subjdir,'forceflag',parms.overwrite_output_flag);
        end;
      end;
    end;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save fit and err sensor waveforms

matname=sprintf('%s/matfiles/%s_fiterr.mat',parms.rootoutdir,parms.prefix);
save(matname,'fit', 'err');

if(parms.write_fif_flag)
  outstem=sprintf('%s/fifs/%s_fit',parms.rootoutdir,parms.prefix);
  fprintf('%s: writing fitted average to fif...\n',mfilename);
  ts_avg2fif(fit,parms.template_fif,outstem,[],parms.overwrite_output_flag);

  outstem=sprintf('%s/fifs/%s_err',parms.rootoutdir,parms.prefix);
  fprintf('%s: writing residual error to fif...\n',mfilename);
  ts_avg2fif(err,parms.template_fif,outstem,[],parms.overwrite_output_flag);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('%s: finished.\n',mfilename);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
