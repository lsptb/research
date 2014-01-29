function [G_xyz]=ts_calc_sph_gainmat(avg_data,lh_dip_info,rh_dip_info,...
  lh_dec_dips,rh_dec_dips,T_mri_head,cen_sph,radii,conductivities,...
  MEG_flag,EEG_flag);
%function [G_xyz]=ts_calc_sph_gainmat(avg_data,lh_dip_info,rh_dip_info,...
%  lh_dec_dips,rh_dec_dips,T_mri_head,cen_sph,radii,conductivities,...
%  MEG_flag,EEG_flag);
%
% Purpose: calculate gain matrix (forward solution) using spherical shell model
%
% Required Input:
%  avg_data: average data structure (see ts_avg_fif_data)
%  lh_dip_info: matrix containing locations and orientations for each dipole
%    in left hemisphere (6 rows, column for each vertex) generated with
%    ts_read_dip_file
%  rh_dip_info: matrix containing locations and orientations for each dipole
%    in right hemisphere (6 rows, column for each vertex)
%  lh_dec_dips: vector containing zeros and ones specifying which left
%    hemisphere dipoles to include can be generated with ts_read_dec_file
%  rh_dec_dips: vector containing zeros and ones specifying which right
%    hemisphere dipoles to include
%  T_mri_head: transformation matrix specifying registration between
%    mri (freesurfer brain space) and head (MEG/EEG sensor space)
%    can be generated with loadtrans function (fiff access) or with ts_pointreg
%  cen_sph: coordinates of center of sphere for spherical head model
%           supply as [x y z] (in mm)
%
% Optional Input:
%  radii: vector of three values specifying radius (in mm) of each shell:
%    Required if EEG_flag=1 -- not required if only MEG_flag=1 is given
%    (1) inner skull
%    (2) outer skull
%    (3) outer scalp
%    {default: []}
%  conductivities: vector of three values specifying conductivity of:
%    (1) brain
%    (2) skull
%    (3) scalp
%    {default: [0.3 0.01 0.3]}
%  MEG_flag: 0 or 1 -- indicating whether MEG gain matrix should be calculated
%    {default: 1}
%  EEG_flag: 0 or 1 -- indicating whether EEG gain matrix should be calculated
%    {default: 1}
%
%  NOTE -- either MEG_flag or EEG_flag must = 1 (both cannot be 0)
%
% Output:
%  G_xyz: gain matrix forward solution with three rows (x,y,z) for each dipole
%         units for gradiometers are fT/cm per nA source
%         units for magnetometers are fT per nA source
%         units for EEG sensors are uV per nA source
%
%         size of G_xyz is 3*ndips x nsensors
%         if EEG_flag=0, EEG columns of G_xyz will contain zeros
%         if MEG_flag=0, MEG columns of G_xyz will contain zeros
%         for non-MEG/EEG columns, G_xyz will contain zeros
%
% created 11/17/05 by Don Hagler (based on code from Mingxiong Huang)
% last modifed 10/05/06 DH
%

if nargin < 8
  help(mfilename);
  return;
end;

G_xyz = [];

scale_meg=1;
scale_eeg=0.001; % scaling factor converting EEG into uV by nAm dipole: 1e-9 (nAm) x 1e6 (uV)  

%%%%%%%%%%%%%%%%%%%%% initialize and check parms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('MEG_flag','var')
  MEG_flag = 1;
end;
if ~exist('EEG_flag','var')
  EEG_flag = 1;
end;
if ~exist('radii','var')
  radii = [];
end;
if ~exist('conductivities','var')
  conductivities = [0.3 0.01 0.3];
end;
if ~MEG_flag & ~EEG_flag
  fprintf('%s: error: both MEG_flag and EEG_flag cannot be zero\n',mfilename);
  return;
end;
nrad = length(radii);
if EEG_flag
  if nrad<3
    fprintf('%s: error: for EEG, must specify 3 radii\n',mfilename);
    return;
  end;
  if length(conductivities)~=nrad
    fprintf('%s: error: conductivities must have same number of values as radii\n',mfilename);
    return;
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% sensors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get integration points from sensor locs
fprintf('%s: preparing sensor info...\n',mfilename);
if MEG_flag
  MEG_info=ts_prep_MEG_info(avg_data);
end;
if EEG_flag
  EEG_info=ts_prep_EEG_info(avg_data);
end;
grad_chans = find(strncmp('grad',lower({avg_data.sensor_info.typestring}),...
             length('grad')));
mag_chans  = find(strcmp('mag',lower({avg_data.sensor_info.typestring})));
MEG_chans = union(grad_chans,mag_chans);
EEG_chans = find(strcmp('eeg',lower({avg_data.sensor_info.typestring})));
nsensors = length(avg_data.sensor_info);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get indices of select dipoles (non-zero in lh_dec_dips and rh_dec_dips)
id_lh_dip=find(lh_dec_dips);
id_rh_dip=find(rh_dec_dips);

% combine info for select dipoles
grid_mri=[lh_dip_info(1:3,id_lh_dip)';rh_dip_info(1:3,id_rh_dip)'];
n_grid=size(grid_mri,1);
fprintf('%s: total number of dipoles in forward = %d\n',mfilename,n_grid);

% apply coordinate transformation to get to sensor-space
grid_head=T_mri_head*[grid_mri'/1000;ones(1,n_grid)]; % in meters
grid_head=grid_head(1:3,:)'*1000; % back to mm and reshape
grid_head_cen=grid_head-ones(n_grid,1)*cen_sph;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% gain matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate gain matrix
if MEG_flag
  % calculate gain matrix for MEG
  fprintf('%s: calculating MEG gain matrix...\n',mfilename);
  Psensor(1).sensor = (MEG_info.intpnt_loc*100-ones(816,1)*cen_sph/10)'; % in cm
  Psensor(1).orient = MEG_info.intpnt_ori';
  Psensor(1).weight = 1000;
  G_xyz_inp=sarvas5(grid_head_cen'/10,Psensor,-1);  % fT and fT/cm

  % convert integration channels into channels
  fprintf('%s: converting integration points to channels...\n',mfilename);
  G_xyz_MEG=gain_intpnt2chan(MEG_info.Coil,G_xyz_inp);
end;
if EEG_flag
  % calculate gain matrix for EEG
  fprintf('%s: calculating EEG gain matrix...\n',mfilename);
  G_xyz_inp=gain3(grid_head_cen/1000,EEG_info.intpnt,...
    fliplr(radii/1000),fliplr(conductivities),1);

  % convert integration channels into channels
  fprintf('%s: converting integration points to channels...\n',mfilename);
  G_xyz_EEG=gain_intpnt2chan(EEG_info.sensor,G_xyz_inp);
end;

% combine MEG and EEG gain matrices
G_xyz=zeros(nsensors,3*n_grid);
if MEG_flag
  G_xyz(MEG_chans,:)=G_xyz_MEG*scale_meg;
end;
if EEG_flag
  G_xyz(EEG_chans,:)=G_xyz_EEG*scale_eeg;
end;

fprintf('%s: finished calculating gain matrix\n',mfilename);

return;
