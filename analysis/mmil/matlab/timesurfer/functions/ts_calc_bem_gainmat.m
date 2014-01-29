function [G_xyz]=ts_calc_bem_gainmat(avg_data,lh_dip_info,rh_dip_info,...
  lh_dec_dips,rh_dec_dips,T_mri_head,bem_matfile,bem_surf_files,...
  conductivities,MEG_flag,EEG_flag)
%function [G_xyz]=ts_calc_bem_gainmat(avg_data,lh_dip_info,rh_dip_info,...
%  lh_dec_dips,rh_dec_dips,T_mri_head,bem_matfile,bem_surf_files,...
%  conductivities,MEG_flag,EEG_flag)
%
% Purpose: calculate gain matrix (forward solution) using boundary element method
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
%  bem_matfile: name of matfile to be created for storage of bem calculations
%  bem_surf_files: cell array of file names for FreeSurfer tri files (text files
%   specifying triangles and faces) with MRI-derived skull and scalp surfaces
%   this should contain three filenames:
%      (1) inner skull
%      (2) outer skull
%      (3) outer scalp
%   Optionally, for single-shell MEG-only gain matrix, a single filename
%     for the inner skull surface can be specified here.
%
% Optional Input:
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
% last modifed 08/09/06 DH
%

if nargin < 8
  help(mfilename);
  return;
end;

G_xyz = [];

%%%%%%%%%%%%%%%%%%%%% initialize and check parms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('MEG_flag','var')
  MEG_flag = 1;
end;
if ~exist('EEG_flag','var')
  EEG_flag = 1;
end;
if ~exist('conductivities','var')
  conductivities = [0.3 0.01 0.3];
end;
if ~MEG_flag & ~EEG_flag
  fprintf('%s: error: both MEG_flag and EEG_flag cannot be zero\n',mfilename);
  return;
end;
if isempty(bem_surf_files)
  fprintf('%s: error: bem_surf_files is empty\n',mfilename);
  return;
end;
if iscell(bem_surf_files)
  num_bem_surfs=length(bem_surf_files);
else
  bem_surf_files={bem_surf_files};
  num_bem_surfs=1;
end;
if EEG_flag & num_bem_surfs<3
  fprintf('%s: error: for EEG, must specify 3 bem surf files\n',mfilename);
  return;
end;
if length(conductivities)~=num_bem_surfs
  fprintf('%s: error: conductivities must same number of values as num bem surfs\n',mfilename);
  return;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bem surfaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load bem surfaces
[verts,faces]=ts_load_bem_surfs(bem_surf_files);
for s=1:num_bem_surfs
  tmp_verts=verts{s};
  maxdist=max(tmp_verts(:,1))-min(tmp_verts(:,1));
  if maxdist > 100 % original unit in mm
      tmp_verts=tmp_verts(:,1:3)/1000; % mm to m
  elseif maxdist > 10 & maxdist < 30 % original unit in cm
      tmp_verts=tmp_verts(:,1:3)/100; % cm to m
  else
      tmp_verts=tmp_verts(:,1:3); % in m
  end
  nvert=size(tmp_verts,1);
  % apply coordinate transformation to get to "head"-space
  tmp_verts=(T_mri_head*[tmp_verts';ones(1,nvert)])';
  verts{s}=tmp_verts(:,1:3);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get indices of select dipoles (non-zero in lh_dec_dips and rh_dec_dips)
id_lh_dip=find(lh_dec_dips==1);
id_rh_dip=find(rh_dec_dips==1);

% combine info for select dipoles
grid_mri=[lh_dip_info(1:3,id_lh_dip)';rh_dip_info(1:3,id_rh_dip)'];
n_grid=size(grid_mri,1);
fprintf('%s: total number of dipoles in forward = %d\n',mfilename,n_grid);

% apply mri2head coordinate transformation to get to "head"-space
grid_head=T_mri_head*[grid_mri'/1000;ones(1,n_grid)]; % mm to meters
grid_head=grid_head(1:3,:)'; % reshape, still in meters
source_locs = grid_head;

%%%%%%%%%%%%%%%%%%%%%%%%%%%% BEM gain matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set BEM input
if num_bem_surfs==1
  sigma=conductivities(1);
else
  sigma=conductivities;
end;
if MEG_flag
  bem_input.R_meg=MEG_info.intpnt_loc; % location of MEG integration points  
  bem_input.O_meg=MEG_info.intpnt_ori; % orientation of MEG integration points
else
  bem_input.R_meg=[];
  bem_input.O_meg=[];
end;
if EEG_flag
  bem_input.R_eeg=EEG_info.intpnt; % location of EEG integration points
else
  bem_input.R_eeg=[];
end;
if MEG_flag & EEG_flag
  bem_input.mode=3; % MEG and EEG
  fprintf('%s: will calculate BEM for MEG and EEG\n',mfilename);
elseif MEG_flag
  bem_input.mode=2; % MEG only
  fprintf('%s: will calculate BEM for MEG only\n',mfilename);
elseif EEG_flag
  bem_input.mode=1; % EEG only
  fprintf('%s: will calculate BEM for EEG only\n',mfilename);
end;

bem_input.vertices=verts; % location for BEM vertices
bem_input.faces=faces; % face connection matrix for BEM mesh
bem_input.sigma=sigma; % conductivity SI unit
bem_input.basis_opt=1; % linear potential function
bem_input.test_opt=0; % collocation
bem_input.ISA=0; % Inhibit Isolated Skull Approach
bem_input.fn_eeg='bem_eeg_transfer';
bem_input.fn_meg='bem_meg_transfer';

% calculate BEM transfer matrix
fprintf('%s: calculating BEM transfer matrix...\n',mfilename);
tic
if num_bem_surfs==1
  return_var = bem_transfer_1shell(bem_input,[],bem_matfile,1);
else
  return_var = bem_transfer_3shell(bem_input,[],bem_matfile,1);
end;
toc

% BEM gain matrix (leadfield) for integration points at the sensors
fprintf('%s: calculating BEM gain matrix...\n',mfilename);
tic
G_bem_xyz_inp=bem_gainmat(grid_head,return_var,0);
toc

% convert integration channels into channels
fprintf('%s: converting integration points to channels...\n',mfilename);
if MEG_flag
  G_xyz_MEG=gain_intpnt2chan(MEG_info.Coil,G_bem_xyz_inp.meg);
end;
if EEG_flag
  G_xyz_EEG=gain_intpnt2chan(EEG_info.sensor,G_bem_xyz_inp.eeg);
end;

% combine MEG and EEG gain matrices
G_xyz=zeros(nsensors,3*n_grid);
if MEG_flag
  G_xyz(MEG_chans,:)=G_xyz_MEG;
end;
if EEG_flag
  % need to reverse polarity of EEG bem gain matrix for some reason (hence negative)
  G_xyz(EEG_chans,:)=-G_xyz_EEG;
end;

fprintf('%s: finished calculating gain matrix\n',mfilename);

return;
