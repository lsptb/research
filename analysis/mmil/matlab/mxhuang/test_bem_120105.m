% This script show and example how to calculate 1-shell BEM model
% (c) M.X. Huang, PhD DEC 1, 2005
% It was designed for UCSD MEG users only
%
%%% load the Vector-View sensor information into patient coordinate system
datafile='mxhuang_vectorview_ucsd_100305_som.fif';
[Coil_306,intpnt_loc_306,intpnt_ori_306]=sensor_306_normal(datafile);
% intpnt_loc_306 contains the locations of 816 integration points for 
%     306 channels: 2 per gradiometer, 4 per magnetometer
% intpnt_ori_306 contains the orientation of 816 integration points 
%     for 306 channels: 2 per gradiometer, 4 per magnetometer

%%% load the 7mm dipole grid from Freesurfer files
% left hemisphere dipoles in meters
lh_dip_name='mxh_anders_scan-7-lh.dip';
lh_dec_name='mxh_anders_scan-7-lh.dec';
[dip_info_lh,dec_dip_lh]=read_dipdec(lh_dip_name,lh_dec_name);

% right hemisphere dipoles in meters
rh_dip_name='mxh_anders_scan-7-rh.dip';
rh_dec_name='mxh_anders_scan-7-rh.dec';
[dip_info_rh,dec_dip_rh]=read_dipdec(rh_dip_name,rh_dec_name);

% now create the mesh grid, Freesurfer in mm
id_lh_dip=find(dec_dip_lh==1);
id_rh_dip=find(dec_dip_rh==1);
grid_mri=[dip_info_lh(1:3,id_lh_dip)';dip_info_rh(1:3,id_rh_dip)']; % in mm
n_grid=size(grid_mri,1);

%%% now load the inner skull tri file
tri_file_name='mxhuang_5mm_skull_mesh.tri';
[mesh_skull,geo_skull]=load_tri_file(tri_file_name);
temp=max(mesh_skull(:,1))-min(mesh_skull(:,1));
if temp > 100 % original unit in mm
    vertices=mesh_skull(:,1:3)/1000; % mm to m
elseif temp > 10 & temp < 30 % original unit in cm
    vertices=mesh_skull(:,1:3)/100; % cm to m
else
    vertices=mesh_skull(:,1:3); % in m
end
n_vert=size(vertices,1);
faces=geo_skull(:,1:3);

%%% load the .fif file for MRI-to-MEG (HEAD) transformation matrix,
%%% generated by "MRILAB" an Neuromag software after register the MEG to
%%% MRI. The matrix contains regid body transformation 
%%% (rotation and translation)
transfer_fif_name = 'COR-mxhuang_100305_som_meg_mri.fif';
T_mri_head=loadtrans(transfer_fif_name,'MRI','HEAD');

% mri to head for the dipole grid
grid_head=T_mri_head*[grid_mri'/1000;ones(1,n_grid)]; % mm to meters
grid_head=grid_head(1:3,:)'; %  in meter

% mri to head for BEM mesh vertices
vertices_head=T_mri_head*[vertices';ones(1,n_vert)];
vertices_head=vertices_head(1:3,:)'; 

%%% plot the dipole grid, BEM mesh, and MEG sensor intergration points
% now draw the BEM mesh, source nodes, and sensor locations
clf
bem_mesh_plot({vertices_head*1000},{faces})
plot3d(grid_head*1000,'b*');
plot3d(intpnt_loc_306*1000,'go')
title(['Sensor and source grid ' datafile]);
xlabel('x (mm)'); 
ylabel('y (mm)'); 
zlabel('z (mm)'); 
hold off
drawnow


%%% Now prepare for the BEM
% put together for BEM
bem_input.R_meg=intpnt_loc_306; % location for 816 integration points  
bem_input.O_meg=intpnt_ori_306; % orientation for 816 integration points
bem_input.R_eeg=[]; % no EEG
bem_input.vertices={vertices_head}; % location for BEM vertices
bem_input.faces={faces}; % face connection matrix for BEM mesh
bem_input.sigma=0.3; % conductivity SI unit
bem_input.mode=2; % MEG only
bem_input.basis_opt=1; % linear potential function
bem_input.test_opt=0; % collocation
bem_input.ISA=0; % Inhibit Isolated Skull Approach
bem_input.fn_eeg='bem_eeg_temp';
bem_input.fn_meg='bem_meg_temp';

% now calculate the BEM transfer matrix
return_var = bem_transfer_1shell(bem_input,[],1);

% BEM gain matrix (leadfield) for 816 integration points at the sensors
G_bem_xyz_inp=bem_gainmat(grid_head,return_var,0);

% finally convert 816 integration channels into 306 channels
G_bem_xyz=gain_intpnt2chan(Coil_306,G_bem_xyz_inp.meg);
% col #1,2,4,5,7,8... are for gradiometers fT/cm per nAm source
% col #3,6,9 are for magnetometers fT per nAm source
%
% The original units for MEG data are in T/m (gradiometers) 
%    and T (magnetometers)


