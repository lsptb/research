clear all; close all;

% get the subject info
Ssdir='/home/nyuproj/subjects';
[Subject_path] = uigetdir(Ssdir,'Select the subject FreeSurfer folder');
ii = findstr(Subject_path,'/');
Ssdir = Subject_path(1:ii(end));
Sdir = Subject_path(ii(end)+1:end);
Sub = regexp(Sdir,'NY+[0-9a-zA-Z]+_','match');
if ~isempty(Sub)
    ss2 = Sdir(length(char(Sub)):end);
else
    ss2 = Sdir;
end

%% read the inital text file
[FileName,PathName] = uigetfile('*.txt','Select the initial text file','/home/ccarlson/loc/');
jj = findstr(PathName,'/');
ss1 = PathName(jj(end-1)+1:jj(end)-1);
Sname = [ss1 ss2]

[dep_img_file dep_img_path] = uigetfile({'*.*';'*.img';'*.nii';'*.nii.gz'},...
        'Select the T1 preop image: ',PathName);

% hdr = load_nifti([dep_img_path dep_img_file],1);
% if ~isequal(hdr.pixdim(2),hdr.pixdim(3),hdr.pixdim(4))
%     disp('Image physical dimension not equal. Will affect the accuracy of distance calculation');
% end
% scale = hdr.pixdim(2);

hdr = MRIread([dep_img_path dep_img_file],1);
scale = hdr.volres(1);

fid = fopen([PathName, FileName]);
elec_all = textscan(fid,'%s %f %f %f %s');

elec_cell = [elec_all{1},num2cell(elec_all{2}),num2cell(elec_all{3}),num2cell(elec_all{4})];

if isempty(char(elec_all{5}(:)))
    g = strmatch('G',upper(elec_all{1}));
    d = strmatch('D',upper(elec_all{1}));
else
    g = strmatch('G',upper(elec_all{5}));
    d = strmatch('D',upper(elec_all{5}));
end

ini_grid = elec_cell(g,:);
ini_depth = elec_cell(d,:);
elec_cell([g;d],:) = [];
if cell2mat(elec_cell(:,2))>0
    sph_s = 'rh';
elseif cell2mat(elec_cell(:,2))<0
    sph_s = 'lh';
else
    sph_s = 'both'; % for x==0, make the elec in the major side
end

%%
% outer-brain surface check and create
ntools_elec_outer_brain(Subject_path)

% get the Grid
[elec_grid grid_stats] = ntools_elec_calc_grid(ini_grid,Subject_path,scale);

% get depth elecs
elec_depth = ntools_elec_calc_depth(ini_depth);

% get the strip
elec_strip = ntools_elec_calc_strip(elec_cell,Subject_path,sph_s);

% elec = [elec_grid; elec_strip; elec_depth];

% save the electrodes locations into a text file
fname_t1 = [PathName Sname '_coor_T1_' datestr(now,'mmddyy') '.txt'];
ntools_elec_savetxt(fname_t1,elec_grid,elec_strip,elec_depth);
[name x y z label] = textread(fname_t1,'%s %f %f %f %s');

% save all into binary image
fname_bin = [PathName,Sname,'_elec_bin_' datestr(now,'mmddyy'), '.nii.gz'];
elec_vox = ntools_elec_savebin([x y z],hdr,fname_bin);

% transform into mni space
elec_mni = ntools_elec_dartel_warp(fname_bin,[dep_img_path,dep_img_file]);
fname_mni = [PathName Sname '_coor_MNI_' datestr(now,'mmddyy') '.txt'];
ntools_elec_savetxt(fname_mni,[name num2cell(elec_mni) label]);

%% Save the surf.mat and plot
% setenv(Ssdir)
if strcmp(sph_s,'both')
    surf_brain_lh = fs_load_subj(Sdir,'lh','pial',0,Ssdir);
    if ~isfield(surf_brain_lh,'coords')
        surf_brain_lh.coords = surf_brain_lh.vertices;
    end
    surf_brain_rh = fs_load_subj(Sdir,'rh','pial',0,Ssdir);
    if ~isfield(surf_brain_rh,'coords')
        surf_brain_rh.coords = surf_brain_rh.vertices;
    end    
%     surf_brain_lh_smoothed = fs_load_subj(Sdir,'lh','pial-outer-smoothed',0,Ssdir);
%     surf_brain_rh_smoothed = fs_load_subj(Sdir,'rh','pial-outer-smoothed',0,Ssdir);
    lh_mat = [PathName,Sname,'_lh_pial_surf.mat'];
    rh_mat = [PathName,Sname,'_rh_pial_surf.mat'];
%     lh_mat_smoothed = [PathName,Sname,'_lh_pial_surf_smoothed.mat'];
%     rh_mat_smoothed = [PathName,Sname,'_rh_pial_surf_smoothed.mat'];
    save(lh_mat,'-struct','surf_brain_lh','coords','faces');
    save(rh_mat,'-struct','surf_brain_rh','coords','faces');
%     save(lh_mat_smoothed,'-struct','surf_brain_lh_smoothed','coords','faces');
%     save(rh_mat_smoothed,'-struct','surf_brain_rh_smoothed','coords','faces');
    ntools_elec_plot(fname_t1,lh_mat,rh_mat);
    clear surf_brain_lh surf_brain_rh surf_brain_lh_smoothed surf_brain_rh_smoothed
else
    surf_brain = fs_load_subj(Sdir,sph_s,'pial',0,Ssdir);
    if ~isfield(surf_brain,'coords')
        surf_brain.coords = surf_brain.vertices;
    end    
%     surf_brain_smoothed = fs_load_subj(Sdir,sph_s,'pial-outer-smoothed',0,Ssdir);
    surf_mat = [PathName,Sname,'_',sph_s,'_pial_surf.mat'];
%     surf_mat_smoothed = [PathName,Sname,'_',sph_s,'_pial_surf_smoothed.mat'];
    save(surf_mat,'-struct','surf_brain','coords','faces');
%     save(surf_mat_smoothed,'-struct','surf_brain_smoothed','coords','faces');
    ntools_elec_plot(fname_t1,surf_mat);
    clear surf_brain surf_brain_smoothed
end

save(['/home/ccarlson/hugh/NY_struct/',Sname,'_',datestr(now,'mmddyy')]);
