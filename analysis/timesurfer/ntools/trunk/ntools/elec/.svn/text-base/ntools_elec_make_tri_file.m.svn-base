function [name,bem_surf_path] = ntools_elec_make_tri_file(startdir)
% make BEM surfaces (.tri files) on PD image
bem = menu('Do you want to make .tri files based on PD or T1 image','PD','T1'); % PD=igor;T1=watershed
% startdir = '/home/halgdev/analysis/MRI_NYU/';

if bem==1
    
% select the PD image in nifti format
%infile = '/space/pogo2/2/mmildev-data/MMILDB/MEGFWD/Containers/MRIPROCESSED_NYn02_fmri_fastwords_20061122.152449.171000_1/loFA_res.mgh';

startdir2 = [startdir '/*'];
[FileName,PathName] = uigetfile(startdir2,'Select PD image');
infile=[PathName FileName];

FileNamex = FileName(1:findstr(FileName,'.')-1);

% % convert to mgh format and save in /bem folder
% str= sprintf('mri_convert %s %s/%s/bem/%s.mgh',nii_file,Ssdir,Sdir,FileNamex);
% system(str);

% outdir = [Ssdir '/' Sdir '/bem/'];
outdir = PathName;
outstem = FileNamex;
cmd = '/space/monkeys/1/home/mmildev/bin/find_bem_surfs';
surfflags = '';


default = menu('Use default values or enter values manually','manual','NYU defaults');

% if default ==1
%   site = menu('Where was the PD scan acquired?','NYU','UCSD');  
%    if site ==1 % NYU defaults 
%     inner_skull_cMRI = 0.8;
%      inner_skull_thresh = 200;
%    elseif site ==2 % UCSD defaults
%     inner_skull_cMRI = 0.7;
%     inner_skull_thresh = 70;   
%    end
% end

if default ==1
  inner_skull_cMRI = input('Enter value for inner_skull_cMRI:');
  inner_skull_thresh = input('Enter value for inner_skull_thresh:');
  str = sprintf('%s %s -inner_skull_nMove 300 -inner_skull_nRest 30 -inner_skull_nRelax 1 -inner_skull_cTang 0.03 -inner_skull_cNorm 0.03 -inner_skull_cMRI %g -inner_skull_cRepell 0.0 -inner_skull_thresh %g -outer_skull_nMove 150 -outer_skull_nRest 30 -outer_skull_nRelax 2 -outer_skull_cTang 0.15 -outer_skull_cNorm 0.15 -outer_skull_cMRI  0.8 -outer_skull_cRepell 1.0 -outer_skull_thresh 150.0 -outer_scalp_nMove 200 -outer_scalp_nRest 30 -outer_scalp_nRelax 1 -outer_scalp_cTang 0.02 -outer_scalp_cNorm 0.01 -outer_scalp_cMRI  0.8 -outer_scalp_cRepell 0.5 -outer_scalp_thresh 100.0 -outdir %s -outstem %s -ico 7 %s',cmd,infile,inner_skull_cMRI,inner_skull_thresh,outdir,outstem,surfflags);
  system(str);
end

% if default ~= 3
% str= sprintf('%s %s -inner_skull_nMove 300 -inner_skull_nRest 30 -inner_skull_nRelax 1 -inner_skull_cTang 0.03 -inner_skull_cNorm 0.03 -inner_skull_cMRI %g -inner_skull_cRepell 0.0 -inner_skull_thresh %g -outer_skull_nMove 150 -outer_skull_nRest 30 -outer_skull_nRelax 2 -outer_skull_cTang 0.15 -outer_skull_cNorm 0.15 -outer_skull_cMRI  0.8 -outer_skull_cRepell 1.0 -outer_skull_thresh 150.0 -outer_scalp_nMove 200 -outer_scalp_nRest 30 -outer_scalp_nRelax 1 -outer_scalp_cTang 0.02 -outer_scalp_cNorm 0.01 -outer_scalp_cMRI  0.8 -outer_scalp_cRepell 0.5 -outer_scalp_thresh 100.0 -outdir %s -outstem %s -ico 7 %s',cmd,infile,inner_skull_cMRI,inner_skull_thresh,outdir,outstem,surfflags)
% system(str);
% end

if default ==2

%NYU values

str= sprintf('%s %s -inner_skull_nMove 300 -inner_skull_nRest 30 -inner_skull_nRelax 1 -inner_skull_cTang 0.03 -inner_skull_cNorm 0.03 -inner_skull_cMRI 0.8 -inner_skull_cRepell 0.0 -inner_skull_thresh 200.0 -outer_skull_nMove 150 -outer_skull_nRest 30 -outer_skull_nRelax 2 -outer_skull_cTang 0.15 -outer_skull_cNorm 0.15 -outer_skull_cMRI  0.8 -outer_skull_cRepell 1.0 -outer_skull_thresh 150.0 -outer_scalp_nMove 200 -outer_scalp_nRest 30 -outer_scalp_nRelax 1 -outer_scalp_cTang 0.02 -outer_scalp_cNorm 0.01 -outer_scalp_cMRI  0.8 -outer_scalp_cRepell 0.5 -outer_scalp_thresh 100.0 -outdir %s -outstem %s -ico 6 %s',cmd,infile,outdir,outstem,surfflags);
system(str);
% else
% %UCSD values
% str= sprintf('%s %s -inner_skull_nMove 100 -inner_skull_nRest 30 -inner_skull_nRelax 1 -inner_skull_cTang 0.07 -inner_skull_cNorm 0.05 -inner_skull_cMRI 0.7 -inner_skull_cRepell 0.0 -inner_skull_thresh 70.0 -outdir %s -outstem %s -ico 4 %s',cmd,infile,outdir,outstem,surfflags)
% str= sprintf('%s %s -inner_skull_nMove 100 -inner_skull_nRest 30 -inner_skull_nRelax 1 -inner_skull_cTang 0.07 -inner_skull_cNorm 0.05 -inner_skull_cMRI 0.7 -inner_skull_cRepell 0.0 -inner_skull_thresh 70.0 -outdir %s -ico 4 %s',cmd,infile,outdir,surfflags)
% system(str);

end
display(pwd)

% change tri file name
str= sprintf('mv %s/%s_inner_skull6.tri %s/inner_skull6_ico10.tri',outdir,outstem,outdir);
system(str);
str= sprintf('mv %s/%s_outer_skull6.tri %s/outer_skull6_ico10.tri',outdir,outstem,outdir);
system(str);
str= sprintf('mv %s/%s_outer_scalp6.tri %s/outer_scalp6_ico10.tri',outdir,outstem,outdir);
system(str);

bem_surf_path = sprintf('%s/inner_skull6_ico10.tri',outdir);

% remove temp files
str= sprintf('rm -f %s/*.par',outdir);
system(str);
% quality check
str= sprintf('tkmedit -f %s -surface %s/inner_skull4.tri',infile,outdir);
system(str);

else
    
% [Ssdir] = uigetdir(startdir,'Select Freesurfer SUBJECTS_DIR directory');
% setenv('SUBJECTS_DIR', Ssdir);
% [Sdir_path] = uigetdir(startdir,'Select Freesurfer SUBJECT directory');
% Sdir = Sdir_path(length(Ssdir)+2:length(Sdir_path));
% setenv('SUBJECT', Sdir);

% str1= sprintf('cd %s/%s',Ssdir,Sdir);
% str2 = sprintf('mri_watershed -atlas -surf %s/%s/bem %s/%s/mri/T1.mgz %s/%s/mri/brainmask.auto.mgz',Ssdir,Sdir,Ssdir,Sdir,Ssdir,Sdir);
% str3 = sprintf('mv %s/%s/bem_brain_surface %s/%s/bem/inner_skull',Ssdir,Sdir,Ssdir,Sdir);
% str4 = sprintf('rm %s/%s/bem_inner_skull_surface',Ssdir,Sdir);
% str5 = sprintf('rm %s/%s/bem_outer_skin_surface',Ssdir,Sdir);
% str6 = sprintf('rm %s/%s/bem_outer_skull_surface',Ssdir,Sdir);

Sdir_path = startdir;
str1= sprintf('cd %s',Sdir_path);
str2 = sprintf('mri_watershed -atlas -surf %s/bem %s/mri/T1.mgz %s/mri/brainmask.auto.mgz',Sdir_path,Sdir_path,Sdir_path);
str3 = sprintf('mv -f %s/bem_brain_surface %s/bem/inner_skull',Sdir_path,Sdir_path);
str4 = sprintf('rm -f %s/bem_inner_skull_surface',Sdir_path);
str5 = sprintf('rm -f %s/bem_outer_skin_surface',Sdir_path);
str6 = sprintf('rm -f %s/bem_outer_skull_surface',Sdir_path);


 system(str1);
 system(str2);
 system(str3);
 system(str4);
 system(str5);
 system(str6);
 
% quality check
str= sprintf('tkmedit -f %s/mri/T1.mgz -surface %s/bem/inner_skull',Sdir_path,Sdir_path);
system(str);

bem_surf_path = sprintf('%s/bem/',Sdir_path);
name = 'inner_skull';

end

display('DONE');

return