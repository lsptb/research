function [elec_mni_ras tm] = ntools_elec_xyz2mni(elec_xyz,imgpath,img,savepath,savename)

% tranform the electrodes' positions from patients' space to the mni
% standard space and save into a binary image

fprintf('Converting electrodes into MNI space...');
tic;

fsl = '/usr/pubsw/packages/fsl/fsl-4.0.0_64/bin';
std = 'MNI152_T1_1mm.nii.gz';
stdpath = '/usr/pubsw/packages/fsl/fsl-4.0.0_64/data/standard/';

hdr = load_nifti([stdpath std],1);
Ssname = [savepath savename];
preop_mni = sprintf('%s/flirt -in %s%s -ref %s%s -omat %s_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear',...
    fsl,imgpath,img,stdpath,std,Ssname);
% apply trans matrix on elec_bin
% mni_bin = sprintf('%s/flirt -in %s_elec_bin.nii.gz -applyxfm -init %s_mni.mat -out %s_elec_bin_mni.nii.gz -paddingsize 0.0 -interp trilinear -ref %s%s',...
%     fsl,Ssname,Ssname,Ssname,stdpath,std);

system(preop_mni);
% system(mni_bin);

 tm = load([Ssname '_mni.mat'],'-ascii');
 
 elec_vox = [elec_xyz ones(length(elec_xyz),1)]';
 elec_mni_vox = tm*elec_vox;
%  elec_mni = elec_mni_vox(1:3,:)';
%  for i=1:3
%     elec_mni(:,i) = (elec_mni(i,:)./elec_mni(4,:))';
%  end

elec_mni_ras = hdr.vox2ras*elec_mni_vox;
elec_mni_ras = elec_mni_ras(1:3,:)';

% ntools_elec_savebin(elec_mni_ras,hdr,[Ssname '_elec_bin_mni_' datestr(now,'mmddyy') '.nii.gz']);

delete([Ssname '_mni.mat']);

 fprintf('Done. (%f seconds) \n\n', toc);
