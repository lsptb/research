clear all; close all;

% Menu
a = menu('Auto Coregitration','MNI Space','Preop Space','Swap L/R','MGZ to ANALYZE');

fsl = '/space/monkeys/1/pubsw/packages/fsl/fsl-4.0.0_64/bin';
% coreg = '/space/md4/1/halgdev/analysis/iEEG_NYU';
coreg = '/home/ccarlson/loc/';
std = 'MNI152_T1_1mm.nii.gz';
stdpath = '/space/monkeys/1/pubsw/packages/fsl/fsl-4.0.0_64/data/standard/';

if a==1 || a==2
    
    t1 = menu('Select the pre-operation image','Pick up my own T1','Select Freesurfer Suject_Dir');
    if t1==1
        [preop preoppath] = uigetfile('*.*', 'Select the pre-operation image',coreg);
    else
        preoppath = uigetdir('/home/nyuproj/subjects/','Select Freesurfer Suject_Dir');
        preoppath = [preoppath '/mri/'];
        if exist([preoppath 'T1.nii.gz'],'file')
            preop = 'T1.nii.gz';
        elseif exist([preoppath 'T1.mgz'],'file')
            preop = 'T1.mgz';
        else
            error('No T1 detected, Please check the freesurfer folder!');
        end
    end
    % [std stdpath] = uigetfile(standard, 'Select the standard mni image');
    [elec elecpath] = uigetfile('*.*', 'Select the post-operation image with electrodes',coreg);
    
    if isnumeric(preop)||isnumeric(elec) 
        return; end
    
    elec_noext = elec(1:findstr(elec,'.')-1);
    preop_noext = preop(1:findstr(preop,'.')-1);
    
    % format convert from mgz to nii and swap x z -y and correct the
    % scanner ras center to the surface ras center
    elec_ext = elec(findstr(elec,'.')+1:end);
    preop_ext = preop(findstr(preop,'.')+1:end);
    if strcmp(elec_ext,'mgz')==1
        mri_elec = MRIread([elecpath elec]);
        convert = sprintf('mri_convert -ic %d %d %d %s%s %s%s.nii.gz',mri_elec.c_r, mri_elec.c_a,mri_elec.c_s,elecpath,elec,elecpath,elec_noext);
        swap_y = sprintf('setFSL414;fslswapdim %s%s.nii.gz x z -y %s%s.nii.gz',elecpath,elec_noext,elecpath,elec_noext);
        system(convert);
        system(swap_y);
        clear mri_elec
        elec = ([elec_noext '.nii.gz']);
        mri_elec = MRIread([elecpath elec]);
        mri_elec.vox2ras0(1,4) = mri_elec.vox2ras0(1,4)-mri_elec.c_r;
        mri_elec.vox2ras0(2,4) = mri_elec.vox2ras0(2,4)-mri_elec.c_a;
        mri_elec.vox2ras0(3,4) = mri_elec.vox2ras0(3,4)-mri_elec.c_s;
        MRIwrite(mri_elec,[elecpath elec]);
        elec_noext = elec(1:findstr(elec,'.')-1);
    end
    
    if strcmp(preop_ext,'mgz')==1
        mri_preop = MRIread([preoppath preop]);
        convert = sprintf('mri_convert -ic %d %d %d %s%s %s%s.nii.gz',mri_preop.c_r, mri_preop.c_a,mri_preop.c_s,preoppath,preop,elecpath,preop_noext);
        swap_lr = sprintf('setFSL414;fslswapdim %s%s.nii.gz x z -y %s%s.nii.gz',preoppath,preop_noext,elecpath,preop_noext);
        system(convert);
        system(swap_lr);
        clear mri_preop
        preop = ([preop_noext '.nii.gz']);
        mri_preop = MRIread([elecpath preop]);
        mri_preop.vox2ras0(1,4) = mri_preop.vox2ras0(1,4)-mri_preop.c_r;
        mri_preop.vox2ras0(2,4) = mri_preop.vox2ras0(2,4)-mri_preop.c_a;
        mri_preop.vox2ras0(3,4) = mri_preop.vox2ras0(3,4)-mri_preop.c_s;
        MRIwrite(mri_preop,[elecpath preop]);
        preop_noext = preop(1:findstr(preop,'.')-1);
    end
    
    if a==1   %% MNI Space

        % preop_mni
        ms_step1 = sprintf('setFSL414;flirt -in %s%s -ref %s%s -out %s%s_mni -omat %s%s_mni.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 12  -interp trilinear',preoppath,preop,stdpath,std,elecpath,preop_noext,elecpath,preop_noext);

        % elec_preop
        ms_step2 = sprintf('setFSL414;flirt -in %s%s -ref %s%s -out %s%s_preop -omat %s%s_preop.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear',elecpath,elec,preoppath,preop,elecpath,elec_noext,elecpath,elec_noext);

        % elec_preop_mni
        ms_step3 = sprintf('setFSL414;flirt -in %s%s_preop.nii.gz -applyxfm -init %s%s_mni.mat -out %s%s_preop_mni.nii.gz -paddingsize 0.0 -interp trilinear -ref %s%s_mni.nii.gz',elecpath,elec_noext,elecpath,preop_noext,elecpath,elec_noext,elecpath,preop_noext);

        % preop_mni_brain_mask
        ms_step4 = sprintf('setFSL414;bet %s%s_mni %s%s_mni_brain  -f 0.5 -g 0 -m',fsl,elecpath,preop_noext,elecpath,preop_noext);

        % multiply
        ms_step5 = sprintf('setFSL414;fslmaths %s%s_preop_mni.nii.gz -mul %s%s_mni_brain_mask.nii.gz %s%s_preop_mni_brain.nii.gz',elecpath,elec_noext,elecpath,preop_noext,elecpath,elec_noext);

        % excute
        system(ms_step1);
        system(ms_step2);
        system(ms_step3);
        system(ms_step4);
        system(ms_step5);

        % delete unnecessary files
        elec_preop = sprintf('%s%s_preop.nii.gz',elecpath,elec_noext);
        mask = sprintf('%s%s_mni_brain_mask.nii.gz',elecpath,preop_noext);
%         preop_mni_mat = sprintf('%s%s_mni.mat',elecpath,preop_noext);
%         elec_preop_mat = sprintf('%s%s_preop.mat',elecpath,elec_noext);

        delete(elec_preop);
        delete(mask);
%         delete(preop_mni_mat);
%         delete(elec_preop_mat);

        % show the result
        mricro_preop = sprintf('mricro %s%s_mni.nii.gz &',elecpath,preop_noext);
        mricro_elec_preop = sprintf('mricro %s%s_preop_mni.nii.gz &',elecpath,elec_noext);
        mricro_preop_brain = sprintf('mricro %s%s_mni_brain.nii.gz &',elecpath,preop_noext);
        mricro_elec_preop_brain = sprintf('mricro %s%s_preop_mni_brain.nii.gz &',elecpath,elec_noext);

        system(mricro_preop);
        system(mricro_elec_preop);
        system(mricro_preop_brain);
        system(mricro_elec_preop_brain);
    
    else     %% Preop Space

        % elec_preop
        ps_step1 = sprintf('setFSL414;flirt -in %s%s -ref %s%s -out %s%s_preop -omat %s%s_preop.mat -bins 256 -cost corratio -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -dof 6  -interp trilinear',elecpath,elec,preoppath,preop,elecpath,elec_noext,elecpath,elec_noext);

        % preop brain mask
        ps_step2 = sprintf('setFSL414;bet %s%s %s%s_brain  -f 0.5 -g 0 -m',preoppath,preop,elecpath,preop_noext);

        % multiply 
        ps_step3 = sprintf('setFSL414;fslmaths %s%s_preop.nii.gz -mul %s%s_brain_mask.nii.gz %s%s_preop_brain.nii.gz',elecpath,elec_noext,elecpath,preop_noext,elecpath,elec_noext);

        % excute
        system(ps_step1);
        system(ps_step2);
        system(ps_step3);

        % delete unnecessary files
        elec_preop_mat = sprintf('%s%s_preop.mat',elecpath,elec_noext);
        mask = sprintf('%s%s_brain_mask.nii.gz',elecpath,preop_noext);
        preop_brain = sprintf('%s%s_brain.nii.gz',elecpath,preop_noext);

%         delete(mask);
        delete(elec_preop_mat);
        delete(mask);
        delete(preop_brain);

        % show the result
%         mricro_preop = sprintf('mricro %s%s &',preoppath,preop);
        mricro_elec_preop = sprintf('mricro %s%s_preop.nii.gz &',elecpath,elec_noext);
%         mricro_preop_brain = sprintf('mricro %s%s_brain.nii.gz &',elecpath,preop_noext);
        mricro_elec_preop_brain = sprintf('mricro %s%s_preop_brain.nii.gz &',elecpath,elec_noext);

%         system(mricro_preop);
        system(mricro_elec_preop);
%         system(mricro_preop_brain);
        system(mricro_elec_preop_brain);
    
    end
    
elseif a==3     %% Swap L/R
    
    [tar tarpath] = uigetfile('*.*', 'Select the img file to swap',coreg);
    if isnumeric(tar) 
        return; end
    tar_noext = tar(1:findstr(tar,'.')-1);
    tar_ext = tar(findstr(tar,'.')+1:end);
    
    if strcmp(tar_ext,'mgz')==1
        convert = sprintf('mri_convert -ot nii %s%s %s%s.nii.gz',tarpath,tar,tarpath,tar_noext);
        swap_lr = sprintf('setFSL414;fslswapdim %s%s.nii.gz x z -y %s%s_swap_y.nii.gz',tarpath,tar_noext,tarpath,tar_noext);
        system(convert);
        system(swap_lr);
        tar = ([tar_noext '.nii.gz']);
    end
    
    swap = sprintf('setFSL414;fslswapdim %s%s -x y z %s%s_swap_x',tarpath,tar,tarpath,tar_noext);
    system(swap);
    
    % show the result
    mricro_swap = sprintf('mricro %s%s_swap_x.nii.gz &',tarpath,tar_noext);
    system(mricro_swap);
   
elseif a==4 % just convert from mgz to analyze
    [mgz mgzpath] = uigetfile('*.mgz', 'Select the mgz file to convert',coreg);
    if isnumeric(mgz) 
        return; end
    mgz_noext = mgz(1:findstr(mgz,'.')-1);
    convert = sprintf('mri_convert -ot nii %s%s %s%s.nii.gz',mgzpath,mgz,mgzpath,mgz_noext);
    swap_y = sprintf('setFSL414;fslswapdim %s%s.nii.gz x z -y %s%s_swap_y.nii.gz',mgzpath,mgz_noext,mgzpath,mgz_noext);
    system(convert);
    system(swap_y);
    
end

display('Co-regstration completed!');

