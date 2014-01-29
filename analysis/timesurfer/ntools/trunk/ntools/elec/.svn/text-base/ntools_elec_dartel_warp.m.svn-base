function elec_mni = ntools_elec_dartel_warp(elec_vox,preop_t1)

fprintf('start DARTEL warping process......\n')

tic
[preop_t1_path,preop_t1_img,ext] = fileparts(preop_t1);
dartel_dir = [preop_t1_path '/dartel/'];

if ~exist(dartel_dir,'dir')
    mkdir(dartel_dir);
end

vox_img = gunzip(elec_vox,dartel_dir);

switch ext
    case '.gz'
        preop_t1 = gunzip([preop_t1_path '/' preop_t1_img ext],dartel_dir);
%         vox_img = gunzip(elec_vox,dartel_dir);
    case '.nii'
        copyfile([preop_t1_path '/' preop_t1_img ext],dartel_dir,'f');
        preop_t1 = cellstr([dartel_dir preop_t1_img ext]);
%         copyfile(elec_vox,dartel_dir,'f');
%         vox_img = cellstr([dartel_dir elec_vox_img]);
    otherwise
        elec_mni = [];
        disp('Please convert the Preop_T1 image into .nii or .nii.gz format first.')
        return
end

%% dartel warping

% cd(dartel_dir)

% List of open inputs
% Segment: Data - cfg_files
% Initial Import: Output Directory - cfg_files
% Deformations: Apply to - cfg_files
nrun = 1; % enter the number of runs here
jobfile = {'/home/ccarlson/hugh/dartel/dartel_warp_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);
for crun = 1:nrun
    inputs{1, crun} = preop_t1; % Segment: Data - cfg_files
    inputs{2, crun} = cellstr(dartel_dir); % Initial Import: Output Directory - cfg_files
    inputs{3, crun} = vox_img; % Deformations: Apply to - cfg_files
end
% spm('defaults', 'FMRI');
% spm_jobman('serial', jobs, '', inputs{:});
job_id = cfg_util('initjob',jobs);
sts = cfg_util('filljob',job_id,inputs{:});
if sts
    cfg_util('run',job_id);
end
cfg_util('deljob',job_id);

toc

%% get elecs' locations

[pathstr, name, ext] = fileparts(char(vox_img)); 
name = ['/w' name];
hdr_ch2 = load_nifti([pathstr, name, ext]);

% elec_list = unique(hdr_ch2.vol);
s = max(unique(hdr_ch2.vol));
elec_ch2_vox = zeros(s,3);

for ii=1:s
    ind = find(hdr_ch2.vol==ii);
    if ~isempty(ind)
        [a b c] = ind2sub(hdr_ch2.dim(2:4)',ind);
        a_avg = mean(a);
        b_avg = mean(b);
        c_avg = mean(c);
        elec_ch2_vox(ii,:) = [a_avg b_avg c_avg];
    else
        elec_ch2_vox(ii,:) = [0 0 0];
        disp(['Electrode #' num2str(ii) ' not found. Probably it shares the same position with anther electrode. Reset to [0 0 0].']);
    end
end

elec_ch2_vox = [elec_ch2_vox ones(s,1)]; 
elec_ch2_ras = hdr_ch2.vox2ras*elec_ch2_vox';
elec_mni = elec_ch2_ras';
elec_mni = elec_mni(:,1:3);
