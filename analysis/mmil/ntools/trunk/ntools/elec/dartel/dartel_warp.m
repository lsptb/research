% List of open inputs
% Segment: Data - cfg_files
% Initial Import: Output Directory - cfg_files
% Deformations: Apply to - cfg_files
nrun = X; % enter the number of runs here
jobfile = {'/home/ccarlson/hugh/dartel/dartel_warp_job.m'};
jobs = repmat(jobfile, 1, nrun);
inputs = cell(3, nrun);
for crun = 1:nrun
    inputs{1, crun} = MATLAB_CODE_TO_FILL_INPUT; % Segment: Data - cfg_files
    inputs{2, crun} = MATLAB_CODE_TO_FILL_INPUT; % Initial Import: Output Directory - cfg_files
    inputs{3, crun} = MATLAB_CODE_TO_FILL_INPUT; % Deformations: Apply to - cfg_files
end
spm('defaults', 'FMRI');
spm_jobman('serial', jobs, '', inputs{:});
