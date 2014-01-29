% movie making script for group average, multiple views

subj = 'fsaverage';
srcfile = 'groupavg_group_ucsd_cond_8_mean';
subjdir = '/space/md4/2/halgdev/analysis/MRI_UCSD/FWIO_MRI';
indir = '/space/emc2/1/halgdev/projects/FWIO/ucsd/MEG/groupavg';
outdir = pwd;

%tfirst = 275;
%tlast = 275;
tfirst = 175;
tlast = 320;
fthresh = 0;
fslope = 0.5;
fmid = 5.5;
image_size = '400x200';

outstem = 'groupavg_multiview';

hemilist =  {'lh' 'rh' 'lh' 'rh' 'lh' 'rh'};
viewlist =  {'lat' 'lat' 'med' 'med' 'ven' 'ven'};
labellist  =  {'lh lateral view' 'rh lateral view'...
               'lh medial view'  'rh medial view'...
               'lh ventral view' 'rh ventral view'};

nrows = 3;
ncols = 2;

ts_surf2movie(...
  subj,...
  srcfile,...
  'SUBJECTS_DIR',subjdir,...
  'label',labellist,...
  'hemi',hemilist,...
  'view',viewlist,...
  'indir',indir,...
  'outdir',outdir,...
  'image_size', image_size,...
  'outstem', outstem,...
  'sparsesmooth',0,...
  'postsmooth',0,...
  'fmid',fmid,...
  'fthresh',fthresh,...
  'fslope',fslope,...
  'tfirst',tfirst,...
  'tlast',tlast,...
  'nrows',nrows,'ncols',ncols,...
  'trim',1,...
  'hashnames',0, ...
  'timestamp','none'...
);

