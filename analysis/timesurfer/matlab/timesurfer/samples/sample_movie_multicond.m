% movie making script for group average, multiple views, multiple conditions

subj = 'fsaverage';
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

outstem = 'groupavg_multicond';

filelist = {...
  'groupavg_group_ucsd_cond_8_mean-lh.mgh'...
  'groupavg_group_ucsd_cond_8_mean-lh.mgh'...
  'groupavg_group_ucsd_cond_8_mean-rh.mgh'...
  'groupavg_group_ucsd_cond_8_mean-rh.mgh'...
  ...
  'groupavg_group_ucsd_cond_6_mean-lh.mgh'...
  'groupavg_group_ucsd_cond_6_mean-lh.mgh'...
  'groupavg_group_ucsd_cond_6_mean-rh.mgh'...
  'groupavg_group_ucsd_cond_6_mean-rh.mgh'...
  ...
  'groupavg_group_ucsd_cond_4_mean-lh.mgh'...
  'groupavg_group_ucsd_cond_4_mean-lh.mgh'...
  'groupavg_group_ucsd_cond_4_mean-rh.mgh'...
  'groupavg_group_ucsd_cond_4_mean-rh.mgh'...
};
            
hemilist =  {...
  'lh' 'lh' 'rh' 'rh'...
  'lh' 'lh' 'rh' 'rh'...
  'lh' 'lh' 'rh' 'rh'...
};

viewlist =  {...
  'lat' 'med' 'lat' 'med'...
  'lat' 'med' 'lat' 'med'...
  'lat' 'med' 'lat' 'med'...
};

labellist = {...
  'cond 8 lh lateral'     'cond 8 lh medial'     'cond 8 rh lateral'     'cond 8 rh medial'... 
  'cond 6 lh lateral'     'cond 6 lh medial'     'cond 6 rh lateral'     'cond 6 rh medial'... 
  'cond 4 lh lateral'     'cond 8vs6 lh medial'  'cond 8vs6 rh lateral'  'cond 8vs6 rh medial'... 
};

nrows = 3;
ncols = 4;

ts_surf2movie(...
  subj,...
  filelist,...
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


