function mmil_aseg_brainmask(FSpath,varargin)
%function mmil_aseg_brainmask(FSpath,[options])
%
% Purpose: create a filled brainmask from a freesurfer aseg volume
%
% Usage: mmil_brainmask_from_aseg(datafile,'key1', value1,...);
%
% Required Input:
%  FSpath: full path of freesurfer recon
%
% Optional Input:
%  'fname_mask': output file name for brain mask
%    {default: 'brainmask.aseg.mgz'}
%  'brain_flag': [0|1] whether output should be a 
%                masked brain (1) or just a mask (0)
%    {default: 1}
%  'edit_flag': [0|1] whether to use aseg.mgz (0) or aseg_edit.mgz (1)
%    {default: 0}
%  'smooth1': smoothing sigma (voxels) for initial fill step
%    (set to 0 for no fill)
%    {default: 10}
%  'thresh1': threshold applied to mask after first smoothing step
%    {default: 0.05}
%  'smooth2': smoothing sigma (voxels) for second dilation step
%    (set to 0 for no dilation)
%    {default: 10}
%  'thresh2': threshold applied to mask after second smoothing step
%    (if smooth2=0, this is ignored)
%    {default: 0.4}
%  'smooth3': smoothing sigma (voxels) for third dilation step
%    (set to 0 for no dilation)
%    {default: 0}
%  'forceflag': [0|1] whether to overwrite fname_mask if it already exists
%    {default: 0}
%
% Created:  02/02/07 by Don Hagler
% Last Mod: 12/10/08 by Don Hagler
%
% NOTE: this is a new version of brainmask_from_aseg
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse input parameters

if (~mmil_check_nargs(nargin,1)) return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_mask','brainmask.aseg.mgz',[],...
  'brain_flag',true,[false true],...
  'edit_flag',false,[false true],...
  'smooth1',10,[0,100],...
  'thresh1',0.5,[0,1],...
  'smooth2',3,[0,100],...
  'thresh2',0.01,[0,1],...
  'smooth3',0,[0,100],...
  'forceflag',false,[false true],...
});

if isempty(parms.fname_mask)
  error('fname_mask is empty');
end;
if exist(parms.fname_mask,'file') && ~parms.forceflag
  return;
end;

if parms.edit_flag
  fname_aseg = sprintf('%s/mri/aseg_edit.mgz',FSpath);
else
  fname_aseg = sprintf('%s/mri/aseg.mgz',FSpath);
end;
if ~exist(fname_aseg,'file')
  error('file %s not found',fname_aseg);
end;

if parms.brain_flag
  fname_T1 =  sprintf('%s/mri/nu.mgz',FSpath);
  if ~exist(fname_T1,'file')
    error('file %s not found',fname_T1);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist(parms.fname_mask,'file') || parms.forceflag
  % load aseg
  [vol_aseg,mr_parms]=ctx_fs_load_mgh(fname_aseg);
  volmask = vol_aseg;
  % select brain
  volmask.imgs = 1.0*(vol_aseg.imgs>0 & vol_aseg.imgs<89);
  volmask.maxI = 1.0;
  vol_aseg = volmask;

  % smooth
  if parms.smooth1>0
    volmask.imgs = ...
      real(smooth3d(volmask.imgs,parms.smooth1,parms.smooth1,parms.smooth1));
    % truncate
    volmask.imgs = 1.0*(volmask.imgs>=parms.thresh1);
    % refresh from aseg
    volmask.imgs(vol_aseg.imgs>0)=1;
  end

  % smooth
  if parms.smooth2>0
    volmask.imgs = ...
      real(smooth3d(volmask.imgs,parms.smooth2,parms.smooth2,parms.smooth2));
    % truncate
    volmask.imgs = 1.0*(volmask.imgs>=parms.thresh2);
    % refresh from aseg
    volmask.imgs(vol_aseg.imgs>0)=1;
  end

  % smooth
  if parms.smooth3>0
    volmask.imgs = max(0,...
      real(smooth3d(volmask.imgs,parms.smooth3,parms.smooth3,parms.smooth3)));
  end;

  if parms.brain_flag
    % create masked T1 image
    [vol_T1,mr_parms] = ctx_fs_load_mgh(fname_T1);
    volmask.imgs = vol_T1.imgs.*volmask.imgs;
  end;
  ctx_fs_save_mgh(volmask,parms.fname_mask,mr_parms);
end;

