function fs_normvol_asegroi(fname_in,fname_out,fname_aseg,roilist,forceflag)
%function fs_normvol_asegroi(fname_in,fname_out,fname_aseg,roilist,forceflag)
%
% Purpose: normalize a volume by the average value inside voxels specified
%  by aseg ROI numbers in a vector 'roilist'
%
% Required Options:
%   fname_in - file name of volume to be normalized
%   fname_out - output file name
%   fname_aseg - file name of freesurfer aseg volume
%     (full paths all)
%   roilist - vector of aseg ROI numbers
%   forceflag - whether to overwrite existing fname_out
%
% Created:  03/09/09 by Matt Erhart
% Last Mod: 03/12/10 by Don Hagler
%

if ~exist(fname_out) || forceflag
  if ~exist(fname_in), error('fname_in does not exist'); end;
  if ~exist(fname_aseg), error('aseg_in does not exist'); end;
  % load fname_in and fname_aseg
  [vol,M]=fs_load_mgh(fname_in);
  [vol_aseg,M_aseg]=fs_load_mgh(fname_aseg);
  %check that Vol sizes match
  if any(size(vol) ~= size(vol_aseg)), error('vol sizes different'); end;
  %check that Ms match
  if any(M~=M_aseg)
    error('M matrices for fname_in and fname_aseg are different');
  end;
  % find voxels in vol_aseg with values matching the values specified in roilist? (a vector of integers)
  k = find((ismember(vol_aseg, roilist)));
  % normalize vol by average of roi voxels
  vol = vol/mean(vol(k));
  % save normalized volume to fname_out
  fs_save_mgh(vol,fname_out,M);
end;

