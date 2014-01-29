function results = fs_aseg_roi(segname,funcname,roilist,roigrouplist,...
  maskname,minval)
%function results = fs_aseg_roi(segname,funcname,[roilist],[roigrouplist],...
% [maskname],[minval])
%
% Required Input:
%  segname: full path to segmentation volume (e.g. aparc+aseg.mgz)
%  funcname: full path to functional volume (mgh/mgz format)
%   (must be registered and resampled to T1)
%
% Optional Input:
%  roilist: (optional) a list of acceptable ROI codes
%   the ROI codes are the values in the segmentation volume
%   and are also found in the FreeSurferColorLUT.txt
%   If this is empty, only the ROI codes in the segmentation volume will be used
%   If an ROI code is specified but not found in the segmentation volume,
%     zero values will be returned for that ROI
%   If an ROI code is specified but not found in the LUT, that ROI code will be ignored
%  roigrouplist: (optional) a list of invidivual group ROI codes (e.g. 9998, 9999) in a cell var. 
%   it will only process these ROIs, and will include regardless if found or
%   not in FreeSurferColorLUT.txt
%  maskname: full path to mask volume (mgh/mgz format)
%  minval: minimum value
%   {default: 10^-5}
%
% created:  10/02/06 by Don Hagler
% last mod: 03/27/09 Don Hagler
%

%% todo: remove roigrouplist

if (~mmil_check_nargs(nargin,2)) return; end;

results = [];

if ~exist('roilist','var'), roilist=[]; end;
if ~exist('roigrouplist','var'), roigrouplist=[]; end;
if ~exist('maskname','var'), maskname=[]; end;
if ~exist('minval','var') || isempty(minval), minval = 10^-5; end;

if ~exist(segname,'file')
  error('segmentation volume file %s not found',segname);
end;
if ~exist(funcname,'file')
  error('functional volume file %s not found',funcname);
end;
if ~isempty(maskname) & ~exist(maskname,'file')
  error('mask volume file %s not found',maskname);
end;

% load segmentation volume
fprintf('%s: loading segmentation volume...\n',mfilename);
[segvol,M,mr_parms,seg_volsz] = fs_load_mgh(segname);

% load functional volume
fprintf('%s: loading functional volume...\n',mfilename);
[funcvol,M,mr_parms,func_volsz] = fs_load_mgh(funcname);
funcvol = reshape(double(funcvol),[prod(func_volsz(1:3)),func_volsz(4)]); % allow multiframe

if length(func_volsz)~=length(seg_volsz) | func_volsz~=seg_volsz
  error('dimensions of seg and func volumes do not match');
end;

% load mask volume
if ~isempty(maskname)
  fprintf('%s: loading mask volume...\n',mfilename);
  [maskvol,M,mr_parms,mask_volsz] = fs_load_mgh(maskname);
  maskvol = 1.0*(maskvol>minval);
  maskvol = reshape(maskvol,[prod(mask_volsz(1:3)),1]); % allow multiframe
  if any(func_volsz(1:3)~=mask_volsz(1:3))
    error('dimensions of func and mask volumes do not match');
  end;
  for f=1:func_volsz(4)
    funcvol(:,f) = funcvol(:,f).*maskvol;
  end;
end;

% load freesurfer LUT
fprintf('%s: loading freesurfer color lut...\n',mfilename);
[roicodes,roinames] = fs_read_fscolorlut;

% find unique roicodes from color lut
[roicodes,ind] = unique(roicodes);
roinames = roinames(ind,:);
nrois = length(roicodes);

% find unique roi numbers in segvol
segroicodes = unique(segvol(:));
segroicodes = segroicodes(find(segroicodes>0));
nsegrois = length(segroicodes);

if isempty(roilist), roilist = segroicodes; end;
[roilist,ia,ib] = intersect(roilist,roicodes);
roicodes=roicodes(ib);
roinames=roinames(ib,:);
fprintf('%s: extracting values...\n',mfilename);


if isempty(roigrouplist)
  nresults = length(roicodes);
else
  nresults = length(roigrouplist);
end;
for i=1:nresults
  invalid_flag = 0;
  if isempty(roigrouplist)
    results(i).roicode = roicodes(i);
    results(i).roiname = deblank(roinames(i,:));
    i_roi = find(segroicodes==results(i).roicode);
    roi = find(segvol==results(i).roicode);
    if isempty(i_roi), invalid_flag = 1; end;
  else
    results(i).roicode = roigrouplist(i).roi;
    results(i).roiname = deblank(roigrouplist(i).roiname);
    i_roi = intersect(segroicodes,roigrouplist(i).roi);
    roi = find(ismember(segvol(:),roigrouplist(i).roi));
    % if aseg file has 0 voxels for any ROI in roigrouplist, do not exclude
    missing_rois = setdiff(roigrouplist(i).roi,segroicodes); 
    if ~isempty(missing_rois)
      invalid_flag = 1; 
    end;
  end;
  
  if invalid_flag
    results(i).vals = [];
    results(i).nvox = 0;
    results(i).nvals = 0;
    results(i).avg = NaN;
    results(i).stdv = NaN;
  else
    raw_vals = funcvol(roi,:);
    ind_good_vals = find(abs(raw_vals(:,1))>minval & ~isnan(raw_vals(:,1)));
    vals = raw_vals(ind_good_vals,:);
    results(i).vals = vals;
    results(i).nvox = size(raw_vals,1);
    results(i).nvals = length(ind_good_vals);
    if results(i).nvals>0
        results(i).avg = mean(vals,1);
    else
        results(i).avg = NaN;
    end;
    if results(i).nvals>1
        results(i).stdv = std(vals,0,1);
    else
        results(i).stdv = NaN;
    end;
  end;
end

