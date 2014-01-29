function results = fs_surf_roi_avg(fname_data,varargin);
%function results = fs_surf_roi_avg(fname_data,[options]);
%
% Purpose: extract average values or time courses for surface ROIs
%
% Usage:
%  fs_surf_roi_avg(fname_data,subject,'key1', value1,...); 
%
% Required Parameters:
%  fname_data: full path to surface data (mgh format)
%
% Optional Input (entered as 'key',value pairs):
%  'fname_label' full path to surface label file (freesurfer label file)
%    can be a cell array (curly bracketed list) of multiple label file names
%  'aparc' - full path to aparc annotation file
%    (e.g. subjdir/subj/label/lh.aparc.annot)
%  'frames' - for multi-frame data, data from only the specified frames
%     (e.g. time points) will be extracted
%     if empty or omitted, will use all available frames
%    {default = []}
%  'minval' - minimum value for inclusion in average
%     if empty or omitted, will accept any value
%    {default = []}
%
% Notes:
%  User must supply fname_label or aparc (or both)
%
%  For multi-frame data, results will contain matrices containing
%   data points for each frame (e.g. time courses)
%
% created:  07/14/08 by Don Hagler
% last mod: 07/14/08 by Don Hagler
%

%% NOTE: this replaces fs_aparc_roi.m
%% todo: test this script

results = [];
if (~mmil_check_nargs(nargin,1)), return; end;
parms = mmil_args2parms(varargin, { ...
  'fname_label',[],[],...
  'aparc',[],[],...
  'frames',[],[],...
  'minval',[],[],...
});

if isempty(parms.fname_label) && isempty(parms.aparc)
  error('must supply fname_label or aparc as input ROIs');
end;

% check input file
if ~exist(fname_data,'file')
  error('error surface data file %s not found\n',fname_data);
end;

% get volsz and check frames
[a,b,c,volsz] = fs_load_mgh(fname_data,[],[],1);
if volsz(2)==1 || volsz(3)==1
  error('input data must be on surface, not volume');
end;
if min(parms.frames)<1
  error('frame numbers must be > 1');
end;
if max(parms.frames)>volsz(4)
  error('frame numbers must be < number of frames (%d) in %s',volsz(4),fname_data);
end;

num_ROIs = 0;
ROI_vertices = {};
ROI_names = {};

% load aparc ROIs
if ~isempty(parms.aparc)
  fprintf('%s: loading aparc annotation...\n',mfilename);
  [aparc_roinums,aparc_names] = fs_read_annotation(aparcname);
  if volsz(1)~=length(aparc_roinums)
    error('number of data points in fname_data (%d) does not match aparc (%d)',...
      volsz(1),length(aparc_roinums));
  end;
  for i=1:length(aparc_names)
    num_ROIs = num_ROIs + 1;  
    ROI_vertices{num_ROIs} = find(aparc_roinums==i);
    ROI_names{num_ROIs} = aparc_names{i};
  end;
end;

% load label files
if ~isempty(parms.fname_label)
  if ~iscell(parms.fname_label), parms.fname_label = {parms.fname_label}; end;
  for i=1:length(parms.fname_label)
    num_ROIs = num_ROIs + 1;  
    fname_label = parms.fname_label{i};
    [tpath,tstem,text] = fileparts(fname_label);
    if ~strcmp(text,'.label')
      error('label file %s does not have .label extension',fname_label);
    end;
    ROI_names{num_ROIs} = tstem;
    if ~exist(fname_label,'file')
      error('label file %s not found',fname_label);
    end;
    ROI_vertices{num_ROIs} = fs_read_label(fname_label);
    if max(ROI_vertices{num_ROIs})>nverts | min(ROI_vertices{num_ROIs})<1
      error('bad label vertex number');
    end;
  end;
end;

% load surface data
fprintf('%s: loading functional surface data from %s...\n',mfilename,fname_data);
surfdata = squeeze(double(fs_load_mgh(fname_data,[],parms.frames)));
nverts = size(surfdata,1);
nframes = size(surfdata,2);

fprintf('%s: extracting values...\n',mfilename);
for i=1:length(ROI_vertices)
  results(i).roiname = ROI_names{i};
  roi = ROI_vertices{i};
  if isempty(roi)
    results(i).vals = [];
    results(i).nverts = 0;
    results(i).nvals = 0;
    results(i).nframes = 0;
    results(i).avg = NaN;
    results(i).stdv = NaN;
  else
    raw_vals = surfdata(roi,:);
    raw_vals(isnan(raw_vals)) = 0;
    if ~isempty(parms.minval)
      min_vals = min(abs(raw_vals),[],2);
      ind = find(min_vals>parms.minval);
      vals = raw_vals(ind,:);
    else
      vals = raw_vals;
    end;
    results(i).vals = vals;
    results(i).nverts = size(raw_vals,1);
    results(i).nvals = size(vals,1);
    if results(i).nvals>0
      results(i).avg = mean(vals,1);
    else
      results(i).avg = NaN;
    end;
    if results(i).nvals>1
      results(i).stdv = std(vals,1);
    else
      results(i).stdv = NaN;
    end;
  end;
end

