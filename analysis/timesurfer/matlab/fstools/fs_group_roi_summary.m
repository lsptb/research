function fs_group_roi_summary(fname_out,varargin);
%function fs_group_roi_summary(fname_out,[options]);
%
% Usage:
%  fs_group_roi_summary(fname_out,'key1', value1,...);
%
% Required input:
%  fname_out: output file name
%    This function will generate a csv (comma separated value) format
%    spreadsheet containing morphological data from Freesurfer aseg
%    and aparc ROIs for a given subjects directory
%
% Optional parameters:
%  'subjdir' - subjects directory (override SUBJECTS_DIR environment variable)
%    {default = $SUBJECTS_DIR}
%  'subjlist' - cell array of subject names
%    output will restricted only to those subjects (must exist in subjdir)
%    if empty or ommitted, all subjects in subjdir will be included
%    {default: []}
%  'grouplist' - cell array of group label for each subject
%    this list must be the same length as subjlist
%    this information will simply be included as an extra column
%    {default: []}
%  'a2005_flag' - [0|1] toggle use of a2005s.stats files
%    (otherwise use aparc.stats)
%    {default: 0}
%  'checkstatus_flag' - [0|1] toggle check that recon is complete
%    {default: 0}
%  'overwrite_flag' - [0|1] toggle overwrite existing output file
%    {default: 1}
%
% Parameters controlling what information is included in output file
%  'aseg' - [0|1] toggle whether to include aseg volumes
%    {default: 1}
%  'aparc' - [0|1] toggle whether to include aparc measures
%    if this value is 0, no aparc measures will be included and the
%      following 8 parameters will be ignored
%    {default: 1}
%  'grayvol' - [0|1] toggle whether to include gray matter volume
%    {default: 1}
%  'surfarea' - [0|1] toggle whether to include surface area
%    {default: 1}
%  'thickavg' - [0|1] toggle whether to include average thickness
%    {default: 1}
%  'thickstd' - [0|1] toggle whether to include standard deviation of thickness
%    {default: 1}
%  'meancurv' - [0|1] toggle whether to include mean curvature
%    {default: 1}
%  'gausscurv' - [0|1] toggle whether to include gaussian curvature
%    {default: 1}
%  'foldind' - [0|1] toggle whether to include folding index
%    {default: 1}
%  'curvind' - [0|1] toggle whether to include curvature index
%    {default: 1}
%
% created:  07/13/07 by Don Hagler
% last mod: 10/11/07 by Don Hagler
%

hemilist = {'Left','Right'};
aseg_roilist = [1:28,40:60,77:79,20001:20003];
exclude_aseg_roilist = [1,6,9,19:23,27,40,45,48,55,56,59];
aseg_roilist = setdiff(aseg_roilist,exclude_aseg_roilist);

num_aseg_rois = length(aseg_roilist);
num_aparc_rois = 36;
num_aparc_rois_2005 = 79;

%% todo: accept cell matrix of extra info


%% todo: keep example aseg and aparc stats in example matfile for writing header
%% todo: accept spreadsheet with subject info
%% todo: accept list of aseg roi numbers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parse options

if nargin<1
  help(mfilename);
  return;
end;

try
  options = varargin;
  for index = 1:length(options)
      if iscell(options{index}) & ~iscell(options{index}{1}), options{index} = { options{index} }; end;
  end;
  if ~isempty( varargin ), g=struct(options{:}); 
  else g = []; end;
catch
  fprintf('%s: calling convention {''key'', value, ... } error\n',mfilename);
  return;
end;

if ~isfield(g,'subjdir'), g.subjdir = []; end;
if ~isfield(g,'subjlist'), g.subjlist = []; end;
if ~isfield(g,'grouplist'), g.grouplist = []; end;
if ~isfield(g,'a2005_flag'), g.a2005_flag = 0; end;
if ~isfield(g,'overwrite_flag'), g.overwrite_flag = 1; end;
if ~isfield(g,'checkstatus_flag'), g.checkstatus_flag = 0; end;
if ~isfield(g,'aseg'), g.aseg = 1; end;
if ~isfield(g,'aparc'), g.aparc = 1; end;
if ~isfield(g,'grayvol'), g.grayvol = 1; end;
if ~isfield(g,'surfarea'), g.surfarea = 1; end;
if ~isfield(g,'thickavg'), g.thickavg = 1; end;
if ~isfield(g,'thickstd'), g.thickstd = 1; end;
if ~isfield(g,'meancurv'), g.meancurv = 1; end;
if ~isfield(g,'gausscurv'), g.gausscurv = 1; end;
if ~isfield(g,'foldind'), g.foldind = 1; end;
if ~isfield(g,'curvind'), g.curvind = 1; end;

gfields = fieldnames(g);
for index=1:length(gfields)
   switch gfields{index}
   case {'subjdir' 'subjlist' 'grouplist'...
    'a2005_flag' 'overwrite_flag' 'checkstatus_flag'...
    'aseg' 'aparc' 'grayvol' 'surfarea' 'thickavg' 'thickstd' 'meancurv',...
    'gausscurv' 'foldind' 'curvind'},;
   otherwise, error([mfilename ': unrecognized option: ''' gfields{index} '''' ]);
   end;
end;

% get rid of options struct
subjdir = g.subjdir;
subjlist = g.subjlist;
grouplist = g.grouplist;
a2005_flag = g.a2005_flag;
overwrite_flag = g.overwrite_flag;
checkstatus_flag = g.checkstatus_flag;
aseg = g.aseg;
aparc = g.aparc;
grayvol = g.grayvol;
surfarea = g.surfarea;
thickavg = g.thickavg;
thickstd = g.thickstd;
meancurv = g.meancurv;
gausscurv = g.gausscurv;
foldind = g.foldind;
curvind = g.curvind;
clear g;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parms

if isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    fprintf('%s: ERROR: SUBJECTS_DIR not defined as an environment variable\n',mfilename);
    return;
  end;
end;
if ~exist(subjdir,'dir')
  fprintf('%s: ERROR: subjects dir %s not found\n',mfilename,subjdir);
  return;
end;

if exist(fname_out,'file') & ~overwrite_flag
  fprintf('%s: output file %s already exists\n',...
    mfilename,fname_out);
  fprintf('      delete first or set overwrite_flag = 1\n');
  return;
end;

aparc_meas_list = {};
if aparc
  if grayvol, aparc_meas_list{end+1} = 'grayvol'; end;
  if surfarea, aparc_meas_list{end+1} = 'surfarea'; end;
  if thickavg, aparc_meas_list{end+1} = 'thickavg'; end;
  if thickstd, aparc_meas_list{end+1} = 'thickstd'; end;
  if meancurv, aparc_meas_list{end+1} = 'meancurv'; end;
  if gausscurv, aparc_meas_list{end+1} = 'gausscurv'; end;
  if foldind, aparc_meas_list{end+1} = 'foldind'; end;
  if curvind, aparc_meas_list{end+1} = 'curvind'; end;
  if isempty(aparc_meas_list), aparc = 0; end;
end;

if isempty(aseg) & isempty(aparc)
  fprintf('%s: ERROR: nothing selected for output\n',mfilename);
  return;
end;

if ~isempty(grouplist) & length(grouplist) ~= length(subjlist)
  fprintf('%s: ERROR: length of grouplist must match subjlist\n',mfilename);
  return;
end;

if a2005_flag
  num_aparc_rois = num_aparc_rois_2005;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid=fopen(fname_out,'wt');
if fid==-1
  fprintf('%s: ERROR: failed to open %s for writing\n',mfilename,fname_out);
  return;
end;
wrote_header_flag = 0;
dirlist = dir(sprintf('%s/*',subjdir));

tmp_subjlist = {dirlist.name};
if isempty(subjlist)
  subjlist = tmp_subjlist;
else
  [subjlist,i_subs] = intersect(subjlist,tmp_subjlist);
  if ~isempty(grouplist)
    grouplist = {grouplist{i_subs}};
  end;
end;
if isempty(subjlist)
  fprintf('%s: ERROR: subjlist contains no valid subjects\n',mfilename);
  return;
end;

for i=1:length(subjlist)
  subjname = subjlist{i};
  if ismember(subjname,{'.','..'}), continue; end;
  subjpath = sprintf('%s/%s',subjdir,subjname);
  if ~isdir(subjpath), continue; end;
  if checkstatus_flag
    [allexist,volexist] = fs_checktouchfiles(subjname,subjdir);
  else
    allexist=1;
    volexist=1;
  end;
  aseg_stats = []; lh_aparc_stats = []; rh_aparc_stats = [];
  if allexist | volexist
    % get cortical thickness, volume for aseg and aparc ROIs
    [aseg_stats,lh_aparc_stats,rh_aparc_stats] = ...
      fs_read_seg_stats(subjname,subjdir,a2005_flag);
    % reduce to subset of ROIs
    if ~isempty(aseg_stats)
      roicodes = cell2mat({aseg_stats.roicode});
      i_roicodes = find(ismember(roicodes,aseg_roilist));
      aseg_stats = aseg_stats(i_roicodes);    
    end;
  else
    fprintf('%s: WARNING: recon incomplete for %s\n',...
      mfilename,subjpath);
    continue;
  end;
  if isempty(aseg_stats) & isempty(lh_aparc_stats) & isempty(rh_aparc_stats)
    fprintf('%s: WARNING: unable to get stats for %s\n',...
      mfilename,subjpath);
    %% todo: output NaN's?
    continue;
  end;
  % check that extracted stats have correct number of ROIs, otherwise exclude from output
  %% todo: could output zeros, but then what if very first subject has wrong number of ROIs
  %%         (must be correct for header)
  if aseg
    tmp_num_aseg_rois = length(aseg_stats);
    if tmp_num_aseg_rois ~= num_aseg_rois
      fprintf('%s: WARNING: number of aseg ROIs (%d) for %s does not match expected (%d)\n',...
        mfilename,tmp_num_aseg_rois,subjname,num_aseg_rois);
      continue;
    end;
  end;
  if aparc
    tmp_num_aparc_rois = length(lh_aparc_stats);
    if tmp_num_aparc_rois ~= num_aparc_rois
      fprintf('%s: WARNING: number of lh aparc ROIs (%d) for %s does not match expected (%d)\n',...
        mfilename,tmp_num_aparc_rois,subjname,num_aparc_rois);
      continue;
    end;
    tmp_num_aparc_rois = length(rh_aparc_stats);
    if tmp_num_aparc_rois ~= num_aparc_rois
      fprintf('%s: WARNING: number of rh aparc ROIs (%d) for %s does not match expected (%d)\n',...
        mfilename,tmp_num_aparc_rois,subjname,num_aparc_rois);
      continue;
    end;
  end;
  if ~wrote_header_flag
    % write header row
    fprintf(fid,'"SubjectID"');
    if ~isempty(grouplist)
      fprintf(fid,',"GroupID"');
    end;
    if aseg
      for k=1:num_aseg_rois
        fprintf(fid,',"%s"',aseg_stats(k).roiname);
      end;
    end;
    if aparc
      for k=1:num_aparc_rois
        for h=1:length(hemilist)
          hemi = hemilist{h};
          switch hemi
            case 'Left'
              tmp_roiname = lh_aparc_stats(k).roiname;
            case 'Right'
              tmp_roiname = rh_aparc_stats(k).roiname;
          end;
          for m=1:length(aparc_meas_list)
            meas = aparc_meas_list{m};
            fprintf(fid,',"%s-ctx-%s-%s"',hemi,tmp_roiname,meas);
          end;
        end;
      end;
    end;
    fprintf(fid,'\n');
    wrote_header_flag=1;
  end;
  % write one row of values
  fprintf(fid,'"%s"',subjname);
  if ~isempty(grouplist)
    fprintf(fid,',"%s"',grouplist{i});
  end;
  if aseg
    for k=1:num_aseg_rois
      fprintf(fid,',%0.6f',aseg_stats(k).volume);
    end;
  end;
  if aparc
    for k=1:num_aparc_rois
      for h=1:length(hemilist)
        hemi = hemilist{h};
        switch hemi
          case 'Left'
            aparc_stats = lh_aparc_stats(k);
          case 'Right'
            aparc_stats = rh_aparc_stats(k);
        end;
        for m=1:length(aparc_meas_list)
          meas = aparc_meas_list{m};
          val = getfield(aparc_stats,meas);
          fprintf(fid,',%0.6f',val);
        end;
      end;
    end;
  end;
  fprintf(fid,'\n');
end;
fclose(fid);

