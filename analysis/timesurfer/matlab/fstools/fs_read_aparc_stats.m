function [lh_aparc_stats,rh_aparc_stats] = fs_read_aparc_stats(subj,subjdir,a2005_flag);
%function [lh_aparc_stats,rh_aparc_stats] = fs_read_aparc_stats(subj,[subjdir],[a2005_flag]);
%
% Required input:
%  subj is a string specifying the subject name
%
% Optional input:
%  subjdir - subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  a2005_flag - [0|1] toggle read aparc.a2005s.stats files instead of aparc.stats
%    {default = 0}
%
% Output:
%   lh_aparc_stats & rh_aparc_stats are struct arrays containing:
%     roiname
%     roicode
%     grayvol
%     surfarea
%     thickavg
%     thickstd
%     meancurv
%     gausscurv
%     foldind
%     curvind
%
% created:  04/01/09 by Don Hagler
% last mod: 05/15/09 by Don Hagler
%

%% todo: roigroups for lobar analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lh_aparc_stats = [];
rh_aparc_stats = [];
if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist('a2005_flag','var') | isempty(a2005_flag)
  a2005_flag = 0;
end;

if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as environment variable');
  end;
end;

if a2005_flag
  lh_aparc_fname = sprintf('%s/%s/stats/lh.aparc.a2005s.stats',subjdir,subj);
else
  lh_aparc_fname = sprintf('%s/%s/stats/lh.aparc.stats',subjdir,subj);
end;
if ~exist(lh_aparc_fname,'file')
  fprintf('%s: WARNING: lh aparc stats file %s not found\n',...
    mfilename,lh_aparc_fname);
  return;
end;

if a2005_flag
  rh_aparc_fname = sprintf('%s/%s/stats/rh.aparc.a2005s.stats',subjdir,subj);
else
  rh_aparc_fname = sprintf('%s/%s/stats/rh.aparc.stats',subjdir,subj);
end;
if ~exist(rh_aparc_fname,'file')
  fprintf('%s: WARNING: rh aparc stats file %s not found\n',...
    mfilename,rh_aparc_fname);
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read lh aparc file
if ~isempty(lh_aparc_fname)
  try
    fid = fopen(lh_aparc_fname);
    tmp_stats = textscan(fid,'%s %d %d %d %f %f %f %f %f %f\n',...
      'commentstyle','#');
    for i=1:length(tmp_stats{1})
      lh_aparc_stats(i).roiname   = char(tmp_stats{1}{i});
      lh_aparc_stats(i).surfarea  = double(tmp_stats{3}(i));
      lh_aparc_stats(i).grayvol   = double(tmp_stats{4}(i));
      lh_aparc_stats(i).thickavg  = double(tmp_stats{5}(i));
      lh_aparc_stats(i).thickstd  = double(tmp_stats{6}(i));
      lh_aparc_stats(i).meancurv  = double(tmp_stats{7}(i));
      lh_aparc_stats(i).gausscurv = double(tmp_stats{8}(i));
      lh_aparc_stats(i).foldind   = double(tmp_stats{9}(i));
      lh_aparc_stats(i).curvind   = double(tmp_stats{10}(i));
    end;
    % go back to beginning of file
    frewind(fid);
    % skip first 20 lines
    for t=1:21
      tmp = fgetl(fid);
    end;
    k = findstr(tmp,', ');
    if isempty(k)
      fprintf('%s: ERROR: unable to get Surface Area\n',mfilename);
    else
      i = i + 1;
      lh_aparc_stats(i).roiname = 'TotalSurfaceArea';
      tmp = tmp(k(3)+2:k(4)-1);
      lh_aparc_stats(i).surfarea = str2double(tmp);
      lh_aparc_stats(i).grayvol   = 0;
      lh_aparc_stats(i).thickavg  = 0;
      lh_aparc_stats(i).thickstd  = 0;
      lh_aparc_stats(i).meancurv  = 0;
      lh_aparc_stats(i).gausscurv = 0;
      lh_aparc_stats(i).foldind   = 0;
      lh_aparc_stats(i).curvind   = 0;
    end;
    fclose(fid);
  catch
    fprintf('%s: ERROR: failed to read lh aparc stats file\n',mfilename);
  end;
end;

  % read rh aparc file
if ~isempty(rh_aparc_fname)
  try
    fid = fopen(rh_aparc_fname);
    tmp_stats = textscan(fid,'%s %d %d %d %f %f %f %f %f %f\n',...
      'commentstyle','#');
    for i=1:length(tmp_stats{1})
      rh_aparc_stats(i).roiname   = char(tmp_stats{1}{i});
      rh_aparc_stats(i).surfarea  = double(tmp_stats{3}(i));
      rh_aparc_stats(i).grayvol   = double(tmp_stats{4}(i));
      rh_aparc_stats(i).thickavg  = double(tmp_stats{5}(i));
      rh_aparc_stats(i).thickstd  = double(tmp_stats{6}(i));
      rh_aparc_stats(i).meancurv  = double(tmp_stats{7}(i));
      rh_aparc_stats(i).gausscurv = double(tmp_stats{8}(i));
      rh_aparc_stats(i).foldind   = double(tmp_stats{9}(i));
      rh_aparc_stats(i).curvind   = double(tmp_stats{10}(i));
    end;
    % go back to beginning of file
    frewind(fid);
    % skip first 20 lines
    for t=1:21
      tmp = fgetl(fid);
    end;
    k = findstr(tmp,', ');
    if isempty(k)
      fprintf('%s: ERROR: unable to get Surface Area\n',mfilename);
    else
      i = i + 1;
      rh_aparc_stats(i).roiname = 'TotalSurfaceArea';
      tmp = tmp(k(3)+2:k(4)-1);
      rh_aparc_stats(i).surfarea = str2double(tmp);
      rh_aparc_stats(i).grayvol   = 0;
      rh_aparc_stats(i).thickavg  = 0;
      rh_aparc_stats(i).thickstd  = 0;
      rh_aparc_stats(i).meancurv  = 0;
      rh_aparc_stats(i).gausscurv = 0;
      rh_aparc_stats(i).foldind   = 0;
      rh_aparc_stats(i).curvind   = 0;
    end;
    fclose(fid);
  catch
    fprintf('%s: ERROR: failed to read rh aparc stats file\n',mfilename);
  end;
end;

