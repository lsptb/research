function aseg_stats = fs_read_aseg_stats(subj,subjdir,roigroups);
%function aseg_stats = fs_read_aseg_stats(subj,[subjdir],[roigroups]);
%
% Required input:
%  subj: string specifying the subject name
%
% Optional input:
%  subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%  roigroups: struct array containing the following fields:
%    roiname - name of new ROI
%    roicode - new ROI code number
%    roicodes - vector of aseg ROI code numbers
%    {default: roigroups = fs_define_aseg_roigroups}
%               (includes 'WholeBrain' and 'AllVentricles')
%
% Output:
%   aseg_stats is a struct array containing:
%     roiname
%     roicode
%     volume
%
% created:  04/01/09 by Don Hagler
% last mod: 05/15/09 by Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aseg_stats = [];
if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    error('SUBJECTS_DIR not defined as environment variable');
  end;
end;

aseg_fname = sprintf('%s/%s/stats/aseg.stats',subjdir,subj);
if ~exist(aseg_fname,'file')
  fprintf('%s: WARNING: aseg stats file %s not found\n',mfilename,aseg_fname);
  return;
end;

if ~exist('roigroups','var') | isempty(roigroups)
  roigroups = fs_define_aseg_roigroups;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% read aseg file
if ~isempty(aseg_fname)
  try
    fid = fopen(aseg_fname);
    tmp_stats = textscan(fid,'%d %d %d %f %s %f %f %f %f %f\n',...
      'commentstyle','#');
    for i=1:length(tmp_stats{1})
      aseg_stats(i).roiname = char(tmp_stats{5}{i});
      aseg_stats(i).roicode = double(tmp_stats{2}(i));
      aseg_stats(i).volume  = double(tmp_stats{4}(i));
    end;

    % add extra ROIs from stats file  
    frewind(fid); % go back to beginning of file
    for t=1:17
      tmp = fgetl(fid); % skip first 16 lines
    end;
    k = findstr(tmp,', ');
    if isempty(k)
      fprintf('%s: WARNING: unable to get Brain Mask Volume\n',mfilename);
    else
      i = i + 1;
      aseg_stats(i).roiname = 'BrainMaskVolume';
      aseg_stats(i).roicode = 20001;
      tmp = tmp(k(3)+2:k(4)-1);
      aseg_stats(i).volume = str2double(tmp);
    end;
    tmp = fgetl(fid);
    tmp = fgetl(fid);
    k = findstr(tmp,', ');
    if isempty(k)
      fprintf('%s: WARNING: unable to get Brain Segmentation Volume\n',mfilename);
    else
      i = i + 1;
      aseg_stats(i).roiname = 'BrainSegmentationVolume';
      aseg_stats(i).roicode = 20002;
      tmp = tmp(k(3)+2:k(4)-1);
      aseg_stats(i).volume = str2double(tmp);
    end;
    tmp = fgetl(fid);
    k = findstr(tmp,', ');
    if isempty(k)
      fprintf('%s: WARNING: unable to get Intracranial Volume\n',mfilename);
    else
      i = i + 1;
      aseg_stats(i).roiname = 'IntracranialVolume';
      aseg_stats(i).roicode = 20003;
      tmp = tmp(k(3)+2:k(4)-1);
      aseg_stats(i).volume = str2double(tmp);
    end;
    fclose(fid);

    % add extra ROIs from roigroups
    all_codes =  cell2mat({aseg_stats.roicode});
    for i=1:length(roigroups)
      [roicodes,ind_group,ind_all] = find(ismember(roigroups(i).roicodes,all_codes));
      if isempty(ind_all)
        error('invalid roicodes in roigroups');
      else
        volume = 0;
        for j=ind_all
          volume = volume + aseg_stats(j).volume;
        end;
        aseg_stats(end+1).roiname = roigroups(i).roiname;
        aseg_stats(end).roicode = roigroups(i).roicode;
        aseg_stats(end).roicodes = roicodes;
        aseg_stats(end).volume = volume;
      end;
    end
  catch
    fprintf('%s: ERROR: failed to read aseg stats file\n',mfilename);
  end;
end;
