function mmil_average_volumes(fname_inlist,fname_mask,fname_out,cleanupflag,forceflag)
%function mmil_average_volumes(fname_inlist,fname_mask,fname_out,[cleanupflag],[forceflag])
%
%  Purpose: Registering and averaging two or more MRI volumes
%            for a single subject, with the same contrast properties
%          e.g. averaging two T1-weighted images together
%
%  Required Input:
%    fname_inlist: cell array of input file names to be registered and averaged
%      input files must be mgh or mgz format
%    fname_mask: name of mgh file containing 3D binary mask registered to 
%      volume specified by first member of fname_inlist
%    fname_out: output file name
%      should have mgh or mgz extension
%
%  Optional Input:
%    cleanupflag: [0|1] whether to remove temporary files created by mriRegister
%      { default: 1 }
%    forceflag: [0|1] whether to run calculations even if output file exists
%      { default: 0 }
%
% Created:  07/09/08 by Don Hagler
% Last Mod: 07/09/08 by Don Hagler
%

if (~mmil_check_nargs(nargin,3)) return; end;
if ~exist('cleanupflag','var') || isempty(cleanupflag), cleanupflag = 1; end;
if ~exist('forceflag','var') || isempty(forceflag), forceflag = 0; end;

fname_inputparms = '/home/mmildev/bin/inputParamsRigid.txt'; % for mriRegister

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check input files

if ~iscell(fname_inlist) || length(fname_inlist)<2
  error('fname_inlist must contain at least two file names');
end;
for i=1:length(fname_inlist)
  fname = fname_inlist{i};
  if ~exist(fname,'file')
    error('input file %s not found',fname);
  end;
  [tmp_path,tmp_stem,tmp_ext]=fileparts(fname);
  if ~ismember(tmp_ext,{'.mgh','.mgz'})
    error('input files must be mgh or mgz format (%s)',fname);
  end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check mask file

if ~exist(fname_mask,'file')
  error('mask file %s not found',fname_mask);
end;
fname = fname_inlist{1};
[vol,M,mr_parms,volsz]=fs_load_mgh(fname,[],[],1);
[vol,M,mr_parms,volsz_mask]=fs_load_mgh(fname_mask,[],[],1);
if any(volsz~=volsz_mask)
  error('mask file %s image dimensions do not match first input file %s',fname_mask,fname);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check output file

[tmp_path,tmp_stem,tmp_ext]=fileparts(fname_out);
if ~ismember(tmp_ext,{'.mgh','.mgz'})
  fname = fname_out;
  fname_out = fullfile(tmp_path,tmp_stem,'.mgh');
  fprintf('%s: WARNING: changing fname_out %s to %s...',...
    mfilename,fname,fname_out);
end;

if exist(fname_out,'file') && ~forceflag, return; end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create temporary directory for mriRegister output

if isempty(tmp_path), tmp_path = pwd; end;
tmpdir = [tmp_path '/tmp_mriregister'];
if ~exist(tmpdir,'dir')
  [success,msg] = mkdir(tmpdir);
  if ~success, error('failed to create tmpdir %s:\n%s',tmpdir,msg); end;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% register volumes and average

fname_A = fname_inlist{1};
fprintf('%s: loading %s...\n',mfilename,fname_A);
[vol_sum,mr_parms] = ctx_fs_load_mgh(fname_inlist{i});
nvols = 1;
for i=2:length(fname_inlist)
  fname_B = fname_inlist{i};
  [tmp_path,tmp_stem,tmp_ext] = fileparts(fname_B);
  fname_tmp = sprintf('%s/%s_TP2_RigidReg_OrigTP1Scale%s',tmpdir,tmp_stem,tmp_ext);
  if ~exist(fname_tmp,'file') || forceflag
    fprintf('%s: registering %s to %s...\n',mfilename,fname_B,fname_A);
    tic
    cmd = sprintf('mriRegister -ip %s -s %s -sm %s -t %s -od %s',...
      fname_inputparms,fname_A,fname_mask,fname_B,tmpdir);
    disp(cmd);
    [status,result] = unix(cmd);
    toc
    disp(result);
    if status
      error('cmd %s failed',cmd);
    end;
    if ~exist(fname_tmp)
      error('file %s not created as expected',fname_tmp);
    end;
  end;
  fprintf('%s: loading %s...\n',mfilename,fname_tmp);
  vol = ctx_fs_load_mgh(fname_tmp);
  vol_sum.imgs = vol_sum.imgs + vol.imgs;
  nvols = nvols + 1;
end;

fprintf('%s: averaging...\n',mfilename);
vol_sum.imgs = vol_sum.imgs / nvols;
ctx_fs_save_mgh(vol,fname_out,mr_parms);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cleanup

if cleanupflag
  cmd = sprintf('rm -r %s',tmpdir);
  [status,result] = unix(cmd);
  if status, error('cmd %s failed:\n%s',cmd,status); end;
end;
