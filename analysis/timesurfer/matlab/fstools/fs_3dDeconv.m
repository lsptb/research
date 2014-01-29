function [fname_out,fname_info] = fs_3dDeconv(fname,stim_fnames,varargin)
%function [fname_out,fname_info] = fs_3dDeconv(fname,stim_fnames,varargin)
%
% Purpose: wrapper around AFNI's 3dDeconvolve
%   input mgh file, get back mgh file as output
%
% Usage:
%  fs_3dDeconv(fname,'key1', value1,...); 
%
% Required Parameters:
%   fname: full or relative path name of input mgh file
%   stim_fnames: string or cell array of strings with full or relative path
%     names of "1D" stimulus time course files
%
% Optional parameters:
%  'fname_motion' - name of 1D file containing motion estimates from 3dvolreg
%    to be used as regressors
%    {default: []}
%  'stim_labels' - string or cell array of strings with stimulus condition names
%    If empty, will use labels such as "cond1", "cond2", etc.
%    {default: []}
%  'contrasts_flag' - [0|1] calculate glt contrasts between each condition
%    {default: 0}
%  'outdir' - output directory
%    If empty, will write output files to directory containing input fname
%    {default: []}
%  'skipTRs' - number of initial repetitions to omit from analysis
%    {default: 0}
%  'norm_flag' - [0|1] whether to normalize input timeseries
%    by mean for each voxel (new mean = 100) before doing Fourier calculations
%    {default: 1}
%  'detrend' - [0|1|2] whether and how to detrend input timeseries
%    0: no detrend (not recommended)
%    1: linear detrend
%    2: quadratic detrend
%    {default: 2}
%  'out_ext' - output file extension ('.mgh' or '.mgz')
%    if empty, will use input file extension
%    {default: []}
%  'forceflag' - [0|1] toggle overwrite existing output files
%    {default: 0}
%
% Output:
%   fname_out: output file name
%   fname_info: name of mat file containing stats_info struct array
%      with information about each frame in fname_out
%      including name, type, and dofs
%
%   Output file will be created with name based on input fname
%
% Created:  09/01/08 Don Hagler
% Last Mod: 07/03/09 Don Hagler
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fname_out = [];

if (~mmil_check_nargs(nargin, 2)), return; end;
parms = mmil_args2parms(varargin, { ...
  'stim_labels',[],[],...
  'contrasts_flag',false,[false true],...
  'outdir',[],[],...
  'skipTRs',0,[0 Inf],...
  'norm_flag',true,[false true],...
  'detrend',2,[0:2],...
  'thresh',10,[0 Inf],...
  'out_ext','.mgh',{'.mgh','.mgz'},...
  'stim_fnames',[],[],...
  'stim_labels',[],[],...
  'minlag',0,[0,10],...
  'maxlag',4,[0,30],...
  'fname_motion',[],[],...
  'glt_fnames',[],[],...
  'glt_labels',[],[],...
  'forceflag',false,[false true],...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check parameters

parms.fname_orig = fname;
hemi = [];
resamp_flag = 0;
if ~exist(fname,'file')
  error('file %s not found',fname);
end;
[fpath,fstem,fext] = fileparts(fname);
if ~ismember(fext,{'.mgh','.mgz'})
  error('input file must be mgh or mgz file type (has %s extension)',...
    fext);
end;
if isempty(parms.out_ext), parms.out_ext = fext; end;
if isempty(parms.outdir), parms.outdir = fpath; end;

if ~isempty(parms.fname_motion)
  nmotion = 6;
  motion_labels = {'Roll' 'Pitch' 'Yaw' 'dS' 'dL' 'dP'};
  if ~exist(parms.fname_motion,'file')
    error('file %s not found',parms.fname_motion);
  end;
else
  nmotion = 0;
  motion_labels = [];
end;

parms.stim_fnames = stim_fnames;
if ~isempty(parms.stim_fnames)
  if ~iscell(parms.stim_fnames), parms.stim_fnames = {parms.stim_fnames}; end;
  for c=1:length(parms.stim_fnames)
    if ~exist(parms.stim_fnames{c},'file')
      error('file %s not found',parms.stim_fnames{c});
    end;
  end;
  nconds = length(parms.stim_fnames);
  if isempty(parms.stim_labels)
    for c=1:nconds
      parms.stim_labels{c} = sprintf('cond%d',c);
    end;
  end;
else
  error('no stimulus files specified');
end;

nstims = nconds + nmotion; % motion regressors

if isempty(parms.glt_fnames)
  outstem = sprintf('%s/%s_3dDeconv',parms.outdir,fstem);
  [parms.glt_fnames,parms.glt_labels] = fs_create_glt_files(outstem,...
    'stim_labels',parms.stim_labels,...
    'contrasts_flag',parms.contrasts_flag,...
    'minlag',parms.minlag,...
    'maxlag',parms.maxlag,...
    'polort',parms.detrend,...
    'motion_flag',(nmotion>0),...
    'forceflag',parms.forceflag);
else
  if isempty(parms.glt_labels)
    for g=1:length(parms.glt_fnames)
      parms.glt_labels{g} = sprintf('glt%d',g);
    end;
  end;
end;
nglt = length(parms.glt_fnames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct output file names, check if they exist

if resamp_flag
  fname_out = sprintf('%s/%s_3dDeconv_resT1%s',...
    parms.outdir,fstem,parms.out_ext);
  fname_info = sprintf('%s/%s_3dDeconv_resT1.mat',...
    parms.outdir,fstem);
else
  fname_out = sprintf('%s/%s_3dDeconv%s',...
    parms.outdir,fstem,parms.out_ext);
  fname_info = sprintf('%s/%s_3dDeconv.mat',...
    parms.outdir,fstem);
end;
if exist(fname_out,'file') && exist(fname_info,'file') && ~parms.forceflag
  return;
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preprocess: skip TRs, normalize
if parms.skipTRs>0 || parms.norm_flag
  tmp_fname = fs_normdetrend(fname,'outdir',parms.outdir,...
    'skipTRs',parms.skipTRs,...
    'detrend',0,'norm_flag',parms.norm_flag,...
    'forceflag',parms.forceflag);
  fname = tmp_fname;
  [fpath,fstem,fext] = fileparts(fname);
end;

% convert to nii
tmp_fname = sprintf('%s/%s.nii',parms.outdir,fstem);
fs_mri_convert(fname,tmp_fname,'forceflag',parms.forceflag);

% check output files
tmp_fstem_out = sprintf('%s/%s_3dDeconv',parms.outdir,fstem);
tmp_fname_out = sprintf('%s+orig.BRIK',tmp_fstem_out);
tmp_fname_hdr = sprintf('%s+orig.HEAD',tmp_fstem_out);
if exist(tmp_fname_out,'file') || exist(tmp_fname_hdr,'file') && parms.forceflag
  delete(tmp_fname_out);
  delete(tmp_fname_hdr);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create shell script to run 3dDeconvolve
deconvscript = sprintf('%s/%s-deconv.sh',parms.outdir,fstem);
fid = fopen(deconvscript,'wt');
if fid==-1
  error('unable to open file %s for writing');
end;
fprintf(fid,'#!/bin/sh\n');
fprintf(fid,'3dDeconvolve \\\n');
fprintf(fid,'  -input %s \\\n',tmp_fname);
fprintf(fid,'  -polort %d \\\n',parms.detrend);
fprintf(fid,'  -num_stimts %d \\\n',nstims);
s=1;
for c=1:nconds
  fprintf(fid,'  -stim_file %d %s  -stim_label %d %s \\\n',...
    s,parms.stim_fnames{c},s,parms.stim_labels{c});
  fprintf(fid,'  -stim_minlag %d %d  -stim_maxlag %d %d \\\n',...
    s,parms.minlag,s,parms.maxlag);
  s=s+1;
end;
for m=1:nmotion
  fprintf(fid,'  -stim_file %d ''%s[%d]'' -stim_label %d %s -stim_base %d \\\n',...
    s,parms.fname_motion,m,s,motion_labels{m},s);
  s=s+1;
end;  
if nglt>0
  fprintf(fid,'  -num_glt %d',nglt);
  for g=1:nglt
    fprintf(fid,'  -glt 1 %s -glt_label %d "%s"\\\n',...
      parms.glt_fnames{g},g,parms.glt_labels{g});
  end;
end;
fprintf(fid,'  -nocout -fout -float \\\n');
fprintf(fid,'  -nox1D \\\n');
fprintf(fid,'  -bucket %s\n',tmp_fstem_out);
fclose(fid);

% run script
[status,result] = unix(sprintf('source %s',deconvscript));
if status
  error('failed to run shell script %s:\n%s',deconvscript,result);
else
  disp(result);
end;

% read HEAD file to get labels and degrees of freedom
stats_info = fs_read_HEAD_statsinfo(tmp_fname_out);
save(fname_info,'stats_info');

% convert output to mgh
fs_BRIK2mgh(tmp_fname_out,fname_out);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [glt_fnames,glt_labels] = fs_create_glt_files(outstem,varargin)
%function [glt_fnames,glt_labels] = fs_create_glt_files(outstem,[options])

  if (~mmil_check_nargs(nargin, 1)), return; end;
  parms = mmil_args2parms(varargin, { ...
    'nconds',1,[1,100],...
    'contrasts_flag',false,[false true],...
    'stim_labels',[],[],...
    'out_ext','.gltmat',[],...
    'minlag',0,[0,10],...
    'maxlag',4,[0,30],...
    'polort',1,[0,1,2],...
    'motion_flag',true,[false true],...
    'forceflag',false,[false true],...
  });

  glt_fnames = [];
  glt_labels = [];

  parms.nlags = parms.maxlag - parms.minlag + 1;
  parms.nregs = parms.nconds * parms.nlags;

  % prepare general linear tests for each condition
  %   (area under curve of hemodynamic response)
  if ~isempty(parms.stim_labels)
    if ~iscell(parms.stim_labels), parms.stim_labels = {parms.stim_labels}; end;
    parms.nconds = length(parms.stim_labels);
  end;
  
  g = 1;
  for c=1:parms.nconds
    if isempty(parms.stim_labels)
      stim_label = sprintf('cond%d',c);
    else
      stim_label = parms.stim_labels{c};
    end;
    fname_out = sprintf('%s_%s%s',outstem,stim_label,parms.out_ext);
    if ~exist(fname_out,'file') || parms.forceflag
      fid = fopen(fname_out,'wt');
      if fid==-1
        error('failed to open file %s for writing',fname_out);
      end;
      %% todo: if multiple runs, need separate baseline regs for each
      for p=0:parms.polort % number of baseline regs depends on polort (detrend)
        fprintf(fid,'0 ');
      end;
      % a '1' for each "lag" for this particular condition
      for tmpc=1:parms.nconds
        if tmpc==c
          val = 1;
        else
          val = 0;
        end;
        for lag = parms.minlag:parms.maxlag
          fprintf(fid,'%d ',val);
        end;
      end;
      % 6 extra zeros at end for motion regressors
      if parms.motion_flag
        for i=1:6
          fprintf(fid,'0 ');
        end;
      end;
      fprintf(fid,'\n');
      fclose(fid);
    end;
    glt_fnames{g} = fname_out;
    glt_labels{g} = stim_label;
    g=g+1;
  end;

  % prepare general linear tests for each contrast
  %   (one condition vs. another)
  if parms.contrasts_flag && parms.nconds>1
    for c1=1:parms.nconds
      if isempty(parms.stim_labels)
        stim_label1 = sprintf('cond%d',c1);
      else
        stim_label1 = parms.stim_labels{c1};
      end;
      for c2=1:parms.nconds
        if c1==c2, continue; end;
        if isempty(parms.stim_labels)
          stim_label2 = sprintf('cond%d',c2);
        else
          stim_label2 = parms.stim_labels{c2};
        end;
        fname_out = sprintf('%s_%sVS%s%s',outstem,...
          stim_label1,stim_label2,parms.out_ext);
        if ~exist(fname_out,'file') || parms.forceflag
          fid = fopen(fname_out,'wt');
          if fid==-1
            error('failed to open file %s for writing',fname_out);
          end;
          %% todo: if multiple runs, need separate baseline regs for each
          for p=0:parms.polort % number of baseline regs depends on polort (detrend)
            fprintf(fid,'0 ');
          end;
          % a '1' for each "lag" for this particular condition
          for tmpc=1:parms.nconds
            if tmpc==c1
              val = 1;
            elseif tmpc==c2
              val = -1;
            else
              val = 0;
            end;
            for lag = parms.minlag:parms.maxlag
              fprintf(fid,'%d ',val);
            end;
          end;
          % 6 extra zeros at end for motion regressors
          if parms.motion_flag
            for i=1:6
              fprintf(fid,'0 ');
            end;
          end;
          fprintf(fid,'\n');
          fclose(fid);
        end;
        glt_fnames{g} = fname_out;
        glt_labels{g} = sprintf('%sVS%s',stim_label1,stim_label2);
        g=g+1;
      end;
    end;
  end;
return;

