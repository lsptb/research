function ts_max_n_avg_stc2mgh(stcfile,mghfile,subj,hemi,matfile,min_lat,max_lat,sparsesmooth,postsmooth,t0,t1,subjdir,mbmask_flag)
% function ts_max_n_avg_stc2mgh(stcfile,mghfile,subj,hemi,matfile,[min_lat],[max_lat],[sparsesmooth],[postsmooth],[t0],[t1],[subjdir],[mbmask_flag])
% % based on ts_stc2mgh by Don Hagler (see below)
% created by Nima Dehghani July 26th 2007
% 
%
% 
% converts a stc (source time course) file to mgh format (readable by tksurfer): solution time course of the mgh file will
% include two arrays. The first time includes max of STC of the desired interval and the second one includes average of dipole
% strengths of the same interval
% optionally applies smoothing on surface
%
%  Required input:
%    stcfile: full pathname of input stc file
%    mghfile: full pathname of output mgh file
%    subj:    subject name
%    hemi:    cortical hemisphere
%    matfile: mat file that includes output_data.epochs.time struc
%    min_lat: minimum latency of desired range for creating the average and finding the max STC (in ms)
%    max_lat: in ms%    
%    
%  Optional parameters:
%    sparsesmooth: number of sparse smoothing steps
%      {default = 0}
%      note: sparse smoothing is a fast way to fill
%            in gaps between sparse vertices
%    postsmooth: number of normal smoothing steps
%      {default = 0}
%      note: postsmoothing is additional nearest-neighbor average
%            smoothing applied after sparse smoothing
%    t0: first time sample to extract
%      {default = 1}
%    t1: last time sample to extract
%      {default = last}
%   subjdir: subjects directory (override SUBJECTS_DIR environment variable)
%    subjdir/subj should contain the freesurfer subject directory
%    {default = $SUBJECTS_DIR}
%   mbmask_flag: [0|1] toggle mask out thalamus, corpus callosum
%    {default: 0}
% NOTICE that mbmask_flag is 0 by DEFAULT (in contrast to ts_stc2mgh)
%
% Created by Nima Dehghani 26th jul 2007
%
%
% 
% see related ts_stc2mgh help
% function ts_stc2mgh(stcfile,mghfile,subj,hemi,[sparsesmooth],[postsmooth],[t0],[t1],[subjdir],[mbmask_flag])
%
% ts_stc2mgh created:       <01/01/07 Don Hagler
% ts_stc2mgh last modified:  04/14/07 Don Hagler
%
% ts_max_n_avg_stc2mgh created : Nima Dehghani 26th july 2007
%%%%

if nargin < 4
  help(mfilename);
  return;
end;

if ~exist('sparsesmooth','var') | isempty(sparsesmooth), sparsesmooth=0; end;
if ~exist('postsmooth','var') | isempty(postsmooth), postsmooth=0; end;

if ~exist('t0','var'), t0=[]; end;
if isempty(t0), t0=0; end;
if ~exist('t1','var'), t1=[]; end;
if isempty(t1), t1=0; end;

if ~exist('subjdir','var') | isempty(subjdir)
  subjdir = getenv('SUBJECTS_DIR');
  if isempty(subjdir)
    fprintf('%s: ERROR: SUBJECTS_DIR not defined as an environment variable\n',mfilename);
    return;
  end;
end;

if ~exist('matfile','var') | isempty(matfile)
	fprintf('%s: ERROR: you need to define the matfile\n',mfilename);
end
	


if ~exist('min_lat','var') | isempty(min_lat)
	fprintf('%s: ERROR: you need to define the min_lat in order to set the time range for max and avg calculations\n',mfilename); 
end


if ~exist('max_lat','var') | isempty(max_lat)
	fprintf('%s: ERROR: you need to define the max_lat in order to set the time range for max and avg calculations\n',mfilename); 
end
 

if ~exist('mbmask_flag','var') | isempty(mbmask_flag), mbmask_flag=0; end;

fprintf('%s: reading stcfile %s...\n',mfilename,stcfile);
[starttime,sample_period,vertices,sol]=ts_read_stc(stcfile);

fprintf('%s: start time = %0.2f ms\n',mfilename,starttime);
fprintf('%s: sample_period = %0.2f ms\n',mfilename,sample_period);



fprintf('%s: loading the file that includes time information: %s ...\n',mfilename,matfile);
load (strvcat(matfile));


fprintf('%s: finding the sample num corresponding to max and min latencies ...\n',mfilename);
time_vec = avg_data.averages.time;
time_tic = round(time_vec*1000);
time_vec_min = match(min_lat,time_tic,1,2);
time_vec_max = match(max_lat,time_tic,1,2);
fprintf('%s: for latency %d found samplepoint %d ...\n',mfilename,min_lat,time_vec_min);
fprintf('%s: for latency %d found samplepoint %d ...\n',mfilename,max_lat,time_vec_max);


fprintf('%s: extracting the desired range of STC from time %d to %d ...\n',mfilename,min_lat,max_lat);
cut_sol = sol(:,time_vec_min:time_vec_max);
[y,x] = size(cut_sol);


fprintf('%s: creating the average of STC from time %d to %d ...\n',mfilename,min_lat,max_lat);
avg_sol = sum(cut_sol,2)/x;

fprintf('%s: finding the max of STC from time %d to %d ...\n',mfilename,min_lat,max_lat);
max_sol = max(cut_sol,[],2);

fprintf('%s: replacing stcfile with max and average of stc...\n',mfilename);
sol = [];
sol = [max_sol avg_sol];
clear cut_sol;


nverts = length(vertices);
[ndips,tpoints] = size(sol);

if t0<1, t0 = 1; end;
if t1<1 | t1>tpoints, t1 = tpoints; end;
if t0>t1, t0=t1; end;

if (nverts ~= ndips)
  fprintf('%s: error: num vertices (%d) does not match num dips (%d)!\n',...
    mfilename,nverts,ndips);
  return;
end;

fprintf('%s: num vertices = %d\n',mfilename,nverts);
fprintf('%s: num time points = %d\n',mfilename,tpoints);

if sparsesmooth | postsmooth
  fprintf('%s: reading %s surface for subject %s...\n',...
    mfilename,hemi,subj);
  surf = fs_load_subj(subj,hemi);
  surf = fs_find_neighbors(surf);
else
  fprintf('%s: reading number of %s verts for subject %s...\n',...
    mfilename,hemi,subj);
  surf = fs_load_subj(subj,hemi,[],1,subjdir);
end;

tpoints = t1-t0+1;
surfstats = zeros(surf.nverts,tpoints);
surfstats(vertices+1,:) = sol(:,t0:t1);

% mask midbrain
if mbmask_flag==1
  surfstats = fs_mask_surfstats_with_aparc(surfstats,subj,hemi,subjdir);
end;

if sparsesmooth>0 | postsmooth>0
  fprintf('%s: smoothing timepoints %d-%d of %s...\n',mfilename,t0,t1,stcfile);
  for t=1:tpoints
    vals = surfstats(:,t);
    if sparsesmooth>0
      vals = fs_smooth_sparse(surf,vals,sparsesmooth);
    end;
    if postsmooth>0
      vals = fs_smooth(surf,vals,postsmooth);
    end;
    surfstats(:,t) = vals;
  end;
end

fprintf('%s: writing mghfile %s...\n',mfilename,mghfile);
vol = reshape(surfstats,[surf.nverts,1,1,tpoints]);
status = fs_save_mgh(vol,mghfile);

if status
  fprintf('%s: error writing mghfile\n',mfilename);
else
  fprintf('%s: finished writing mghfile\n',mfilename);
end;
