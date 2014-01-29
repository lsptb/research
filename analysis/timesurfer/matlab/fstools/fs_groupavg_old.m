function results = fs_groupavg(fnamelist,gnamelist,condition_weights,frame,offset,verbose,contrasts_flag)
%function results = fs_groupavg(fnamelist,[gnamelist],[condition_weights],[frame],...
%   [offset],[verbose],[contrasts_flag])
%
% Required Input:
%   fnamelist: cell array of input full path file names
%     all input files should be mgh format surface or volume files
%     registered to common space (e.g. spherical surface)
%     
%     fnamelist can also be a cell matrix of file names
%       Each row should correspond to a single subject
%       Each column should correspond to a condition
%     Specify a linear combination of these conditions with
%       condition_weights
%     For a paired t-test between two conditions, supply
%       fnamelist with two columns and condition_weights = [1 -1]
%
%  Optional Input:
%   gnamelist: cell array of group names
%     one for each file name
%     If ommitted, will calculate average of all files
%       as single group
%   condition_weights: vector of linear weights for each
%     condition (each column of fnamelist)
%     {default = [1 1 1 ...]}
%   frame: frame number (e.g. for files with multiple time points)
%     {default = 1}
%   offset: value subtracted from all data points before calculating means and t-stats
%     {default = 0}
%   verbose: [0|1] whether to print status messages
%     {default = 1}
%   contrasts_flag: [0|1] whether to calculate contrasts between groups
%     {default = 1}
%
% Output:
%   results: structure containing group averages, standard deviations
%     t-stats, and significance maps (-log10(p))
%
%   Separate maps will be generated for each group
%     and for each pairwise comparison between groups (if contrasts_flag=1).
%
% Created:  08/16/07 by Don Hagler
% Last Mod: 07/23/08 by Don Hagler
%

%% todo: allow input of vector of scaling values for each entry in fnamelist? or negflags?
%% todo: use args2parms

results = [];

if (~mmil_check_nargs(nargin,1)) return; end;

if ~exist('offset','var') || isempty(offset), offset = 0; end;
if ~exist('verbose','var') || isempty(verbose), verbose = 1; end;
if ~exist('frame','var') || isempty(frame), frame = 1; end;
if ~exist('contrasts_flag','var') || isempty(contrasts_flag), contrasts_flag = 1; end;

if isempty(fnamelist)
  error('input file name list is empty');
elseif ~iscell(fnamelist)
  fnamelist = {fnamelist};
end;
[nsubs,nconds] = size(fnamelist);
if nsubs==1 && nconds==1
  error('input file name list contains only one file');
elseif nsubs==1
  fnamelist = reshape(fnamelist,[nconds,nsubs]);
  [nsubs,nconds] = size(fnamelist);
elseif nsubs<nconds
  if verbose
    fprintf('%s: WARNING: fnamelist has more columns (conditions) than rows (subjects)\n',mfilename);
  end;
end;

if ~exist('gnamelist','var') | isempty(gnamelist)
  gnamelist = {};
  for s=1:nsubs
    gnamelist{s} = 'group';
  end;
else
  if ~iscell(gnamelist)
    gnamelist = {gnamelist};
  end;
  if length(gnamelist)==1
    for s=1:nsubs
      gnamelist{s} = gnamelist{1};
    end;
  end;
  if length(gnamelist)~=nsubs
    error('length of gnamelist (%d) does not match number of rows (subjects) of fnamelist (%d)',...
      length(gnamelist),nsubs);
  end;
end;

if length(frame)>1
  error('frame must be a single number, not a vector');
end;

if ~exist('condition_weights','var') | isempty(condition_weights)
  condition_weights = ones(nconds,1);
end;
condition_weights = reshape(condition_weights,[prod(size(condition_weights)),1]);
if size(condition_weights,1) ~= nconds
  error('number of condition_weights (%d) does not match number of columns (conditions) of fnamelist (%d)',...
    size(condition_weights,1),nconds);
end;


if verbose
  fprintf('%s: checking input files...\n',mfilename);
end;
% check files
orig_volsz = [];
rejectlist = zeros(nsubs,1);
for s=1:nsubs
  for c=1:nconds
    if isempty(fnamelist{s,c})
      rejectlist(s) = 1;
      continue;
    end;
    [vol,M,mrparms,volsz] = fs_load_mgh(fnamelist{s,c},[],[],1); % header-only
    if length(volsz)==3, volsz = [volsz 1]; end;
    if frame>volsz(4)
      error('only %d frames in file %s',volsz(4),fnamelist{s,c});
    end;
    if isempty(orig_volsz)
      orig_volsz = volsz;
    end;
    if any(orig_volsz~=volsz)
      error('size of input volume for\n%s (%d,%d,%d,%d)\ndoes not match that of\n%s (%d,%d,%d,%d)',...
        fnamelist{s,c},volsz(1),volsz(2),volsz(3),volsz(4),...
        fnamelist{1,1},orig_volsz(1),orig_volsz(2),orig_volsz(3),orig_volsz(4));
    end;
  end;
end;
nvals = prod(orig_volsz(1:3));
gnamelist = {gnamelist{~rejectlist}};

results.volsz = orig_volsz(1:3);
results.nvals = nvals;

nsubs_keep = length(find(~rejectlist));
tmp_fnamelist = cell(nsubs_keep,nconds);
k=1;
for s=1:nsubs
  if ~rejectlist(s)
    for c=1:nconds
      tmp_fnamelist{k,c} = fnamelist{s,c};
    end;
    k=k+1;
  end;
end;
fnamelist = tmp_fnamelist;
nsubs = nsubs_keep;

if contrasts_flag
  % prepare contrasts
  groupNames = unique(gnamelist);
  numGroups = length(groupNames);
  gnumlist = zeros(length(gnamelist));
  c=1;
  results.groups = [];
  results.contrasts = [];
  for i=1:length(groupNames)
    gnumlist(find(strcmp(groupNames{i},gnamelist)))=i;
    results.groups(i).name = groupNames{i};
    results.groups(i).n = 0;
    results.groups(i).mean = zeros(nvals,1);
    results.groups(i).stdv = zeros(nvals,1);
    results.groups(i).df = [];
    results.groups(i).tstat = [];
    results.groups(i).sig = [];
    for j=1:length(groupNames)
      if i==j, continue; end;
      results.contrasts(c).name = sprintf('%s-VS-%s',groupNames{i},groupNames{j});
      results.contrasts(c).n = [];
      results.contrasts(c).df = [];
      results.contrasts(c).mean = [];
      results.contrasts(c).stdv = [];
      results.contrasts(c).tstat = [];
      results.contrasts(c).sig = [];
      c = c + 1;
    end;
  end;
end;

if verbose
  fprintf('%s: summing data...\n',mfilename);
end;
% load data and calculate sums and sums of squares
for s=1:nsubs
  data = zeros(nvals,nconds);
  for c=1:nconds
    vec = fs_load_mgh(fnamelist{s,c},[],frame);
    vec = reshape(vec,[prod(size(vec)) 1]);
    data(:,c) = vec - offset;
  end;
  vec = data*condition_weights;
  gnum = gnumlist(s);
  results.groups(gnum).n = results.groups(gnum).n + 1;
  results.groups(gnum).mean = results.groups(gnum).mean + vec;
  results.groups(gnum).stdv = results.groups(gnum).stdv + vec.^2;
end;

if verbose
  fprintf('%s: calculating stats...\n',mfilename);
end;
% calculate mean, stdv, tstats, sig for single conditions
for i=1:length(groupNames)
  n = results.groups(i).n;
  if n==0 % should not happen
    error('group %s has 0 subjects',groupNames{i});
  elseif n==1
    if verbose
      fprintf('%s: WARNING: group %s has 1 subject\n',mfilename,groupNames{i});
    end;
    results.stdv = zeros(nvals,1);
    continue;
  end;
  results.groups(i).stdv = sqrt((n*results.groups(i).stdv - results.groups(i).mean.^2)./(eps+n*(n-1)));
  results.groups(i).mean = results.groups(i).mean / (eps+n);
  results.groups(i).tstat = results.groups(i).mean ./ (results.groups(i).stdv ./ sqrt(eps+n));
  results.groups(i).tstat(results.groups(i).stdv==0)=0;
  results.groups(i).df = n-1;
  results.groups(i).sig = -log10(2*tcdf(-abs(results.groups(i).tstat),results.groups(i).df)); % 2-tailed
end;

if contrasts_flag
  % calculate mean diffs, pooled stdv, two-sample tstats, sig for contrasts
  c=1;
  for i=1:length(groupNames)
    for j=1:length(groupNames)
      if i==j, continue; end;
      n1 = results.groups(i).n;
      n2 = results.groups(j).n;
      n = n1+n2;
      results.contrasts(c).n = n;
      results.contrasts(c).df = results.groups(i).df + results.groups(j).df;
      results.contrasts(c).mean = results.groups(i).mean - results.groups(j).mean;
      results.contrasts(c).stdv = sqrt(((n1-1)*results.groups(i).stdv.^2+...
                                        (n2-1)*results.groups(j).stdv.^2)/(n-2));
      results.contrasts(c).tstat = results.contrasts(c).mean ./ (results.contrasts(c).stdv ./ sqrt(eps+n));
      results.contrasts(c).tstat(results.contrasts(c).stdv==0)=0;
      results.contrasts(i).sig = -log10(2*tcdf(-abs(results.contrasts(c).tstat),results.contrasts(c).df)); % 2-tailed
      c = c + 1;
    end;
  end;
end;
