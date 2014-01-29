function results = fs_groupavg_glm(fnamelist,varargin)
%function results = fs_groupavg_glm(fnamelist,varargin)
%
%  Purpose: calculate group averages and t-stats from surface stats
%   in mgh format, resampled to common space of FreeSurfer's icosahedral sphere
%
% Usage:
%  fs_groupavg_glm(fnamelist, 'key1', value1,...);
%
% Required Input:
%   fnamelist: cell array of input full path file names
%     all input files should be mgh format surface stats files
%     registered to common (sphere) space
%
%     fnamelist can also be a cell matrix of file names
%       Each row should correspond to a single subject
%       Each column should correspond to a condition
%
%  Optional Input:
%   'condnames' - cell array of condition names
%     Length should match number of columns in fnamelist
%     If ommitted, will assign names like '1', '2', etc.
%   'regressors' - matrix containing regressors for each entry in fnamelist
%     Each column should contain a different numerical regressor
%     Each row should contain values for a different subject
%     If ommitted, will simply calculate average and t-stats
%   'regnames' - cell array of regressor names
%     Length should match number of columns in regressors
%     If ommitted, will assign names like reg1, reg2, reg3, etc.
%   'contrast_vectors': cell array of contrast vectors
%     e.g. {[0 0 1 0],[0 0 0 1],[0 1 -1 0]}
%     Number of entries in contrast vectors should equal:
%        1 + (number of columns in fnamelist) + (number of columns in regressors)
%     First entry corresponds to baseline (intercept)
%     Next entries correspond to conditions
%     Last entries correspond to regressors
%     Results (means, t-stats, p-vals, etc.) will be generated for
%       each contrast vector
%     If ommitted, will generate a contrast vector for each condition and regressor
%   'contrast_names' - cell array of contrast names
%     Must have one for each contrast_vector in contrast_vectors
%     If contrast_vectors and contrast_names are both ommitted,
%       will use condnames and regnames as contrast names
%     If contrast_vectors is supplied but contrast_names is not,
%       will use cont1, cont2, cont3, etc. for contrast names
%   'offset' - value subtracted from all data points before calculating means and t-stats
%     {default = 0}
%   'frame' - frame number (e.g. for files with multiple time points)
%     {default = 1}
%   'data_matfile' - name of matlab mat file to store data in fnamelist
%     If file exists, will load data from it instead of files in fnamelist
%       (but fnamelist is still required)
%     If file does not exist, will save data matrix into it
%   'verbose' - [0|1] whether to print status information to the screen
%     {default: 1}
%
% Output:
%   results: structure containing group averages, standard deviations
%     t-stats, and significance maps (-log10(p))
%
%   Separate maps will be generated for each contrast vector
%
% Created:  08/16/07 by Don Hagler
% Last Mod: 04/14/09 by Don Hagler

results = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parse arguments
if (~mmil_check_nargs(nargin, 1)), return; end;
parms = mmil_args2parms( varargin, { ...
  'condnames',[],[],...
  'regressors',[],[],...
  'regnames',[],[],...
  'contrast_vectors',[],[],...
  'contrast_names',[],[],...
  'offset',0,[],...
  'frame',1,[1 Inf],...
  'data_matfile',[],[],...
  'verbose',true,sort([false true]),...
});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check parms
if isempty(fnamelist)
  error('input file name list is empty');
elseif ~iscell(fnamelist)
  fnamelist = {fnamelist};
end;
[nsubs,nconds] = size(fnamelist);
if nsubs==1
  error('input file name list has data for only one subject');
end;

if ~isempty(parms.condnames)
  if ~iscell(parms.condnames), parms.condnames = {parms.condnames}; end;
  if length(parms.condnames)~=nconds
    error('number of condnames (%d) does not match number of conditions (%d)',...
      length(parms.condnames),nconds);
  end;
else
  for i=1:nconds
    parms.condnames{i} = sprintf('cond%d',i);
  end;
end;

if ~isempty(parms.regressors)
  [tmp_nsubs,nregs] = size(parms.regressors);
  if tmp_nsubs~=nsubs
    error('number of rows in regressors (%d) does not match fnamelist (%d)',...
      tmp_nsubs,nsubs);
  end;
  if ~isempty(parms.regnames)
    if ~iscell(parms.regnames), parms.regnames = {parms.regnames}; end;
    if length(parms.regnames)~=nregs
      error('number of regnames (%d) does not match number of regressors (%d)',...
        length(parms.regnames),nregs);
    end;
  else
    for i=1:nregs
      parms.regnames{i} = sprintf('reg%d',i);
    end;
  end;
else
  nregs = 0;
end;

if parms.verbose
  fprintf('%s: %d subjects, %d conditions, %d regressors\n',...
    mfilename,nsubs,nconds,nregs);
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create design matrix
%% WARNING: this is a fixed effects, non-repeated measures design
if nconds > 1
  ncols = 1 + nconds + nregs;
else
  ncols = 1 + nregs;
end
nrows = nsubs*nconds;

X = zeros(nrows,ncols);
X(:,1) = ones(nrows,1); % baseline
j = 1;
if nconds > 1
  for i=1:nconds
    for s=1:nsubs
      X(j,1+i) = 1;
      j = j + 1;
    end;
  end;
  for i=1:nregs
    X(:,1+nconds+i) = parms.regressors(:,i);
  end;
else
  for i=1:nregs
    X(:,1+i) = parms.regressors(:,i);
  end;
end

% contrast vectors
if ~isempty(parms.contrast_vectors)
  % make sure it is cell array
  if ~iscell(parms.contrast_vectors), parms.contrast_vectors = {parms.contrast_vectors}; end;
  % check that length of parms.contrast_names matches
  if ~isempty(parms.contrast_names)
    if ~iscell(parms.contrast_names), parms.contrast_names = {parms.contrast_names}; end;
    if length(parms.contrast_vectors)~=length(parms.contrast_names)
      error('length of contrast_names must match contrast_vectors');
    end;
  end;
  % check size of contrast vectors
  for i=1:length(parms.contrast_vectors)
    cvec = parms.contrast_vectors{i};
    if any(size(cvec)~=[1 ncols])
      error('contrast vector %d has wrong size (should be 1 x %d)',i,ncols);
    end;
  end;
else
  % generate contrast_vectors (one for each condition and regressor)
  k = 1;
  parms.contrast_names = [];
  for i=1:nconds
    tmp = zeros(1,ncols);
    tmp(k) = 1;
    parms.contrast_vectors{k} = tmp;
    parms.contrast_names{k} = parms.condnames{i};
    k = k + 1;
  end;
  for i=1:nregs
    tmp = zeros(1,ncols);
    tmp(k) = 1;
    parms.contrast_vectors{k} = tmp;
    parms.contrast_names{k} = parms.regnames{i};
    k = k + 1;
  end;
end;

% check files
if parms.verbose
  fprintf('%s: checking input files...\n',mfilename);
end;
nverts = -1;
rejectlist = zeros(nsubs,1);
for s=1:nsubs
  for c=1:nconds
    if isempty(fnamelist{s,c})
      if parms.verbose
        fprintf('%s: WARNING: rejecting subject %d because of missing file\n',...
          mfilename,s);
      end;
      rejectlist(s) = 1;
      continue;
    end;
    [vol,M,mrparms,volsz] = fs_load_mgh(fnamelist{s,c},[],[],1); % header-only
    if parms.frame>volsz(4)
      error('only %d frames in file %s',volsz(4),fnamelist{s,c});
    end;
    if volsz(2)>1 | volsz(3)>1
      error('this function only accepts surface (not volume) data (check %s)',...
        fnamelist{s,c});
    end;
    tmp_nverts = volsz(1);
    if nverts==-1
      nverts = tmp_nverts;
    elseif nverts~=tmp_nverts
      error('number of vertices in %s (%d) does not match that in %s (%d)\n',...
        fnamelist{s,c},tmp_nverts,fnamelist{1,1},nverts);
    end;
  end;
end;

% exclude bad subjects
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
j = 1;
ind_keep_rows = [];
for s=1:nsubs
  for c=1:nconds
    if ~rejectlist(s)
      ind_keep_rows = [ind_keep_rows,j];
    end;
    j=j+1;
  end;
end;
if parms.verbose
  fprintf('%s: excluding %d subjects (empty entries in fnamelist)\n',...
    mfilename,nsubs-nsubs_keep);
end;
fnamelist = tmp_fnamelist;
X = X(ind_keep_rows,:);
nsubs = nsubs_keep;
nrows = length(ind_keep_rows);

% load data
tic;
datamat = [];
if ~isempty(parms.data_matfile) & exist(parms.data_matfile,'file')
  fprintf('%s: loading data from matfile %s...\n',mfilename,parms.data_matfile);
  load(parms.data_matfile);
else
  fprintf('%s: loading data from %d files in fnamelist...\n',mfilename,nsubs);
  datamat = zeros(nrows,nverts);
  j = 1;
  for s=1:nsubs
    for c=1:nconds
      vec = fs_load_mgh(fnamelist{s,c},[],parms.frame);
      vec = reshape(vec,[prod(size(vec)) 1]);
      datamat(j,:) = vec - parms.offset;
      j=j+1;
    end;
  end;
  if ~isempty(parms.data_matfile)
    fprintf('%s: saving data to matfile %s...\n',mfilename,parms.data_matfile);
    save(parms.data_matfile,'datamat');
  end;
end;
fprintf('%s: finished loading data after %0.1f seconds\n',mfilename,toc);

% calc glm
fprintf('%s: running GLM calculations...\n',mfilename);
tic;
results = calc_voxelwise_glm(X,datamat,parms.contrast_vectors,parms.contrast_names);
fprintf('%s: finished GLM calculations after %0.1f seconds\n',mfilename,toc);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function results = calc_voxelwise_glm(X,datamat,cont_vecs,cont_names)
  results = [];

  results.nvols = size(X,1);
  results.nparams = size(X,2);
  % datamat must be nvols x nvox
  if size(datamat,1) ~= results.nvols
    error('size of datamat does not match X');
  end;
  results.nvox = size(datamat,2);

  Xi = pinv(X);
  C = Xi*Xi';

  results.betas = Xi*datamat;
  datamat_hat = X*results.betas;
  sserr = sum((datamat_hat-datamat).^2,1);
  results.sigma = sqrt(sserr/(results.nvols-results.nparams));
  results.dof = results.nvols - results.nparams; % should this vary by contrast?

  results.contrasts = [];
  for i=1:length(cont_vecs)
    contrast_vec = cont_vecs{i};
    results.contrasts(i).name = cont_names{i};
    sigma_contrast = sqrt(contrast_vec*C*contrast_vec');

    results.contrasts(i).mean = contrast_vec*results.betas;
    results.contrasts(i).stdv = sigma_contrast*results.sigma;
    results.contrasts(i).mean(find(abs(results.contrasts(i).mean)<10*eps)) = 0;
    results.contrasts(i).tstats = results.contrasts(i).mean./(sigma_contrast*results.sigma+eps);
    results.contrasts(i).pvals = (results.contrasts(i).tstats<0).*...
                                 log10(2*tcdf(results.contrasts(i).tstats,results.dof)) +...
                                 (results.contrasts(i).tstats>0).*...
                                 (-log10(2*tcdf(-results.contrasts(i).tstats,results.dof)));
  end;

return;
