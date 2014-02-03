res = load('/home/jsherfey/svn/dev/onestream/test_20110125/sub16_recon/bem/lh_white.dip');


mghlh1 = '/home/jsherfey/svn/dev/onestream/test_20110125/mghfiles/test_dSPM_WordNPNW_grad_bem_cond01-spsm10-sm10-ico7-sm3-lh.mgh';
mghlh2 = '/home/jsherfey/svn/dev/onestream/test_20110125/mghfiles/test_dSPM_WordNPNW_grad_bem_cond01-spsm10-sm10-lh.mgh';
mghrh1 = '/home/jsherfey/svn/dev/onestream/test_20110125/mghfiles/test_dSPM_WordNPNW_grad_bem_cond01-spsm10-sm10-ico7-sm3-rh.mgh';
mghrh2 = '/home/jsherfey/svn/dev/onestream/test_20110125/mghfiles/test_dSPM_WordNPNW_grad_bem_cond01-spsm10-sm10-rh.mgh';

mghfile = mghlh1;
[vol, M, mr_parms, volsz] = fs_load_mgh(mghfile);


%% group averaging

% 1) code copied from ts_groupavg_old.m

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


% 2) code copied from ts_groupavg.m

% load precalculated group means and standard deviations
volA_mean = fs_load_mgh(innameA_mean);
volA_stdv = fs_load_mgh(innameA_stdv);
volB_mean = fs_load_mgh(innameB_mean);
volB_stdv = fs_load_mgh(innameB_stdv);
% load num subjects
clear n; load(innameA_n); nA = n;
clear n; load(innameB_n); nB = n;
n = nA + nB;
df = n-2;
% calculate difference of means
outvol_mean = volA_mean - volB_mean;
% calculate pooled standard deviations
outvol_stdv = sqrt(((nA-1)*volA_stdv.^2+...
                    (nB-1)*volB_stdv.^2)/(n-2));
% calculate tstats
outvol_tstat = outvol_mean ./ (outvol_stdv ./ sqrt(eps+n) + eps);
% calculate pvals
outvol_pval = -log10(2*tcdf(-abs(outvol_tstat),df)); % 2-tailed



