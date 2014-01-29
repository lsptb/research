function results = fs_multi_roi(fname_data,fname_roi_list,minval,weighted_avg_flag);
%function results = fs_multi_roi(fname_data,fname_roi,[minval],[weighted_avg_flag]);
%
% Required Input:
%  fname_data: full path to data volume (mgh/mgz format)
%  fname_roi: full path to ROI volume (mgh/mgz format)
%    -must be registered to and have same resolution as data volume
%    -can be a cell array (curly bracketed list) of multiple ROI file names
%
% Optional Input:
%   minval: minimum value
%     {default: 10^-5}
%   weighted_avg_flag: [0|1] toggle calculation of weighted average
%     if 1, use values in ROI volumes to weight the data volume
%     {default: 0}
%
% Created:  04/05/07 by Don Hagler
% Last mod: 07/14/08 by Don Hagler
% Last mod: 11/06/08 by Don Hagler
%

results = [];
if (~mmil_check_nargs(nargin, 2)), return; end;

if ~exist('minval','var') | isempty(minval), minval = 10^-5; end;
if ~exist('weighted_avg_flag','var') | isempty(weighted_avg_flag)
  weighted_avg_flag = 0;
end;

if ~iscell(fname_roi_list), fname_roi_list = {fname_roi_list}; end;

if ~exist(fname_data,'file')
  error('data file %s not found',fname_data);
end;
[vol_data,M,mr_parms,data_volsz] = fs_load_mgh(fname_data,[],[],1);
for i=1:length(fname_roi_list)
  fname_roi = fname_roi_list{i};
  if ~isempty(fname_roi)
    if ~exist(fname_roi,'file')
      error('ROI file %s not found',fname_roi);
    end;
    [vol_roi,M,mr_parms,roi_volsz] = fs_load_mgh(fname_roi,[],[],1);
    if any(data_volsz(1:3)~=roi_volsz(1:3))
      error('ROI vol size for file %s does not match data',fname_roi);
    end;
    if roi_volsz(4)>1
      error('ROI volumes cannot be multi-frame (%s has %d frames)',...
        fname_roi,roi_volsz(4));
    end;
  end;
end;

% load data file
fprintf('%s: loading data...\n',mfilename);
[vol_data,M,mr_parms,data_volsz] = fs_load_mgh(fname_data);
vol_data = reshape(vol_data,[prod(data_volsz(1:3)),data_volsz(4)]); % allow multiframe

fprintf('%s: extracting values...\n',mfilename);
for i=1:length(fname_roi_list)
  fname_roi = fname_roi_list{i};
  if isempty(fname_roi)
    roi = [];
  else
    [vol_roi,M,mr_parms,roi_volsz] = fs_load_mgh(fname_roi);
    roi = find(vol_roi>0);
  end;
  if isempty(roi)
    results(i).raw_vals = [];
    results(i).vals = [];
    results(i).nvox = 0;
    results(i).nvals = 0;
    results(i).avg = NaN;
    results(i).stdv = NaN;
  else
    raw_vals = vol_data(roi,:);
    results(i).raw_vals = raw_vals;
    nframes = size(raw_vals,2); % allow multiframe
    ind_good_vals = find(abs(raw_vals(:,1))>minval & ~isnan(raw_vals(:,1)));
    nvals = length(ind_good_vals);
    results(i).vals = zeros(nvals,nframes);
    results(i).avg = zeros(1,nframes);
    results(i).stdv = zeros(1,nframes);
    if weighted_avg_flag
      weights = vol_roi(roi);
      results(i).weights = weights;
      weights = weights*ones(1,nframes);
      raw_vals_sqrd = raw_vals.^2;
      wtd_vals_sqrd = raw_vals_sqrd .* weights;
      wtd_vals = raw_vals .* weights;
      for f=1:nframes
        tmp_vals = raw_vals(ind_good_vals,f);
        tmp_wtd_vals = wtd_vals(ind_good_vals,f);
        tmp_wtd_vals_sqrd = wtd_vals_sqrd(ind_good_vals,f);
        tmp_weights = weights(ind_good_vals,f);

        wtd_sum = sum(tmp_wtd_vals);
        wtd_sumsq = sum(tmp_wtd_vals_sqrd);
        sum_weights = sum(tmp_weights);
        sumsq_weights = sum(tmp_weights.^2);
        mean_weights = mean(tmp_weights);
        wtd_mean = wtd_sum / sum_weights;
        nvals = length(tmp_wtd_vals);
        results(i).vals(:,f) = tmp_wtd_vals;
        results(i).avg(f) = wtd_mean;
        % from: Bland, JM and Kerry, SM. "Weighted comparison of means". BMJ. 1998. 316:129
  %      wtd_stdv = sqrt(((wtd_sumsq / mean_weights) - nvals*(wtd_mean.^2))/(nvals-1));
        % from: http://pygsl.sourceforge.net/reference/pygsl/node36.html
        wtd_stdv = sqrt((sum_weights/(sum_weights.^2 - sumsq_weights))*sum(tmp_weights.*(tmp_vals-wtd_mean).^2));
        results(i).stdv(f) = wtd_stdv;
      end;
    else
      tmp_vals = raw_vals(ind_good_vals,:);
      results(i).vals = tmp_vals;
      results(i).avg = mean(tmp_vals,1);
      results(i).stdv = std(tmp_vals,0,1);
    end;
    results(i).nvox = size(raw_vals,1);
    results(i).nvals = nvals;
  end;
end

