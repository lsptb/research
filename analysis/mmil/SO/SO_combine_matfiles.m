function DATA = SO_combine_matfiles(matfiles,params,fid,nogap,chantype)
if nargin < 4, nogap = 0; end
if nargin < 5, chantype = []; end
if ~iscell(matfiles), matfiles = {matfiles}; end
for f = 1:length(matfiles)
  load(matfiles{f}); % data
  if ~isempty(chantype)
    data = ts_data_selection(data,'chantype',chantype);
  end
  % if preproc params are provided, preprocess the data
  if nargin > 1 && isstruct(params)
    try
      if nargin < 3, fid = 1; end
      fprintf(fid,'Preprocessing data from single MAT-file %g of %g\n',f,length(matfiles));
      data = ts_preproc(data, 'dsfact',     params.ICA_dsfact,      ...
                              'bpfilter',   'yes',  'bpfreq',   params.ICA_bpfreq, 'bandpass_detrend_flag',0,...
                              'notch_flag', params.ICA_notchflag,      ...
                              'blc',        params.ICA_blcflag,  'blcwindow',params.ICA_blcwindow  );  
    end
  end
  if nogap
    data.epochs.time = data.epochs.time - min(data.epochs.time);
  end
  if f == 1
    DATA = data;
  else
    t0 = DATA.epochs.time;
    tn = data.epochs.time;
    DATA.epochs.time = [t0 tn+t0(end)+(1/data.sfreq)];
    DATA.epochs.data = cat(2,DATA.epochs.data,data.epochs.data);
  end
  DATA.epochs.matfiles{f} = matfiles{f};
  DATA.epochs.duration(f) = data.epochs.time(end);
  clear data t0 tn
end
