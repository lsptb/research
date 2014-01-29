function data = calc_analytic_amplitude(data)

[datatype,datafield,dataparam] = ts_object_info(data);

ncond = length(data.(datafield));
nchan = data.num_sensors;
ndim  = ndims(data.(datafield)(1).(dataparam{1}));

for c = 1:ncond
  if ndim == 2
    data.(datafield)(c).(dataparam{1}) = abs(hilbert(data.(datafield)(c).(dataparam{1})'))';
  elseif ndim == 3
    for ch = 1:nchan
      dat = squeeze(data.(datafield)(c).(dataparam{1})(ch,:,:)); % time x trial
      dat = abs(hilbert(dat));
      data.(datafield)(c).(dataparam{1})(ch,:,:) = dat;
    end
  end
end