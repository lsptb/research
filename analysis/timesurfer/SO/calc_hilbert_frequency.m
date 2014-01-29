function data = calc_hilbert_frequency(data)

[datatype,datafield,dataparam] = ts_object_info(data);

ncond = length(data.(datafield));
nchan = data.num_sensors;
ndim  = ndims(data.(datafield)(1).(dataparam{1}));
t     = data.(datafield)(1).time;

for c = 1:ncond
  if ndim == 2
    dat = data.(datafield)(c).(dataparam{1})'; % time x chan
    dat = angle(hilbert(dat));
    dat = (1/(2*pi)) * diff(unwrap(dat))./diff(t);
    data.(datafield)(c).(dataparam{1}) = dat';
  elseif ndim == 3
    for ch = 1:nchan
      dat = squeeze(data.(datafield)(c).(dataparam{1})(ch,:,:)); % time x trial
      dat = angle(hilbert(dat));
      dat = (1/(2*pi)) * diff(unwrap(dat))./diff(t);
      data.(datafield)(c).(dataparam{1})(ch,:,:) = dat;
    end
  end
end