function data = calc_hilbert_phase(data)

[datatype,datafield,dataparam] = ts_object_info(data);

ncond = length(data.(datafield));
nchan = data.num_sensors;
ndim  = ndims(data.(datafield)(1).(dataparam{1}));

for c = 1:ncond
  if ndim == 2
    data.(datafield)(c).(dataparam{1}) = angle(hilbert(data.(datafield)(c).(dataparam{1})'))';
  elseif ndim == 3
    for ch = 1:nchan
      dat = squeeze(data.(datafield)(c).(dataparam{1})(ch,:,:)); % time x trial
      dat = angle(hilbert(dat));
      data.(datafield)(c).(dataparam{1})(ch,:,:) = dat;
    end
  end
end