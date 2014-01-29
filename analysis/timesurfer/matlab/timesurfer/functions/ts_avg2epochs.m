function outdata = ts_avg2epochs(indata)

[datatype,datafield,dataparam] = ts_object_info(indata);
if strcmp(datatype,'avg_data')
  outfield = 'epochs';
elseif strcmp(datatype,'timefreq_data')
  outfield = 'timefreq';
else
  warning('data structure not recognized');
  return;
end

outdata = indata;
outdata.(outfield) = indata.(datafield);
outdata            = rmfield(outdata,datafield);

nc = length(outdata.(outfield));
[outdata.(outfield)(1:nc).num_trials] = deal(1);

if isfield(outdata.(outfield),'stdev')
  outdata.(outfield) = rmfield(outdata.(outfield),'stdev');
end