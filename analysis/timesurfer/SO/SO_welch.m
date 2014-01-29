function outdata = SO_welch(data,frqs,dsfact)
if nargin < 3, dsfact = 1;   end
if nargin < 2, frqs   = 128; end

if dsfact > 1
  data.sfreq = data.sfreq / dsfact;
end
Fs = data.sfreq;

[datatype datafield]    = ts_object_info(data);
outdata                 = rmfield(data,datafield);
outdata.averages        = data.(datafield);
[outdata.averages.data] = deal([]);

for c = 1:length(data.(datafield))
  dat = data.(datafield)(c).data;
  if dsfact > 1
    dat = decimate(dat,dsfact);
  end
  clear P
  for ch = 1:size(dat,1)
    x = dat(ch,:);
    [p f] = pwelch(x,[],[],frqs,Fs);
    P(ch,:) = log10(p);
  end
  outdata.averages(c).time  = f;
  outdata.averages(c).data  = P;
  outdata.averages(c).stdev = nan(size(P));
end
outdata.sfreq = 1 / (f(2)-f(1));
