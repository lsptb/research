function lfp = lpfilter(data,lowpass_freq,orig_sampling_rate,new_sampling_rate);
%   
%   lfp = lpfilter(data,lowpass_freq,orig_sampling_rate,new_sampling_rate);
%

lfp = (mtfilter(data,[0.005,lowpass_freq],orig_sampling_rate,0));
decimate_factor = round(orig_sampling_rate./new_sampling_rate);
lfp = lfp(:,1:decimate_factor:end);
