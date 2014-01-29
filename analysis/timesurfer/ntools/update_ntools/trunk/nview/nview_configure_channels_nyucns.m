%
%
% configures channels to map 1-128 -> 1-128 and duplicate to 129 -> 256 for
% nyucns nspike hardware
i=1;
for twice = 0:1
    for group = 0:3
        for division = 0:7
            for channel = 0:3
                hw_num = (division * 16) + (group * 4) + channel;
                experiment.channels(i).hardware_number = hw_num; 
                i = i + 1;
             end
        end
    end
end 


for i = 1:128
    experiment.channels(i).name = num2str(i);
    experiment.channels(i).lowpass_cutoff = 11000;
    experiment.channels(i).highpass_cutoff = 1;
end
for i = 129:256
    experiment.channels(i).name = num2str(i);
    experiment.channels(i).lowpass_cutoff = 11000;
    experiment.channels(i).highpass_cutoff = 300;
end

