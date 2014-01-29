%
%
% this ugly for loop maps nspike channels 1-256 onto 1-256
function experiment = nview_configure_channels_256(experiment)
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
end