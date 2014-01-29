% [channels] = ntools_get_channels(data_type)
%
% For a given file's data type, computes the file's number of channelS

function [channels] = ntools_get_channels(data_type)
    warning('Deprecated! Will probably be removed soon!');

    switch(data_type)
        case {'nspike','decieeg','ieeg', 'cleanieeg','cleandecieeg'},
            channels = 256;
        case {'audio','videosync'}
            channels = 1;
        otherwise
            error('Unknown data_type: "%s"', data_type);
    end
end
