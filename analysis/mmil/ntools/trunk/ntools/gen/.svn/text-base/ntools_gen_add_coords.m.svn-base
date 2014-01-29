% epoch_data = ntools_gen_add_coords(epoch_data, coords)
% avg_data = ntools_gen_add_coords(avg_data, coords)

function epoch_data = ntools_gen_add_coords(epoch_data, coords, wh_channels_to_keep)
    if(~exist('coords', 'var') || isempty(coords))
        [coord_filename, coord_pathname] = uigetfile('*', 'Select a file containing the coordinates, either a .txt file, or a .mat file with a coords variable'); if(coord_filename == 0), return; end
        [d f e] = fileparts(coord_filename);
        switch(lower(e))
            case '.mat'
                load(fullfile(coord_pathname, coord_filename));
            case '.txt'
                coords = load(fullfile(coord_pathname, coord_filename));
        end
    end
    
    if(exist('wh_channels_to_keep', 'var') && length(wh_channels_to_keep) ~= epoch_data.num_sensors)
        error('The number of channels to keep must match');
    elseif(size(coords, 1) ~= epoch_data.num_sensors)
        error('The number of coords must equal the number of electrodes');
    end
    
    if(~exist('wh_channels_to_keep', 'var') || isempty(wh_channels_to_keep))
        wh_channels_to_keep = 1:epoch_data.num_sensors;
    end
    
    for i=1:epoch_data.num_sensors
        epoch_data.sensor_info(i).coords = coords(wh_channels_to_keep(i),:);
    end
end