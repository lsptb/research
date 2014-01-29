%  ft_data = ntools_gen_ft_data(dat_file, dio_file, experiment, ['channels', wh_channels_to_keep], ...
%       ['max_epochs', max_epochs], ['is_continuous', 0|1])
%
%  Generates a fieldtrip data structure from an experiment data_structure
%  Data structure contains all experimental conditions
%  Trial event codes are stored in ft_data.cfg.trl(:,4)
%
%  Parameters:
%       dat_file: path to a .dat file
%       dio_file: path to a .dio file
%       experiment: experiment structure
%
%  Optional keyword parameters:
%       channels: which channels to keep (default: all)
%       max_epochs: max num epochs (default: no limit)
%       is_continuous: treats all data in file as one long trial (default:
%       no) - also see ntools_gen_cont_data 

function ft_data = ntools_gen_ft_data(dat_file, dio_file, experiment, varargin)
    [sfreq, num_channels, num_samples, precision, chan_names] = ntools_hdrxml_read(dat_file);

    max_epochs = Inf;
    wh_channels_to_keep = 1:num_channels;
    is_continuous = 0;
    
    for i = 1:length(varargin)
        if(ischar(varargin{i}))
            switch(lower(varargin{i}))
                case {'max_epochs', 'max epochs'}
                    max_epochs = varargin{i+1};
                case {'wh_channels_to_keep', 'channels'}
                    wh_channels_to_keep = varargin{i+1};
                case {'is_continuous', 'is_cont_data'}
                    is_continuous = varargin{i+1};
            end
        end
    end        

    ft_data = [];
    ft_data.cfg = [];
    ft_data.fsample = sfreq; % assuming this is in Hz
    ft_data = set_sensor_info(ft_data, chan_names, wh_channels_to_keep);
    [ft_data] = set_epochs(ft_data, dat_file, dio_file, experiment, sfreq, num_channels, num_samples, precision, max_epochs, wh_channels_to_keep, is_continuous);
end

function ft_data = set_sensor_info(ft_data, chan_names, wh_channels_to_keep)
    ft_data.label = chan_names;
    for i = wh_channels_to_keep
        ft_data.cfg.sensor_info(i).typestring = 'eeg';
        ft_data.cfg.sensor_info(i).type = 1;
        ft_data.cfg.sensor_info(i).kind = 2;
        ft_data.cfg.sensor_info(i).badchan = 0;
        ft_data.cfg.sensor_info(i).lognum = i;
        ft_data.cfg.sensor_info(i).loc = eye(4);
    end
end

function [ft_data] = set_epochs(ft_data, dat_file, dio_file, experiment, sfreq, num_channels, ...
        num_samples, precision, max_epochs, wh_channels_to_keep, is_continuous)
    NSPIKE_SAMPLING_RATE = 30000;  %  NSpike timestamps are in Nspike system sampling rate
    sfreq_to_nspike_sfreq_ratio = sfreq / NSPIKE_SAMPLING_RATE;
    samples_per_ms = sfreq / 1000;
    
    if(is_continuous)
        ft_data.cfg.name = 'Whole recording';
        %ft_data.cfg.num_rejects = struct('mag', 0, 'grad', 0, 'eeg', 0, 'eog', 0, 'manual', 0, 'skip', 0);
        ft_data.cfg.event = [];            
        ft_data.cfg.num_trials = 1;
        data = ntools_load_samples(dat_file, precision, num_channels);
        ft_data.cfg.trl(1,[1 3])  = 1;            
        ft_data.cfg.trl(1,2) = size(data,2); 
        ft_data.cfg.trl(1,4) = 0;            
        ft_data.time{1} = linspace(0, num_samples-1, num_samples) / sfreq;
        ft_data.trial{1} = double(data(wh_channels_to_keep, :));
    else
        [time_codes, A, B, C, D] = ntools_dio_parse(dio_file); %#ok;
        itrl = 0; 
        ft_data.cfg.event = experiment.event;            
        ft_data.cfg.num_conditions = length(experiment.event);            
   
        for i = 1:length(experiment.event) 
            experiment.event(i).timestamps = time_codes(C == experiment.event(i).source_code);
            
            ft_data.cfg.condition{i} = experiment.event(i).name;
            %ft_data.cfg.num_rejects(i) = struct('mag', 0, 'grad', 0, 'eeg', 0, 'eog', 0, 'manual', 0, 'skip', 0);

            % offsets in ms
            start_offset = experiment.event(i).start_offset;
            stop_offset = experiment.event(i).stop_offset;
            num_samples = round((stop_offset - start_offset) * samples_per_ms);

            num_trials = min(length(experiment.event(i).timestamps), max_epochs);
            ft_data.cfg.num_trials(i) = num_trials;
            
            if(num_trials == 0)
               error('No trials/epochs found in recording. Use ntools_gen_cont_data for continuous data.');
            end
            
            for j = 1:num_trials
                itrl = itrl+1;
                ft_data.time{itrl} = linspace(start_offset, stop_offset, num_samples) / 1000;
                timestamp = experiment.event(i).timestamps(j); % in 30 kHz
                start = round(timestamp * sfreq_to_nspike_sfreq_ratio + start_offset * samples_per_ms) + 1;
                stop = round(timestamp * sfreq_to_nspike_sfreq_ratio + stop_offset * samples_per_ms);
                ft_data.cfg.trl(itrl,1) = start; 
                ft_data.cfg.trl(itrl,2) = stop; 
                ft_data.cfg.trl(itrl,3) = start_offset; 
                ft_data.cfg.trl(itrl,4) = experiment.event(i).source_code;
                [data, eof] = ntools_load_samples(dat_file, precision, num_channels, start, stop);
                if(eof)
                    error('Trial %d with code %d extends beyond length of file', j, ft_data.epochs(i).event_code);
                end
                ft_data.trial{itrl} = double(data(wh_channels_to_keep, :));
            end
        end
    end
end

% function matlab_data_class = get_matlab_data_class(fread_data_class)
%     switch(fread_data_class)
%         case 'float32'
%             matlab_data_class = 'single';
%         case 'float64'
%             matlab_data_class = 'double';
%         otherwise
%             matlab_data_class = fread_data_class;
%     end
% end