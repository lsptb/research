% recording = ntools_gen_recording_data(dat_file, experiment)
%
% processes digitalio events for a given recording and generates a recording
% data structure, containing timestamps of known events.
%
% recording_filename_root: root for recording filename, eg rec001
%
% E.g.:
% recording = ntools_gen_recording_data(dat_file, experiment);
%
% examine the experiment.events array for information on events for more
% information

function recording = ntools_gen_recording_data(dat_file, experiment)

    if isempty(dat_file)
        error('Missing dat_file parameter');
    end
    if isempty(experiment)
        error('Missing experiment parameter');
    end

    if(~exist(dat_file, 'file'))
        error('Unable to find .dat file: %s\n', dat_file);
    else
        % try to parse out base from actual filename
        exts = {'.low.nspike.dat', '.high.nspike.dat', '.nspike.dat'}; %order is important
        matches = cellfun(@(x)strfind(dat_file, x), exts, 'UniformOutput', 0);
        base_ends = [matches{~cellfun(@isempty, matches)}] - 1;
        if(isempty(base_ends)), error('Unidentifiable extension in dat file %s\n', dat_file); end
        recording_filename_root = dat_file(1:base_ends(1));
    end  
    
    % copy the pertinent data from the experiment definition file
    recording.settings                = experiment.recording;
    recording.settings.channels       = experiment.channels;
    recording.settings.data_glove     = experiment.data_glove;
    recording.event                   = experiment.event;
    recording.recording_filename_root = recording_filename_root;
    recording.recording_path          = fileparts(experiment.recording.recording_path_base);


    % load digital io data
    dio_data = load([recording.recording_filename_root '.dio.txt']);

    %
    % parse out event codes for DIO and stick them in a structure
    % or "hash table" for quick lookup.
    num_events = size(recording.event,2);
    for i = 1:num_events
        if strcmp(recording.event(i).source, 'dio')
            event_lookup_table.([recording.event(i).source_port num2str(recording.event(i).source_code)]).code = i;
            recording.event(i).timestamps = [];
        end
    end

    % loop over each recorded digital io event, parse the text representation
    % to an integer.  then add them to the recording data structure if they're
    % meaningful to us
    %
    for dio_event_index = 1:size(dio_data,1)
        timestamp = dio_data(dio_event_index,1);
        a_value = parse_8_bits(dio_data(dio_event_index,2:9));
        b_value = parse_8_bits(dio_data(dio_event_index,10:17));
        c_value = parse_8_bits(dio_data(dio_event_index,18:25));
        d_value = parse_8_bits(dio_data(dio_event_index,26:33));

        % now check codes against lookup table for each dio port
        % if we find it, stick it in the experiment defn structure
        event_index = -1;
        if (a_value ~= 0) && (isfield(event_lookup_table, ['A' num2str(a_value)]))
            event_index = event_lookup_table.(['A' num2str(a_value)]).code;
        end
        if (event_index ~= -1)
            recording.event(event_index).timestamps = [recording.event(event_index).timestamps timestamp];
            event_index = -1;
        end

        if (b_value ~= 0) && (isfield(event_lookup_table, ['B' num2str(b_value)]))
            event_index = event_lookup_table.(['B' num2str(b_value)]).code;
        end
        if (event_index ~= -1)
            recording.event(event_index).timestamps = [recording.event(event_index).timestamps timestamp];
            event_index = -1;
        end

        if (c_value ~= 0) && (isfield(event_lookup_table, ['C' num2str(c_value)]))
            event_index = event_lookup_table.(['C' num2str(c_value)]).code;
        end
        if (event_index ~= -1)
            recording.event(event_index).timestamps = [recording.event(event_index).timestamps timestamp];
            event_index = -1;
        end

        if (d_value ~= 0) && (isfield(event_lookup_table, ['D' num2str(d_value)]))
            event_index = event_lookup_table.(['D' num2str(d_value)]).code;
        end
        if (event_index ~= -1)
            recording.event(event_index).timestamps = [recording.event(event_index).timestamps timestamp];
            event_index = -1;
        end
    end
end
