%
%   init_exp_defn(experiment, nstream, exp_defn_file, subj_code, montage_file)
%
function init_exp_defn(nstream, exp_defn_file, patient_defn_file, hw_defn_file, subj_code, montage_file)
    global experiment
    
    run(exp_defn_file);

    % This needs to be folded into a check_defn_file function
    if ~isfield(experiment, 'recording')
        error('Experiment definition file needs experiment.recording field')
    end

    if ~isfield(experiment.recording, 'recording_path_base')
        experiment.recording.recording_path_base = sprintf('/mnt/esata/%s', subj_code);  %Default location
    end
    experiment.recording.recording_day = datestr(now, 'yymmdd');
    experiment.recording.recording_path = fullfile(experiment.recording.recording_path_base, experiment.recording.recording_day);
    if ~isdir(experiment.recording.recording_path)
        disp('Creating recording path for recording day');
        disp(sprintf('%s',experiment.recording.recording_path));
        mkdir(experiment.recording.recording_path);
    end

    experiment.initialization.patient_defn_file = patient_defn_file;
    experiment.initialization.hw_defn_file = hw_defn_file;
    experiment.initialization.exp_defn_file = exp_defn_file;
    experiment.initialization.montage_file = montage_file;
   
    for i=1:length(nstream.hardware.nspike.masterdsps)
        nstream_add_masterdsp(nstream.hardware.nspike.masterdsps(i).ip);
        disp(sprintf('Adding MasterDSP: %s',nstream.hardware.nspike.masterdsps(i).ip));
    end
 
    for i=1:16
        nstream_add_auxdsp(nstream.hardware.nspike.auxdsps(i).ip);
    end

    experiment = nview_configure_channels_256(experiment);
    
    for i=1:256
        nstream_set_channel(i, experiment.channels(i).hardware_number);
        nstream_set_filters(i, experiment.channels(i).highpass_cutoff, experiment.channels(i).lowpass_cutoff);
    end
    
    % Read in channel names
    chan_names = textread(montage_file,'%s','delimiter','\n');
    chan_names(cellfun(@isempty, chan_names)) = [];
    num_channels = length(chan_names);
    for i=1:num_channels
        experiment.channels(i).name = chan_names{i};
    end

    % Remove extraneous channels from recording
    experiment.recording.nspike_num_channels_to_write_low = num_channels;


    if isfield(experiment, 'data_glove') && isfield(experiment.data_glove, 'enable') && experiment.data_glove.enable == 1
        nstream_enable_dataglove(nstream.hardware.dataglove(1).serial_port, experiment.data_glove.ticks_between_samples);
    end

    nstream_start_acquire();


    %  Initialize DAC channel and gain
    nstream_set_dac_channel(0, 1);
    nstream_set_dac_channel(1, 1);
    nstream_set_dac_gain(0, 50);
    nstream_set_dac_gain(1, 50);
 
    if length(nstream.hardware.nspike.masterdsps) > 1
        nstream_set_dac_channel(2, 1);
        nstream_set_dac_channel(3, 1);
        nstream_set_dac_gain(2, 50);
        nstream_set_dac_gain(3, 50);
    end
end
