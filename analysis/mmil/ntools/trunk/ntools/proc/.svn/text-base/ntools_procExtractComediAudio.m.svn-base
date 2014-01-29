% ntools_procExtractComediAudio(recording_filename_root)
%   processes raw comedi file to extract Audio channel (2).
%   at specified sampling rate in experiment.processing.audio.sample_rate
%
%  Inputs: recording_filename_root - recording root without the suffixes.
%
% E.g.:
%   recording_filename_root = '/mnt/esata/NY140/SP/081030/rec005';
%   ntools_procExtractComediAudio(recording_filename_root);
%
function ntools_procExtractComediAudio(recording_filename_root, fs)

    global experiment

    recording_filename_root = make_filename_root(recording_filename_root); %removes suffix if necess

    comedi_sample_format = experiment.recording.comedi.sample_format;
    audio_sample_format = experiment.processing.comedi.audio.sample_format;
    SIZEOF_SAMPLE = experiment.recording.comedi.sample_size;
    CH = experiment.recording.comedi_num_channels_to_write;
    decimate_factor = round(experiment.recording.comedi.sample_rate./experiment.processing.comedi.audio.sample_rate);
    AUDIOFS = experiment.recording.comedi.sample_rate./decimate_factor;
    AUDIOCH = experiment.recording.comedi.audio_channel_index;  % Channel number for audio data
    
    T = round(1.*1024.^2/(AUDIOFS*SIZEOF_SAMPLE));
    comedi_file = [recording_filename_root '.comedi.dat'];
    if length(dir(comedi_file))
        comedi_fid = fopen(comedi_file, 'r');
        audio_fid = fopen([recording_filename_root '.comedi_audio.dat'], 'w');
        if audio_fid < 0
            ferror(audio_fid);
        end
        status = fseek(comedi_fid,SIZEOF_SAMPLE*(AUDIOCH-1),'bof');
        if status == -1 ferror(comedi_fid); end
        chk=1; N=0;
        audio = zeros(1, AUDIOFS*T);
        while(chk)
            N = N+1;
            disp(['AUDIO: Loop ' num2str(N)]);
            audio = fread(comedi_fid, AUDIOFS*T, comedi_sample_format, (decimate_factor*CH-1)*SIZEOF_SAMPLE);
            audio = audio-32768;
            fwrite(audio_fid, audio, audio_sample_format);
            if (size(audio, 1) < AUDIOFS*T); chk = 0; end
        end
        fclose(audio_fid);
        fclose(comedi_fid);
    else
        disp('No comedi file found.  Skipping.')
    end
end
