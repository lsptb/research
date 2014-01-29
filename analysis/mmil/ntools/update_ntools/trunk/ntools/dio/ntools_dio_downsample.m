function ntools_dio_downsample(dio_file, new_sample_rate, orig_sample_rate)
    NSPIKE_SAMPLE_RATE = 30000;
    
    if(~exist('orig_sample_rate', 'var') || isempty(orig_sample_rate))
        orig_sample_rate = NSPIKE_SAMPLE_RATE;
    end
    if(new_sample_rate == orig_sample_rate)
        warning('New sampling rate is the same as the old sampling rate!');
    end
    
    % Convert timestamps
    [timestamps port_codes] = textread(dio_file, '%d %[^\n]');
    timestamps = max(round(timestamps / orig_sample_rate * new_sample_rate), 1);
   
    % Write out file
    [d f e] = fileparts(dio_file);
    [dummy dio_base ee] = fileparts(f); %#ok
    new_dio_file = fullfile(d, sprintf('%s-%dHz%s%s', dio_base, new_sample_rate, ee, e));
    fid = fopen(new_dio_file, 'w');
    
    for i=1:length(timestamps)
        fprintf(fid, '%d %s\n', timestamps(i), port_codes{i});
    end
    
    fclose(fid);
end