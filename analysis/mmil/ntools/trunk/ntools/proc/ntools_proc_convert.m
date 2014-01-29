% ntools_proc_convert(in_file, new_precision, new_sfreq)
%   converts .dat file to a different sampling rate and/or data type
%
% Inputs: 
%   in_file - .dat file
%   new_precision - 'int16', 'float32', etc. (defined by fread)
%   new_sfreq - new sampling rate - must be an integer divisor of old sfreq
%
% E.g.:
%   in_file = '/mnt/esata/NY140/SP/081030/rec001.dat';
%   ntools_proc_convert(in_file, 'float32', 1000);

function ntools_proc_convert(in_file, new_precision, new_sfreq)
    [sfreq, num_channels, num_samples, precision, chan_names] = ntools_hdrxml_read(in_file);

    decimate_factor = sfreq/new_sfreq;
    if(mod(decimate_factor, 1) ~= 0)
        error('New sampling rate must be an integer factor of the old rate.');
    end

    block_size = 5;

    lowpass_freq = .4 * new_sfreq;
    filt = gen_filter(sfreq, lowpass_freq);
    Nf = length(filt);
    
    if exist(in_file, 'file')
        in_fid = fopen(in_file, 'r');
        in_size = file_size(in_file);

        [out_fid out_filename] = setup_out_file(in_file);

        out_data = zeros(num_channels, sfreq*block_size+Nf-1, 'single');
        progbarh = waitbar(0, 'Converting...');
        while(~feof(in_fid))
            [data, size_data_read] = fread(in_fid, [num_channels, sfreq*block_size], [precision '=>single']);

            if(size_data_read > 0)
                for j = 1:size(data, 1)
                    out_data(j, 1:size(data, 2)+Nf-1) = conv(data(j,:), filt);
                end
                wh_dec_samples = (Nf/2):decimate_factor:(size(data, 2) + Nf/2 - 1);
                fwrite(out_fid, out_data(:, wh_dec_samples), new_precision);
            end
            waitbar(ftell(in_fid)/in_size, progbarh);            
        end
        close(progbarh);
        fclose(out_fid);
        fclose(in_fid);
        
        ntools_hdrxml_write(out_filename, new_sfreq, num_channels, [], new_precision, chan_names)
    else
        error('File %s not found!', in_file);
    end
end

function filt = gen_filter(sfreq, lowpass_freq)
%     tapers = [0.0025, lowpass_freq];
    tapers = [0.01, lowpass_freq];
    n = tapers(1);
    w = tapers(2);
    p = n*w;
    k = floor(2*p-1);
    tapers = [n, p, k];
    tapers(1) = tapers(1).*sfreq;
    tapers = dpsschk(tapers);
    filt = mtfilt(tapers, sfreq, 0);
    filt = single(filt./sum(filt));
end

function [out_fid out_fullfile] = setup_out_file(in_file)
    out_filename = 0;
    while(out_filename == 0)
        [out_filename, out_pathname] = uiputfile('*.dat', 'Save .dat file as', in_file);
    end
    [d f e] = fileparts(out_filename); %#ok
    if(~strcmp(e, '.dat'))
        out_fullfile = fullfile(out_pathname, [out_filename '.dat']);
    else
        out_fullfile = fullfile(out_pathname, out_filename);
    end
    out_fid = fopen(out_fullfile, 'w');
end