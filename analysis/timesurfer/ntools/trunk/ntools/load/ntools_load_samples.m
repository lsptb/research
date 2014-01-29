% [data, eof] = ntools_load_samples(filename, precision, num_channels, start_samp, stop_samp)
%
% Loads data from a given datafile starting with start_samp (inclusive)
% and ending at end_samp (exclusive) into a matlab matrix in [channel,sample]
% format. If start_samp is empty, it loads the entire file. Samples range
% from 1 to n.
%
% Examples: 
%   data = ntools_load_samples('rec001.low.nspike.dat', 'int16', 256); %loads whole file

function [data, eof] = ntools_load_samples(filename, precision, num_channels, start_samp, stop_samp)
    if(exist('start_samp', 'var') && start_samp < 1)
        error('Start_samp must be >= 1');
    end
    
    [fid, message] = fopen(filename); 
    if(fid == -1), error(message); end

    sample_size = ntools_get_sample_size(precision);
    fread_precision = [precision '=>' precision];

    if(~exist('start_samp', 'var') || isempty(start_samp))
        start_byte = 0;
        data_range = Inf;
    else
        % we want to include start_samp
        start_byte = (start_samp-1) * num_channels * sample_size;
        if start_byte < 0
            start_byte = 0;
        end
        data_range = stop_samp - start_samp + 1;
    end

    status = fseek(fid, start_byte, 'bof');
    if(status == -1), error(ferror(fid)); end

    data = fread(fid, [num_channels, data_range], [num2str(num_channels) '*' fread_precision]);
    eof = feof(fid);
    
    fclose(fid);
end
