% byte_size = ntools_get_sample_size(precision)
%
% For a given file's data type, computes the file's internal precision, the
% output precision of fread, and the size of each sample in bytes.

function byte_size = ntools_get_sample_size(precision)
    switch(precision)
        case {'int16', 'uint16'}
            byte_size = 2;
        case {'int32', 'uint32', 'float32', 'single'}
            byte_size = 4;
        case {'int64', 'uint64', 'float64', 'double'}
            byte_size = 8;
        otherwise
            error('Unknown precision: %s', precision);
    end
end
