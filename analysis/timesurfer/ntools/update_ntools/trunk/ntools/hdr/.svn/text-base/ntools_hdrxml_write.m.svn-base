% ntools_hdrxml_write(hdr_name, sfreq, num_channels, num_samples, precision, chan_names)
% 
% Generates a .hdrxml file containing meta-information about a .dat file.
%
% Parameters:
%   hdr_name - new .hdrxml file - will also work with the .dat file
%   sfreq - sampling frequency
%   num_channels - number of channels
%   num_samples - number of channels - if an empty array is passed in, they will be automatically computed
%   precision - data type, e.g., 'int16'
%   chan_names - list of channel names

function ntools_hdrxml_write(hdr_name, sfreq, num_channels, num_samples, precision, chan_names, varargin)
    [hdr_name d f e] = get_hdr_name(hdr_name);
    dat_name = fullfile(d, [f '.dat']);

    verify_inputs();
    
    %Just a placeholder for future keyword args
    %     if(~isempty(varargin))
    %         for i = 1:2:length(varargin)
    %             switch(varargin{i})
    %                 otherwise
    %             end
    %         end
    %     end

    doc_node = com.mathworks.xml.XMLUtils.createDocument('hdr');

    root_node = doc_node.getDocumentElement();
    root_node.setAttribute('fname', [f '.dat']);
    root_node.appendChild(doc_node.createComment('See ntools_hdrxml_write.m and ntools_hdrxml_read.m for info on deciphering this file.'));

    add_elt('sfreq', sfreq);
    add_elt('num_channels', num_channels);
    add_elt('num_samples', num_samples);
    add_elt('data_type', precision);
    
    dur_elt = doc_node.createElement('duration');
    dur_elt.setTextContent(sprintf('%.2f', num_samples/sfreq));
    dur_elt.setAttribute('units', 'seconds');
    root_node.appendChild(dur_elt);

    chan_names_elt = doc_node.createElement('channel_names');
    root_node.appendChild(chan_names_elt);


    for i = 1:num_channels
        if(exist('chan_names', 'var') && ~isempty(chan_names))
            curr_chan_name = chan_names{i};
        else
            curr_chan_name = sprintf('iEEG %03d', i);
        end
        
        elt = doc_node.createElement('channel_name');
        elt.setTextContent(curr_chan_name);
        elt.setAttribute('num', num2str(i));
        chan_names_elt.appendChild(elt);
    end

    % Save the sample XML document.
    xmlwrite(hdr_name, doc_node);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inline functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function verify_inputs()
        if(isempty(num_samples))
            if(num_channels == 0)
                num_samples = 0;
            else
                num_samples = file_size(dat_name) / num_channels / ntools_get_sample_size(precision);
            end
        end

        p = inputParser();
        p.addRequired('sfreq', @validate_sfreq);
        p.addRequired('num_channels', @validate_num_channels);
        p.addRequired('precision', @validate_precision);
        p.addRequired('chan_names', @(chan_names)validate_chan_names(chan_names, num_channels));        
        p.addRequired('num_samples', @(num_samples)validate_num_samples(num_samples, dat_name, num_channels, precision));
        p.parse(sfreq, num_channels, precision, chan_names, num_samples);
    end

    function add_elt(tag_name, elt_value)
        root_node = add_elt_to(root_node, tag_name, elt_value);
    end

    function node = add_elt_to(node, tag_name, elt_value)
        if(isnumeric(elt_value))
            elt_value = num2str(elt_value);
        end

        elt = doc_node.createElement(tag_name);
        elt.setTextContent(elt_value);
        node.appendChild(elt);
    end
end

function is_valid = validate_sfreq(sfreq)
    is_valid = is_integer_value(sfreq) && sfreq <= 30000;
end

function is_valid = validate_num_channels(num_channels)
    is_valid = is_integer_value(num_channels) && num_channels >= 0 && num_channels <= 256;
end

function is_valid = validate_num_samples(num_samples, fname, num_channels, precision)
    byte_size = ntools_get_sample_size(precision);
    is_valid = is_integer_value(num_samples) && (file_size(fname) == num_samples * num_channels * byte_size);
end

function is_valid = validate_precision(precision)
    is_valid = ischar(precision) && any(strcmp(precision, ...
        {'int16', 'int32', 'int64', 'uint16', 'uint32', 'uint64', 'float32', 'float64', 'single', 'double'}));
end

function is_valid = validate_chan_names(chan_names, num_channels)
    is_valid = isempty(chan_names) || (iscellstr(chan_names) && length(chan_names) == num_channels);
end

function is_int = is_integer_value(val)
    is_int = isreal(val) && (mod(val, 1) == 0);
end