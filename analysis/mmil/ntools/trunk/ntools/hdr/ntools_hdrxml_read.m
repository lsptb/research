%[sfreq, num_channels, num_samples, precision, chan_names, varargout] = ntools_hdrxml_read(file_name, varargin)

function [sfreq, num_channels, num_samples, precision, chan_names, varargout] = ntools_hdrxml_read(hdr_name, varargin)
    [hdr_name d f e] = get_hdr_name(hdr_name);
    doc_node = xmlread(hdr_name);
    
    chan_names = {};
    varargout = {};
    
    root_node = doc_node.getDocumentElement();

    sfreq = read_elt('sfreq', 1);
    num_channels = read_elt('num_channels', 1);
    num_samples = read_elt('num_samples', 1);
    precision = read_elt('data_type', 0);

    
    chan_names_list = root_node.getElementsByTagName('channel_name');

    for i = 0:(chan_names_list.getLength() - 1)
        curr_chan = str2double(chan_names_list.item(i).getAttribute('num'));
        chan_names{curr_chan} = char(chan_names_list.item(i).getTextContent());
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inline functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function value = read_elt(tag_name, is_number)
        value = read_elt_from(root_node, tag_name, is_number);
    end

    function value = read_elt_from(node, tag_name, is_number)
        elt_list = node.getElementsByTagName(tag_name);
        if(elt_list.getLength() == 0)
            error('Missing tag %s.', tag_name);
        elseif(elt_list.getLength() > 1)
            error('Too many instances of tag %s', tag_name);
        else
            elt = elt_list.item(0); % 0 is first elt in Java indexing

            value = char(elt.getTextContent());
            if(is_number)
                value = str2double(value);
            end
        end
    end
end
