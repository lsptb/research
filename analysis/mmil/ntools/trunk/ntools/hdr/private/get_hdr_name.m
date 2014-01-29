function [hdr_name d f e] = get_hdr_name(hdr_name)
    [d f e] = fileparts(hdr_name);
    if(~strcmpi(e, '.hdrxml'))
        if(strcmpi(e, '.dat'))
            hdr_name = fullfile(d, [f '.hdrxml']);
            e = '.hdrxml';
        else
            hdr_name = fullfile(d, [f e '.hdrxml']);
            f = [f e];
            e = '.hdrxml';
        end
    end
end