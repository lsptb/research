% ntools_dio_view(dio_file, [displaying_zeros])
%
% Displays information from a DIO file in a human-readable format
%
% Inputs:
%   dio_file - path to the .dio.txt file
%
% Optional inputs:
%   displaying zeros - whether or not to display lines that are all zero (default: off)
%
% Example:
%   dio_file = '~/Desktop/NY000_task1.dio.txt';
%   ntools_dio_view(dio_file);

function ntools_dio_view(dio_file, displaying_zeros)
    NSPIKE_SAMP_RATE = 30000;

    if(~exist('displaying_zeros', 'var') || isempty(displaying_zeros))
        displaying_zeros = 0;
    end

    [time_codes, A, B, C, D] = ntools_dio_parse(dio_file, displaying_zeros);

    if(displaying_zeros)
        fprintf('Displaying all %d entries.\n', size(time_codes, 1));
    else
        fprintf('Displaying %d non-zero entries.\n', size(time_codes, 1));
    end

    for i=1:length(time_codes)
        fprintf('%11.3fs\t\tA: %d\tB: %d\tC: %d\tD: %d\n', time_codes(i)/NSPIKE_SAMP_RATE, A(i), B(i), C(i), D(i));
    end
    
    fprintf('\n\n** Summary **\n');
    fprintf('A: ');
    disp_port_summary(A);
    fprintf('B: ');
    disp_port_summary(B);
    fprintf('C: ');
    disp_port_summary(C);
    fprintf('D: ');
    disp_port_summary(D);
end

function disp_port_summary(data)
    if(any(data ~= 0))
        fprintf('\n');
        [counts values] = hist(data, (0:max(data)));
        wh_not_empty = find(counts > 0);
        for i = 1:length(wh_not_empty)
            fprintf('\t%3d:\t%4d\n', values(wh_not_empty(i)), counts(wh_not_empty(i)));
        end
    else
        fprintf('nothing\n');
    end
end