% [time_codes, A, B, C, D] = ntools_dio_parse(dio_file, keeping_zeros)
%
% Parses dio.txt file and returns extracted information

function [time_codes, A, B, C, D] = ntools_dio_parse(dio_file, keeping_zeros)
    if(~exist('keeping_zeros', 'var') || isempty(keeping_zeros))
        keeping_zeros = 0;
    end

    dio_data = load(dio_file);
    total_num_lines = size(dio_data, 1);


    if(keeping_zeros)
        wh_lines = 1:total_num_lines;
    else
        wh_lines = (sum(dio_data(:,2:end), 2) ~= 0);
    end

    time_codes = dio_data(wh_lines, 1);
    A = sum(bsxfun(@times, dio_data(wh_lines, 2:9), realpow(2, 7:-1:0)), 2);
    B = sum(bsxfun(@times, dio_data(wh_lines, 10:17), realpow(2, 7:-1:0)), 2);
    C = sum(bsxfun(@times, dio_data(wh_lines, 18:25), realpow(2, 7:-1:0)), 2);
    D = sum(bsxfun(@times, dio_data(wh_lines, 26:33), realpow(2, 7:-1:0)), 2);
end