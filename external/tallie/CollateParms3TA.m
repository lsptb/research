function X = CollateParms3TA(Parms,fieldname,NumStr)
% For getting cell array data out of individual fields of Parms into
% matrices after using ParamXls1TA.m.
%
% INPUTS:
% Parms     - the output from ParamxXls1TA.m
% fieldname - a fieldname from Parms
% NumStr    - whether the field value should be numeric or string ('num'
%             'str').
%
% OUTPUT:
% X         - If only single numeric values are being collated, this will
%             be a 1D matrix of values for each file. If multiple values
%             exist for each file, or the values are strings, this will be
%             a cell array, with each cell representing the value(s) for
%             that file.
%

if strcmp('str',NumStr);
    % cells of strings that are strings
    X = arrayfun(@(x) x.(fieldname),Parms,'Uni',0);
else
    % cells of strings that are actually numbers
    s2n1 = arrayfun(@(x) x.(fieldname),Parms,'Uni',0);
    if any(cellfun(@iscell,s2n1))
        for i = 1:length(s2n1)
            X{i} = str2double(s2n1{i});
        end
    else
        for i = 1:length(s2n1)
            X(i) = s2n1{i};
        end
    end
end

end

